"""Verify saved medium-preparation artifacts without the new helper or runner.

Reconstructs historical parents/media, derives an explicit boundary oracle,
checks exact signatures and saved primal vectors. Never optimizes. Infeasible
statuses are retained as solver reports, not independently proved without rays.
Every verification attempt requires a new directory and snapshots this source.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import asdict
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.util
import json
import logging
import math
from pathlib import Path, PurePosixPath
import re
import sys
import time
import warnings

import cobra

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gembench.fitness_browser import base_medium
from gembench.media import Medium, apply_medium
from gembench.protocols.carbon_fitness_generic import complete_medium_transport

SOURCE_LOADER = 'results/medium_preparation_2026_09_06/evidence/inventory_medium_structure.py'
DEFAULT_RUN = ROOT / 'results/medium_preparation_2026_09_06/runs/main_v2'
LABELS = {'Putida', 'Btheta', 'MR1', 'Smeli', 'Putida_iJN1463'}


def require(value, message):
    if not value:
        raise ValueError(message)


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def unique_pairs(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f'Duplicate JSON key: {key}')
        result[key] = value
    return result


def loads(text):
    def invalid(value):
        raise ValueError(f'Nonfinite JSON token: {value}')
    return json.loads(text, object_pairs_hook=unique_pairs, parse_constant=invalid)


def read(path):
    return loads(path.read_text())


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def json_digest(value, *, unicode=False):
    raw = json.dumps(value, sort_keys=True, separators=(',', ':'),
                     ensure_ascii=not unicode, allow_nan=False).encode()
    return hashlib.sha256(raw).hexdigest()


def resolve_input(root, relative):
    require(isinstance(relative, str) and relative and '\\' not in relative and '\0' not in relative,
            'Invalid input path')
    path = PurePosixPath(relative)
    require(not path.is_absolute() and '..' not in path.parts and path.as_posix() == relative,
            'Unsafe/noncanonical input path')
    require(path.parts[0] not in {'.git', '.venv'}, 'Reserved input path')
    resolved = (root / relative).resolve(strict=True)
    require(resolved.is_relative_to(root.resolve()) and resolved.is_file(), 'Input escapes root or is not a file')
    return resolved


def verify_manifest(root, manifest):
    require(manifest.get('schema_version') == 1 and manifest.get('manifest_type') == 'gembench-study-freeze',
            'Unexpected manifest schema')
    plan, records = manifest['plan'], manifest['files']
    require(isinstance(records, list) and records, 'Empty/missing manifest file records')
    paths, resolved = [], set()
    for row in records:
        require(set(row) == {'path', 'size_bytes', 'sha256'}, 'Malformed file record')
        require(type(row['size_bytes']) is int and row['size_bytes'] >= 0, 'Invalid file size')
        require(isinstance(row['sha256'], str) and re.fullmatch('[0-9a-f]{64}', row['sha256']), 'Invalid file digest')
        path = resolve_input(root, row['path'])
        require(path not in resolved, 'Duplicate resolved input')
        require(path.stat().st_size == row['size_bytes'] and sha(path) == row['sha256'], f'Changed input: {row["path"]}')
        paths.append(row['path']); resolved.add(path)
    require(len(set(paths)) == len(paths) and sorted(paths) == plan['paths'], 'Missing/duplicate manifest input')
    fingerprint = json_digest({'schema_version': 1, 'manifest_type': 'gembench-study-freeze',
                               'plan': plan, 'files': sorted(records, key=lambda row: row['path'])}, unicode=True)
    require(fingerprint == manifest['content_fingerprint'], 'Fingerprint differs from plan and file records')
    return {row['path']: row['sha256'] for row in records}


def canonical(model):
    """Independently assemble the published signature schema from the legacy LP."""
    reactions = {}
    for reaction in model.reactions:
        require(reaction.id not in reactions, 'Duplicate reaction identity')
        reactions[reaction.id] = {'stoichiometry': {met.id: value for met, value in reaction.metabolites.items()},
                                  'bounds': [reaction.lower_bound, reaction.upper_bound],
                                  'gpr': reaction.gene_reaction_rule}
    chemicals = {met.id: {'formula': met.formula, 'charge': met.charge} for met in model.metabolites}
    require(len(chemicals) == len(model.metabolites), 'Duplicate metabolite identity')
    objective = {reaction.id: coefficient for reaction, coefficient in
                 cobra.util.solver.linear_reaction_coefficients(model).items()}
    require(model.objective_direction in {'max', 'min'}, 'Unexpected objective direction')
    return {'reactions': reactions, 'metabolites': chemicals, 'objective': objective,
            'direction': model.objective_direction}


def differences(left, right):
    result = {}
    for kind in left:
        if left[kind] != right[kind]:
            if isinstance(left[kind], dict):
                result[kind] = {key: {'before': left[kind].get(key), 'after': right[kind].get(key)}
                               for key in sorted(set(left[kind]) | set(right[kind]))
                               if left[kind].get(key) != right[kind].get(key)}
            else:
                result[kind] = {'before': left[kind], 'after': right[kind]}
    return result


def policy_record(item):
    policy = {'external': None, 'cytoplasm': None, 'external_aliases': [], 'cytoplasm_aliases': [],
              'completion_exclude': ['pnto__R', 'fol', 'hco3'], 'secretion_capacity': 1000.0,
              'exchange_metabolites': {}}
    require(not set(item) - set(policy), 'Unsupported policy field')
    policy.update(item)
    require(policy['external'] != policy['cytoplasm'] and policy['secretion_capacity'] == 1000,
            'Oracle supports the frozen distinct compartments and legacy capacity 1000 only')
    return policy


def boundary_record(model, policy):
    exchanges, boundaries, supply, isolated = [], {}, [], []
    bindings = policy['exchange_metabolites']
    require(isinstance(bindings, dict), 'Identity bindings must be explicit mappings')
    require(all(rid in model.reactions and isinstance(mid, str) and mid.endswith('_e')
                for rid, mid in bindings.items()), 'Unknown exchange identity binding')
    for reaction in model.reactions:
        if reaction.id.startswith('EX_'):
            require(len(reaction.metabolites) == 1, 'Unsupported exchange shape')
            met, coefficient = next(iter(reaction.metabolites.items()))
            require(reaction.id.startswith('EX_') and reaction.id.endswith('_e') and
                    met.id == bindings.get(reaction.id, reaction.id[3:]) and coefficient == -1,
                    'Unsupported exchange identity/sign')
            require(met.compartment == policy['external'], 'Unexpected external compartment')
            exchanges.append(reaction.id)
            if not any(not other.boundary for other in met.reactions):
                isolated.append(reaction.id)
    exchange_set = set(exchanges)
    for reaction in sorted(model.reactions, key=lambda r: r.id):
        if not reaction.boundary:
            continue
        external = reaction.id in exchange_set
        boundaries[reaction.id] = {'bounds': list(reaction.bounds),
            'stoichiometry': {m.id: c for m, c in reaction.metabolites.items()},
            'role': 'exchange' if external else 'internal_boundary'}
        if not external and any(c * bound > 0 for c in reaction.metabolites.values() for bound in reaction.bounds):
            supply.append(reaction.id)
    return {'exchange_ids': sorted(exchanges), 'boundary_reactions': boundaries,
            'internal_supply_capable_boundaries': sorted(set(supply)),
            'exchanges_without_network_connection': sorted(isolated)}


def independent_target(parent, medium, media_names, policy):
    """Oracle using old completion plus explicit, independently applied reset.

    No new helper import. Existing gene/stoichiometric/chemical metadata remains
    the legacy source; only declared labels and environmental bounds are set.
    """
    target = parent.copy()
    renamed = []
    for met in target.metabolites:
        before = met.compartment
        if before in policy['external_aliases']:
            require(met.id.endswith('_e'), 'Invalid source external alias')
            met.compartment = policy['external']
        if before in policy['cytoplasm_aliases']:
            require(met.id.endswith('_c'), 'Invalid source cytoplasm alias')
            met.compartment = policy['cytoplasm']
        if met.compartment != before:
            renamed.append({'metabolite': met.id, 'before': before, 'after': met.compartment})
    requested = sorted(set().union(*(base_medium(name).uptakes for name in media_names)))
    skipped, chemistry_pairs = [], []
    for rid in requested:
        stem = rid[3:-2]
        if rid in parent.reactions:
            skipped.append({'exchange': rid, 'reason': 'exchange_already_exists'})
            if 'MEDt_' + stem in parent.reactions:
                chemistry_pairs.append((stem + '_e', stem + '_c'))
        elif stem in policy['completion_exclude']:
            skipped.append({'exchange': rid, 'reason': 'excluded'})
        elif stem + '_c' not in parent.metabolites:
            skipped.append({'exchange': rid, 'reason': 'no_cytoplasmic_metabolite'})
        else:
            chemistry_pairs.append((stem + '_e', stem + '_c'))
    added = complete_medium_transport(target, media_names, policy['completion_exclude'])
    for met in target.metabolites:
        if met.id.endswith('_e'):
            met.compartment = policy['external']
    for reaction in target.reactions:
        if reaction.id.startswith('EX_'):
            reaction.bounds = (0., policy['secretion_capacity'])
    missing = []
    for rid, bound in medium.uptakes.items():
        if rid in target.reactions:
            target.reactions.get_by_id(rid).lower_bound = bound
        else:
            missing.append(rid)
    chemistry = []
    for external, internal in chemistry_pairs:
        e, c = target.metabolites.get_by_id(external), target.metabolites.get_by_id(internal)
        require(not (e.elements and c.elements) or e.elements == c.elements, 'Known chemistry conflict in oracle')
        require(e.charge is None or c.charge is None or e.charge == c.charge, 'Known charge conflict in oracle')
        chemistry.append({'extracellular': external, 'cytoplasm': internal,
                          'formula_metadata_known': bool(e.elements and c.elements),
                          'charge_metadata_known': e.charge is not None and c.charge is not None})
    expected = {'schema_version': 1, 'policy': policy, 'missing_policy': 'report',
                'medium': asdict(medium), 'completion_media': [asdict(base_medium(name)) for name in media_names],
                'added_exchanges': sorted(added), 'completion_skips': skipped,
                'normalized_compartments': sorted(renamed, key=lambda row: row['metabolite']),
                'completion_chemistry': chemistry, 'missing_medium': sorted(missing), **boundary_record(target, policy)}
    return target, expected


def check_report(actual, expected):
    require(set(actual) == set(expected) | {'scope'}, 'Missing or extra medium-report fields')
    for key, value in expected.items():
        observed = actual[key]
        if key == 'completion_chemistry':
            require(all(set(row) == {'extracellular', 'cytoplasm', 'formula_metadata_known',
                                     'charge_metadata_known', 'note'} for row in observed), 'Malformed chemistry report')
            observed = [{k: v for k, v in row.items() if k != 'note'} for row in observed]
        require(observed == value, f'Medium-report mismatch: {key}')
    require('No biological transport validation' in actual['scope'], 'Boundary report loses interpretation limit')


def primal_certificate(structure, solve, tolerance):
    require(solve['status'] in {'optimal', 'infeasible'}, 'Unresolved saved solver status')
    if solve['status'] == 'infeasible':
        require(solve.get('objective') is None and 'full_primal_fluxes' not in solve and
                'primal_certificate' not in solve, 'Infeasibility represented as a numerical outcome')
        return None
    reactions, fluxes = structure['reactions'], solve['full_primal_fluxes']
    require(set(fluxes) == set(reactions) and all(finite(v) for v in fluxes.values()),
            'Incomplete/nonfinite full primal vector')
    require(finite(solve['objective']), 'Nonfinite objective')
    terms = defaultdict(list)
    violation, bound_count = 0.0, 0
    for rid, row in reactions.items():
        flux = fluxes[rid]
        for mid, coefficient in row['stoichiometry'].items():
            terms[mid].append(coefficient * flux)
        lower, upper = row['bounds']
        current = max(0., lower - flux, flux - upper)
        violation = max(violation, current)
        bound_count += int(current > 1e-6)
    residual_values = [abs(math.fsum(values)) for values in terms.values()]
    residual = max(residual_values, default=0.)
    objective = math.fsum(coefficient * fluxes[rid] for rid, coefficient in structure['objective'].items())
    require(abs(objective - solve['objective']) <= tolerance, 'Objective does not match primal fluxes')
    require(residual <= tolerance and violation <= tolerance, 'Reconstructed primal feasibility failed')
    expected = {'max_abs_S_residual': residual, 'max_bound_violation': violation,
                'n_S_rows_over_1e-6': sum(r > 1e-6 for r in residual_values), 'n_bounds_over_1e-6': bound_count}
    reported = solve['primal_certificate']
    require(set(reported) == set(expected), 'Unexpected certificate fields')
    for key, value in expected.items():
        require(finite(reported[key]) and reported[key] >= 0, 'Invalid certificate value')
        if key.startswith('n_'):
            require(reported[key] == value, 'Wrong residual violation count')
        else:
            # fsum and the runner's sparse dot product need not round identically.
            require(abs(reported[key] - value) <= 1e-9, 'Stored certificate disagrees with independent summation')
    return {**expected, 'objective_recomputed': objective}


def medium_for(condition):
    medium = base_medium(condition.media)
    for exchange in condition.exchanges:
        medium = medium.with_carbon_source(exchange, -10.)
    return medium


def load_source():
    spec = importlib.util.spec_from_file_location('independent_medium_parent_loader', ROOT / SOURCE_LOADER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def retained_attempt(plan, frozen):
    """The v1 abort is frozen evidence of a partial attempt, never a success."""
    prefix = 'results/medium_preparation_2026_09_06/runs/main/'
    paths = [prefix + name for name in ('manifest.json', 'failure.json', 'cases.jsonl.gz')]
    require(all(path in frozen and path in plan['input_paths'] for path in paths),
            'The retained v1 failure artifacts were not frozen')
    manifest, failure = read(ROOT / paths[0]), read(ROOT / paths[1])
    verify_manifest(ROOT, manifest)
    require(manifest['plan']['recipe']['study_id'] == 'explicit-medium-preparation-v1',
            'Retained attempt has a different declaration')
    require(failure['type'] == 'ValueError' and
            failure['message'] == 'Noncanonical exchange identity/orientation: EX_AEP_e',
            'Retained failure is not the documented source-identity refusal')
    require(not (ROOT / prefix / 'summary.json').exists(), 'Incomplete attempt unexpectedly has a completed summary')
    with gzip.open(ROOT / paths[2], 'rt') as handle:
        rows = [loads(line) for line in handle if line.strip()]
    indexed = {(row['model'], row['condition']): row for row in rows}
    require(len(indexed) == len(rows) == 113 and
            {row['model'] for row in rows} == LABELS - {'Putida_iJN1463'},
            'Retained incomplete attempt inventory changed')
    return indexed, {'scope': 'Incomplete v1 attempt; all four drafts completed before the curated-model identity guard refused EX_AEP_e.',
                     'cases': len(rows), 'finite_primal_vectors_verified': 0,
                     'infeasible_solver_reports_retained': 0,
                     'artifact_sha256': {path: frozen[path] for path in paths},
                     'study_fingerprint': manifest['content_fingerprint']}


def strict_prefix(cases, inventory):
    """Require ordered, unique coverage of an incomplete declared prefix."""
    actual = [(row['model'], row['condition']) for row in cases]
    require(0 < len(actual) < len(inventory), 'A partial attempt must be nonempty and incomplete')
    require(len(set(actual)) == len(actual) and actual == inventory[:len(actual)],
            'Saved cases are not the unique, ordered declaration prefix')
    return set(actual)


def verify_partial(run):
    """Check the two preserved 113-case attempts without inventing a summary."""
    names = ('manifest.json', 'failure.json', 'cases.jsonl.gz')
    paths = [run / name for name in names]
    require(all(path.is_file() for path in paths) and not (run / 'summary.json').exists(),
            'Partial verification requires failure and saved cases, with no completion summary')
    initial = {str(path.relative_to(ROOT)): sha(path) for path in paths}
    manifest, failure = read(paths[0]), read(paths[1])
    frozen = verify_manifest(ROOT, manifest)
    plan = manifest['plan']['recipe']
    version = plan['study_id']
    failures = {'explicit-medium-preparation-v1': ('ValueError', 'Noncanonical exchange identity/orientation: EX_AEP_e'),
                'explicit-medium-preparation-v2': ('RuntimeError', 'Numerical certificate exceeds tolerance')}
    require(version in failures and (failure['type'], failure['message']) == failures[version],
            'Partial attempt does not have its documented failure')
    require(len(plan['models']) == len(LABELS) and {m['label'] for m in plan['models']} == LABELS,
            'Changed or duplicate model panel')
    source = load_source()
    conditions = {m['label']: source.metadata_conditions(m['organism']) for m in plan['models']}
    inventory = [(m['label'], condition.key) for m in plan['models'] for condition in conditions[m['label']]]
    require(len(inventory) == plan['case_count_from_source_inventory'] == 156, 'Declared case inventory differs')
    with gzip.open(paths[2], 'rt') as handle:
        cases = [loads(line) for line in handle if line.strip()]
    observed = strict_prefix(cases, inventory)
    require(len(cases) == 113 and {label for label, _ in observed} == LABELS - {'Putida_iJN1463'},
            'This audit covers the documented four-draft, 113-case partial attempts')
    indexed = {(row['model'], row['condition']): row for row in cases}
    finite_solves, infeasible_solves, value_count, row_count, boundary_count = 0, 0, 0, 0, 0
    maxima = {'S_residual': 0., 'bound_violation': 0., 'solver_objective_difference': 0.}
    models = []
    for item in plan['models']:
        label = item['label']
        selected = [condition for condition in conditions[label] if (label, condition.key) in observed]
        if not selected:
            continue
        base, provenance = source.build_parent(label)
        for path, digest in provenance['input_sha256'].items():
            require(frozen.get(path) == digest, f'Unfrozen parent reconstruction dependency: {path}')
        media_names = sorted({condition.media for condition in conditions[label]})
        policy = policy_record(item['compartments'])
        require(policy['exchange_metabolites'] == {}, 'Draft panel unexpectedly has exchange identity overrides')
        legacy = base.copy()
        added = complete_medium_transport(legacy, media_names, policy['completion_exclude'])
        for reaction in legacy.exchanges:
            reaction.bounds = (0., 1000.)
        counts = Counter(conditions=0, optimal=0, infeasible=0, finite_primal_vectors=0,
                         primal_values=0, mass_balance_rows=0, boundary_records=0)
        sizes = set()
        for condition in selected:
            case = indexed[(label, condition.key)]
            medium = medium_for(condition)
            before = legacy.copy()
            apply_medium(before, medium, close_all=True)
            signature = canonical(before)
            digest = json_digest(signature)
            require(case['legacy_signature_sha256'] == case['algebra_and_chemistry_sha256'] == digest,
                    f'Partial saved signature differs from independent legacy reconstruction: {label}/{condition.key}')
            target, report = independent_target(base, medium, media_names, policy)
            require(canonical(target) == signature, 'Independent explicit preparation changes the saved legacy LP')
            if version.endswith('-v1'):
                report['policy'] = {k: value for k, value in policy.items() if k != 'exchange_metabolites'}
            check_report(case['medium_report'], report)
            nr, nm = len(signature['reactions']), len(signature['metabolites'])
            sizes.add((nr, nm))
            counts['boundary_records'] += len(report['boundary_reactions'])
            boundary_count += len(report['boundary_reactions'])
            solves = case['solves']
            require(len(solves) == 2 and {s['solver'] for s in solves} == {'glpk', 'highs'},
                    'Partial case has missing or duplicate solvers')
            require(solves[0]['status'] == solves[1]['status'], 'Partial solver feasibility reports disagree')
            for solve in solves:
                certificate = primal_certificate(signature, solve, plan['residual_tolerance'])
                if certificate is None:
                    infeasible_solves += 1
                else:
                    finite_solves += 1
                    value_count += nr
                    row_count += nm
                    counts.update(finite_primal_vectors=1, primal_values=nr, mass_balance_rows=nm)
                    maxima['S_residual'] = max(maxima['S_residual'], certificate['max_abs_S_residual'])
                    maxima['bound_violation'] = max(maxima['bound_violation'], certificate['max_bound_violation'])
            counts['conditions'] += 1
            counts[solves[0]['status']] += 1
            if solves[0]['status'] == 'optimal':
                a, b = (solve['objective'] for solve in solves)
                require(math.isclose(a, b, rel_tol=plan['objective_tolerance'], abs_tol=plan['objective_tolerance']),
                        'Partial solver objectives disagree')
                maxima['solver_objective_difference'] = max(maxima['solver_objective_difference'], abs(a - b))
        models.append({'label': label, **counts, 'added_exchanges': added,
                       'prepared_sizes': [{'reactions': nr, 'metabolites': nm} for nr, nm in sorted(sizes)],
                       'source_provenance': provenance})
    require(frozen == verify_manifest(ROOT, manifest), 'Inputs changed during partial verification')
    require(initial == {str(path.relative_to(ROOT)): sha(path) for path in paths}, 'Partial artifacts changed during verification')
    return {'status': 'passed_partial', 'coverage': 'Only saved cases verified; the 156-case study remains incomplete.',
            'cases_verified': len(cases), 'declared_cases': len(inventory), 'unsaved_cases': len(inventory) - len(cases),
            'finite_primal_vectors_verified': finite_solves, 'infeasible_solver_reports_retained': infeasible_solves,
            'reaction_flux_values_verified': value_count, 'mass_balance_rows_verified': row_count,
            'boundary_reports_verified': len(cases), 'boundary_records_verified': boundary_count,
            'strict_prefix_inventory_verified': True, 'failure_preserved': failure, 'models': models, 'maxima': maxima,
            'frozen_input_files_verified': len(frozen), 'study_fingerprint': manifest['content_fingerprint'],
            'frozen_input_sha256': frozen, 'direct_artifact_sha256': initial, 'source_sha256': sha(Path(__file__).resolve()),
            'method': 'Independent legacy parent/media reconstruction, explicit reset oracle, complete boundary reports, exact canonical signatures and full-primal math.fsum residuals; no new helper/runner imported and no optimizer called.',
            'limitations': ['The failing curated-model case has no saved vector in these attempts and is not verified here.',
                'Synthetic reuse reports were not saved before the failure; this partial audit cannot verify them.',
                'Primal feasibility is not a dual optimality certificate; an infeasible solver status without a Farkas ray is not independently proved.',
                'No mutant fitness, measured nutrient availability or biological transport validation is evaluated.']}


def verify(run):
    paths = [run / name for name in ('manifest.json', 'summary.json', 'cases.jsonl.gz')]
    require(all(path.is_file() for path in paths), 'Completed artifacts required before verification')
    initial = {str(path.relative_to(ROOT)): sha(path) for path in paths}
    manifest, summary = read(paths[0]), read(paths[1])
    frozen = verify_manifest(ROOT, manifest)
    plan = summary['study']
    require(plan == manifest['plan']['recipe'], 'Summary recipe differs from frozen declaration')
    require(summary['study_fingerprint'] == manifest['content_fingerprint'], 'Summary fingerprint mismatch')
    require(summary['input_verification']['valid'], 'Main run did not verify its inputs')
    require(plan['study_id'] == 'explicit-medium-preparation-v2', 'This checker targets the v2 declaration')
    require({row['label'] for row in plan['models']} == LABELS and len(plan['models']) == len(LABELS), 'Changed/duplicate model panel')
    source = load_source()
    partial, partial_result = retained_attempt(plan, frozen)
    partial_checked = set()
    with gzip.open(paths[2], 'rt') as handle:
        cases = [loads(line) for line in handle if line.strip()]
    indexed = {(row['model'], row['condition']): row for row in cases}
    require(len(indexed) == len(cases), 'Duplicate case identities')
    model_stats = {row['label']: row for row in summary['models']}
    reuses = {row['label']: row for row in summary['reuse_stress_tests']}
    require(len(model_stats) == len(summary['models']) == len(LABELS) and set(model_stats) == LABELS, 'Missing model summary')
    require(len(reuses) == len(summary['reuse_stress_tests']) == len(LABELS) and set(reuses) == LABELS, 'Missing reuse summary')
    inventory, outcomes = set(), []
    maxima = {'S_residual': 0., 'bound_violation': 0., 'solver_objective_difference': 0.}
    finite_solves, infeasible_solves, reports = 0, 0, 0
    for item in plan['models']:
        label = item['label']
        base, provenance = source.build_parent(label)
        require(provenance == model_stats[label]['provenance'], 'Parent reconstruction provenance differs')
        for path, digest in provenance['input_sha256'].items():
            require(frozen.get(path) == digest, f'Reconstruction dependency not frozen: {path}')
        conditions = source.metadata_conditions(item['organism'])
        media_names = sorted({c.media for c in conditions})
        policy = policy_record(item['compartments'])
        require(policy['exchange_metabolites'] == ({'EX_AEP_e': '2ameph_e'} if label == 'Putida_iJN1463' else {}),
                'Exchange binding differs from the source-audited v2 exception')
        legacy = base.copy()
        added = complete_medium_transport(legacy, media_names, policy['completion_exclude'])
        for reaction in legacy.exchanges:
            reaction.bounds = (0., 1000.)
        counts = Counter(conditions=len(conditions))
        boundary_counts = Counter()
        for condition in conditions:
            key = (label, condition.key)
            inventory.add(key)
            require(key in indexed, f'Missing declared case: {key}')
            case = indexed[key]
            medium = medium_for(condition)
            before = legacy.copy()
            apply_medium(before, medium, close_all=True)
            signature = canonical(before)
            digest = json_digest(signature)
            require(case['legacy_signature_sha256'] == case['algebra_and_chemistry_sha256'] == digest,
                    f'Saved fresh signature differs from reconstructed legacy: {key}')
            target, expected_report = independent_target(base, medium, media_names, policy)
            require(canonical(target) == signature, f'Independent boundary reset changes legacy fresh LP: {key}')
            check_report(case['medium_report'], expected_report)
            if key in partial:
                prior = partial[key]
                require(prior['legacy_signature_sha256'] == prior['algebra_and_chemistry_sha256'] == digest,
                        f'Retained v1 LP differs: {key}')
                prior_report = {**expected_report, 'policy': {k: v for k, v in policy.items()
                                                             if k != 'exchange_metabolites'}}
                check_report(prior['medium_report'], prior_report)
                require(len(prior['solves']) == 2 and {s['solver'] for s in prior['solves']} == {'glpk', 'highs'},
                        'Retained attempt has missing/duplicate solvers')
                require(prior['solves'][0]['status'] == prior['solves'][1]['status'],
                        'Retained solver feasibility classifications disagree')
                for solve in prior['solves']:
                    certificate = primal_certificate(signature, solve, plan['residual_tolerance'])
                    count = ('infeasible_solver_reports_retained' if certificate is None
                             else 'finite_primal_vectors_verified')
                    partial_result[count] += 1
                partial_checked.add(key)
            reports += 1
            boundary_counts.update({'boundary_records_checked': len(expected_report['boundary_reactions']),
                                    'missing_components_reported': len(expected_report['missing_medium'])})
            solves = case['solves']
            require(len(solves) == 2 and {s['solver'] for s in solves} == {'glpk', 'highs'}, 'Missing/duplicate solver')
            require(solves[0]['status'] == solves[1]['status'], 'Saved feasibility classifications differ')
            for solve in solves:
                certificate = primal_certificate(signature, solve, plan['residual_tolerance'])
                if certificate is None:
                    infeasible_solves += 1
                else:
                    finite_solves += 1
                    maxima['S_residual'] = max(maxima['S_residual'], certificate['max_abs_S_residual'])
                    maxima['bound_violation'] = max(maxima['bound_violation'], certificate['max_bound_violation'])
            status = solves[0]['status']
            counts[status] += 1
            counts['positive_objective'] += int(status == 'optimal' and solves[0]['objective'] >= plan['growth_threshold'])
            if status == 'optimal':
                difference = abs(solves[0]['objective'] - solves[1]['objective'])
                require(math.isclose(solves[0]['objective'], solves[1]['objective'], rel_tol=plan['objective_tolerance'],
                                     abs_tol=plan['objective_tolerance']), 'Saved solver objectives differ')
                maxima['solver_objective_difference'] = max(maxima['solver_objective_difference'], difference)
        for field in ('conditions', 'optimal', 'infeasible', 'positive_objective'):
            require(model_stats[label][field] == counts[field], f'Model count mismatch: {label}/{field}')
        require(model_stats[label]['added_exchanges'] == added, 'Completion inventory mismatch')
        # Independent reuse oracle: old inference-driven reset versus explicit EX reset.
        reference = next((c for c in conditions if 'EX_glc__D_e' in c.exchanges), conditions[0])
        medium = medium_for(reference)
        rich = Medium('completion_withdrawal_start', 'Synthetic bound-reset stress test',
                      {**medium.uptakes, **{rid: -.001 for rid in added}})
        withdrawn = Medium('completion_withdrawal_target', 'Remove all newly completed components',
                           {rid: lower for rid, lower in medium.uptakes.items() if rid not in added})
        old_reuse, old_fresh = legacy.copy(), legacy.copy()
        apply_medium(old_reuse, rich); apply_medium(old_reuse, withdrawn); apply_medium(old_fresh, withdrawn)
        start, _ = independent_target(base, rich, media_names, policy)
        strict_reuse, reuse_report = independent_target(start, withdrawn, media_names, policy)
        strict_fresh, fresh_report = independent_target(base, withdrawn, media_names, policy)
        reuse = reuses[label]
        require(reuse['reference_condition'] == reference.key, 'Changed reuse reference condition')
        require(reuse['scope'] == 'Synthetic withdrawal; bounds-only comparison, no growth or phenotype claim.', 'Reuse scope changed')
        for name, model in [('legacy_reuse', old_reuse), ('legacy_fresh', old_fresh),
                            ('strict_reuse', strict_reuse), ('strict_fresh', strict_fresh)]:
            require(reuse[name + '_signature_sha256'] == json_digest(canonical(model)), f'Reuse signature mismatch: {label}/{name}')
        require(canonical(strict_reuse) == canonical(strict_fresh) and reuse['strict_reuse_equals_fresh'] is True,
                'Strict withdrawal state depends on history')
        check_report(reuse['strict_reuse_report'], reuse_report)
        check_report(reuse['strict_fresh_report'], fresh_report)
        reports += 2
        retained = {rid: list(old_reuse.reactions.get_by_id(rid).bounds) for rid in added
                    if old_reuse.reactions.get_by_id(rid).lower_bound < 0}
        closed = {rid: list(strict_reuse.reactions.get_by_id(rid).bounds) for rid in added}
        require(reuse['legacy_retained_uptakes'] == retained, 'Wrong legacy retained uptake report')
        require(reuse['legacy_reuse_vs_fresh'] == differences(canonical(old_fresh), canonical(old_reuse)), 'Wrong legacy reuse difference')
        require(reuse['strict_withdrawn_bounds'] == closed and all(value == [0., 1000.] for value in closed.values()),
                'Withdrawn completed component is not closed')
        outcomes.append({'label': label, **counts, **boundary_counts, 'legacy_retained_uptake_ids': sorted(retained),
                         'strict_closed_completed_exchange_ids': sorted(closed), 'fresh_signatures_verified': len(conditions),
                         'reuse_full_signatures_and_reports_verified': True})
    require(set(indexed) == inventory, 'Unexpected extra case')
    require(set(partial) == partial_checked == {key for key in inventory if key[0] != 'Putida_iJN1463'},
            'Retained attempt does not cover exactly the four draft panels')
    require(len(cases) == summary['n_conditions'] == plan['case_count_from_source_inventory'], 'Case count mismatch')
    require(2 * len(cases) == summary['n_solves'] == finite_solves + infeasible_solves, 'Solve count mismatch')
    require(frozen == verify_manifest(ROOT, manifest), 'Inputs changed during checker')
    require(initial == {str(path.relative_to(ROOT)): sha(path) for path in paths}, 'Outputs changed during checker')
    return {'status': 'passed', 'method': 'Legacy parent/media reconstruction, independent explicit boundary oracle, canonical JSON hashes and math.fsum primal checks; no new helper/runner imported and no optimizer called.',
            'models': outcomes, 'cases_verified': len(cases), 'boundary_reports_verified': reports,
            'finite_primal_vectors_verified': finite_solves, 'infeasible_solver_reports_retained': infeasible_solves,
            'maxima': maxima, 'frozen_input_files_verified': len(frozen), 'study_fingerprint': manifest['content_fingerprint'],
            'frozen_input_sha256': frozen, 'direct_artifact_sha256': initial,
            'source_sha256': sha(Path(__file__).resolve()),
            'retained_incomplete_attempt': partial_result,
            'limitations': ['No independent proof of infeasibility without a Farkas certificate; paired solver infeasible reports remain distinct from zero objective.',
                'Primal feasibility and agreeing reported optima are not dual optimality certificates.',
                'The saved new signature/report is checked against reconstructed algebra; the new implementation is deliberately not rerun by this checker.',
                'Internal supply boundaries remain inherited and disclosed; this audit does not establish their biological validity.',
                'Synthetic withdrawal concerns state reset only, not measured nutrient availability, mutant fitness or a biological correction.']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, default=DEFAULT_RUN)
    parser.add_argument('--out', type=Path, help='Fresh attempt directory; default <run>/independent_verification')
    parser.add_argument('--partial', action='store_true', help='Verify the documented 113-case prefix of a failed v1/v2 attempt')
    args = parser.parse_args()
    args.run = args.run.resolve()
    require(args.run.is_relative_to(ROOT), 'Run artifacts must be within the repository root')
    out = args.out or args.run / 'independent_verification'
    out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    source = Path(__file__).resolve()
    with (out / 'checker_source.py').open('xb') as handle:
        handle.write(source.read_bytes())
    with (out / 'attempt.json').open('x') as handle:
        json.dump({'started_utc': datetime.now(timezone.utc).isoformat(), 'source': str(source.relative_to(ROOT)),
                   'source_sha256': sha(source), 'scope': 'Additive independent checker, separate from frozen experiment inputs.'}, handle, indent=2)
        handle.write('\n')
    logging.getLogger('cobra').setLevel(logging.ERROR)
    warnings.filterwarnings('ignore', category=FutureWarning, module=r'cobra\.medium\.boundary_types')
    def forbidden(*args, **kwargs):
        raise AssertionError('Independent artifact checker must not optimize')
    cobra.Model.optimize = forbidden
    cobra.Model.slim_optimize = forbidden
    try:
        result = verify_partial(args.run) if args.partial else verify(args.run)
        result['elapsed_seconds'] = time.time() - started
        with (out / 'result.json').open('x') as handle:
            json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write('\n')
        print(json.dumps({key: result[key] for key in ('status', 'cases_verified', 'finite_primal_vectors_verified',
                                                       'infeasible_solver_reports_retained', 'maxima')}))
    except Exception as error:
        with (out / 'failure.json').open('x') as handle:
            json.dump({'type': type(error).__name__, 'message': str(error), 'elapsed_seconds': time.time() - started,
                       'scope': 'Checker attempt preserved; primary artifacts remain unchanged.'}, handle, indent=2)
            handle.write('\n')
        raise


if __name__ == '__main__':
    main()
