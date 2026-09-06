"""Freeze and compare explicit medium preparation on exposed development models.

No numeric fitness values, gene scoring or biological corrections. Exact
stoichiometry/bounds/objective/GPR equality is stronger than matching a few
growth rates: it identifies the same LP and the same reaction-deletion problems.
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import math
from pathlib import Path
import sys
import time

import cobra
import highspy
import swiglpk

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gembench.cards import sha256_of, software_versions
from gembench.fitness_browser import base_medium
from gembench.media import Medium, apply_medium
from gembench.medium_preparation_v2 import CompartmentPolicy, prepare_medium
from gembench.protocols.carbon_fitness_generic import complete_medium_transport
from gembench.study import freeze_study, verify_study
from scripts.diagnose_quinone_producibility import matrix_model, W
from scripts.verify_quinone_repair_solvers import assert_ordinary_lp

PLAN = ROOT / 'data/studies/medium_preparation_v2.json'
INVENTORY = ROOT / 'results/medium_preparation_2026_09_06/evidence/inventory_medium_structure.py'


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def inventory_module():
    spec = importlib.util.spec_from_file_location('medium_source_inventory', INVENTORY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def signature(model):
    """Full algebra and chemical metadata; compartments recorded separately."""
    assert_ordinary_lp(model)
    return {'reactions': {r.id: {'stoichiometry': {m.id: v for m, v in r.metabolites.items()},
            'bounds': list(r.bounds), 'gpr': r.gene_reaction_rule} for r in model.reactions},
            'metabolites': {m.id: {'formula': m.formula, 'charge': m.charge} for m in model.metabolites},
            'objective': {r.id: v for r, v in cobra.util.solver.linear_reaction_coefficients(model).items()},
            'direction': model.objective_direction}


def differences(first, second):
    result = {}
    for kind in first:
        if first[kind] != second[kind]:
            if isinstance(first[kind], dict):
                result[kind] = {key: {'before': first[kind].get(key), 'after': second[kind].get(key)}
                    for key in sorted(set(first[kind]) | set(second[kind])) if first[kind].get(key) != second[kind].get(key)}
            else:
                result[kind] = {'before': first[kind], 'after': second[kind]}
    return result


def digest(value):
    import hashlib
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def medium_for(condition):
    medium = base_medium(condition.media)
    for rid in condition.exchanges:
        medium = medium.with_carbon_source(rid, -10.)
    return medium


def checked_solves(model, recipe):
    assert_ordinary_lp(model)
    matrix = matrix_model(model)
    results = []
    for solver in ('glpk', 'highs'):
        if solver == 'glpk':
            model.solver = 'glpk'
            model.solver.configuration.tolerances.feasibility = recipe['solver_tolerance']
            solution = model.optimize()
            status = model.solver.status.lower()
            values = solution.fluxes.reindex(matrix.rxns).to_numpy() if status == 'optimal' else None
            objective = solution.objective_value if status == 'optimal' else None
        else:
            answer = W.solve_highs(
                matrix, method='simplex', threads=0, feas_tol=recipe['solver_tolerance'],
                opt_tol=recipe['solver_tolerance'], time_limit=60., extra_options={'parallel': 'off'})
            status = answer.status.lower()
            values, objective = (answer.x, answer.objective) if status == 'optimal' else (None, None)
        if status not in ('optimal', 'infeasible'):
            raise RuntimeError(f'Unresolved solver status: {solver}/{status}')
        record = {'solver': solver, 'status': status, 'objective': objective}
        if status == 'optimal':
            fluxes = {str(rid): float(v) for rid, v in zip(matrix.rxns, values)}
            certificate = W.certify(matrix, values)
            if not math.isfinite(objective) or not all(math.isfinite(v) for v in fluxes.values()):
                raise RuntimeError('Nonfinite objective or primal vector')
            if any(not math.isfinite(certificate[key]) or certificate[key] < 0 or certificate[key] > recipe['residual_tolerance']
                   for key in ('max_abs_S_residual', 'max_bound_violation')):
                raise RuntimeError('Numerical certificate exceeds tolerance')
            record.update(full_primal_fluxes=fluxes, primal_certificate=certificate)
        results.append(record)
    if results[0]['status'] != results[1]['status']:
        raise RuntimeError('Solver feasibility classifications disagree')
    if results[0]['status'] == 'optimal' and not math.isclose(results[0]['objective'], results[1]['objective'],
            abs_tol=recipe['objective_tolerance'], rel_tol=recipe['objective_tolerance']):
        raise RuntimeError('Solver objectives disagree')
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    try:
        plan = read(PLAN)
        source = inventory_module()
        evidence = read(ROOT / plan['inventory_record'])
        inventory_by_label = {row['label']: row for row in evidence['models']}
        if set(inventory_by_label) != {row['label'] for row in plan['models']} or len(inventory_by_label) != len(plan['models']):
            raise ValueError('Declared model panel differs from source inventory')
        # Retain the preceding exact source chain as well as this new input set.
        previous = read(ROOT / plan['parent_manifest'])
        old = verify_study(ROOT, previous)
        if not old['valid']:
            raise ValueError(f'Previous inputs changed: {old}')
        paths = {r['path'] for r in previous['files']}
        paths.add(plan['parent_manifest'])
        paths.update(plan['input_paths'])
        paths.update(str(p.relative_to(ROOT)) for p in INVENTORY.parent.rglob('*') if p.is_file())
        paths.update(evidence['input_sha256'])
        for path, expected in evidence['input_sha256'].items():
            if sha256_of(ROOT / path) != expected:
                raise ValueError(f'Source inventory input changed: {path}')
        exposure = read(ROOT / 'data/studies/exposure_registry_v1.json')
        manifest = freeze_study(ROOT, {'schema_version': 1, 'study_id': plan['study_id'],
            'evaluation_role': 'development', 'development_organisms': exposure['development_organisms'],
            'quarantined_organisms': exposure['quarantined_organisms'], 'paths': sorted(paths),
            'recipe': plan, 'runtime': {**software_versions(), 'glpk': swiglpk.glp_version(), 'highs': highspy.Highs().version()}})
        write(args.out / 'manifest.json', manifest)
        total, stats, reuses = 0, [], []
        with gzip.open(args.out / 'cases.jsonl.gz', 'xt') as stream:
            for item in plan['models']:
                label, org = item['label'], item['organism']
                base, provenance = source.build_parent(label)
                policy = CompartmentPolicy(**item['compartments'])
                conditions = source.metadata_conditions(org)
                keys = [{'key': c.key, 'name': c.name, 'media': c.media, 'exchanges': c.exchanges} for c in conditions]
                if keys != inventory_by_label[label]['conditions']:
                    raise ValueError(f'Condition inventory changed: {label}')
                media_names = sorted({c.media for c in conditions})
                completion_media = [base_medium(name) for name in media_names]
                legacy = base.copy()
                added = complete_medium_transport(legacy, media_names, list(policy.completion_exclude))
                # Completion occurs on a closed baseline, before per-condition contexts.
                for r in legacy.exchanges:
                    r.bounds = (0., 1000.)
                label_stats = {'label': label, 'conditions': 0, 'optimal': 0, 'infeasible': 0,
                               'positive_objective': 0, 'added_exchanges': added, 'provenance': provenance}
                for condition in conditions:
                    medium = medium_for(condition)
                    before = legacy.copy()
                    apply_medium(before, medium, close_all=True)
                    prepared = prepare_medium(base, medium, policy=policy,
                                              completion_media=completion_media, missing_policy='report')
                    first, second = signature(before), signature(prepared.model)
                    diff = differences(first, second)
                    if diff:
                        write(args.out / f'{label}_unexpected_difference.json', {'condition': condition.key, 'differences': diff})
                        raise ValueError('Fresh-condition LP or chemical metadata differs; preserve and review before a new plan')
                    # Reusing a correctly prepared model must give exactly this target LP.
                    again = prepare_medium(prepared.model, medium, policy=policy,
                                           completion_media=completion_media, missing_policy='report')
                    if signature(again.model) != second:
                        raise ValueError('Preparation is not idempotent')
                    solves = checked_solves(prepared.model, plan)
                    row = {'model': label, 'condition': condition.key, 'medium_report': prepared.report,
                           'algebra_and_chemistry_sha256': digest(second), 'legacy_signature_sha256': digest(first),
                           'solves': solves}
                    stream.write(json.dumps(row, allow_nan=False) + '\n'); stream.flush()
                    label_stats['conditions'] += 1
                    label_stats[solves[0]['status']] += 1
                    label_stats['positive_objective'] += int(solves[0]['status'] == 'optimal' and solves[0]['objective'] >= plan['growth_threshold'])
                    total += 1
                # Synthetic withdrawal stress test; not an experimental medium or phenotype.
                reference = next((c for c in conditions if 'EX_glc__D_e' in c.exchanges), conditions[0])
                reference_medium = medium_for(reference)
                rich = Medium('completion_withdrawal_start', 'Synthetic bound-reset stress test',
                              {**reference_medium.uptakes, **{rid: -.001 for rid in added}})
                withdrawn = Medium('completion_withdrawal_target', 'Remove all newly completed components',
                                  {rid: lb for rid, lb in reference_medium.uptakes.items() if rid not in added})
                before = legacy.copy(); apply_medium(before, rich); apply_medium(before, withdrawn)
                fresh = legacy.copy(); apply_medium(fresh, withdrawn)
                strict_start = prepare_medium(base, rich, policy=policy, completion_media=completion_media, missing_policy='report')
                strict_reuse = prepare_medium(strict_start.model, withdrawn, policy=policy, completion_media=completion_media, missing_policy='report')
                strict_fresh = prepare_medium(base, withdrawn, policy=policy, completion_media=completion_media, missing_policy='report')
                if signature(strict_reuse.model) != signature(strict_fresh.model):
                    raise ValueError('Nutrient withdrawal depends on prior preparation')
                reuses.append({'label': label, 'reference_condition': reference.key,
                    'legacy_retained_uptakes': {rid: list(before.reactions.get_by_id(rid).bounds) for rid in added if before.reactions.get_by_id(rid).lower_bound < 0},
                    'legacy_reuse_vs_fresh': differences(signature(fresh), signature(before)),
                    'strict_reuse_equals_fresh': True,
                    'strict_reuse_signature_sha256': digest(signature(strict_reuse.model)),
                    'strict_fresh_signature_sha256': digest(signature(strict_fresh.model)),
                    'legacy_reuse_signature_sha256': digest(signature(before)),
                    'legacy_fresh_signature_sha256': digest(signature(fresh)),
                    'strict_reuse_report': strict_reuse.report,
                    'strict_fresh_report': strict_fresh.report,
                    'strict_withdrawn_bounds': {rid: list(strict_reuse.model.reactions.get_by_id(rid).bounds) for rid in added},
                    'scope': 'Synthetic withdrawal; bounds-only comparison, no growth or phenotype claim.'})
                stats.append(label_stats)
                print('Completed:', label, label_stats['conditions'], flush=True)
        verified = verify_study(ROOT, manifest)
        if total != plan['case_count_from_source_inventory']:
            raise ValueError('Case inventory differs from the source declaration')
        if not verified['valid']:
            raise ValueError('Frozen inputs changed during run')
        write(args.out / 'summary.json', {'study': plan, 'study_fingerprint': manifest['content_fingerprint'],
            'models': stats, 'n_conditions': total, 'n_solves': total * 2, 'reuse_stress_tests': reuses,
            'input_verification': verified, 'seconds': time.time() - started,
            'scope': 'Fresh LP/GPR/chemical equality plus solver checks, not a fitness rescore or independent biological validation.'})
    except Exception as error:
        write(args.out / 'failure.json', {'type': type(error).__name__, 'message': str(error), 'seconds': time.time() - started})
        raise


if __name__ == '__main__':
    main()
