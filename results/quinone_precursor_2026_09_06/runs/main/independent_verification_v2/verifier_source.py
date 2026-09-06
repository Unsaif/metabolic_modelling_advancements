"""Independently audit saved precursor-supply artifacts without solving an LP.

Reconstruct media, gene restrictions and full primal vectors for both solvers.
Check declared inventory, exact row sums and necessary supply bounds. Solver
agreement plus primal feasibility is not a dual/optimality proof or evidence of
biological sharing. Only experiment metadata and fitness column names are read.
"""
from __future__ import annotations

import argparse
import ast
import csv
from datetime import datetime, timezone
from fractions import Fraction
import gzip
import hashlib
from importlib.metadata import version
import json
import math
from pathlib import Path
import platform
import shutil
import sys

import cobra
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EX = 'EX_4hbz_e'
COUM = 'EX_T4hcinnm_e'
POOL = ('4hbz_c', '4hbz_p', '4hbz_e', '4hbzcoa_c', '3ophb_c', '2oph_c',
        '2ohph_c', '2omph_c', '2ombzl_c', '2ommbl_c', '2omhmbl_c', 'q8_c',
        'q8h2_c', 'T4hcinnm_c', 'T4hcinnm_p', 'T4hcinnm_e', 'coucoa_c', '4hbald_c')
METALS = {'fe2', 'fe3', 'mn2', 'zn2', 'cu2', 'cobalt2', 'mobd', 'ni2', 'sel', 'slnt', 'tungs'}
EXCLUDED_COMPLETION = {'pnto__R', 'fol', 'hco3'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def near(a, b, tolerance=1e-8):
    return finite(a) and finite(b) and math.isclose(a, b, abs_tol=tolerance, rel_tol=tolerance)


def pairs_unique(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f'Duplicate JSON field: {key}')
        result[key] = value
    return result


def read_json(path):
    return json.loads(path.read_text(), object_pairs_hook=pairs_unique,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError(f'Nonfinite JSON: {value}')))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write('\n')


def table(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def verify_manifest(root, manifest):
    records = manifest['files']
    paths = [record['path'] for record in records]
    require(len(set(paths)) == len(paths) and sorted(paths) == sorted(manifest['plan']['paths']), 'Manifest path inventory differs')
    payload = {key: manifest[key] for key in ('schema_version', 'manifest_type', 'plan')}
    payload['files'] = sorted(records, key=lambda record: record['path'])
    fingerprint = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False).encode()).hexdigest()
    require(fingerprint == manifest['content_fingerprint'], 'Manifest fingerprint differs')
    for record in records:
        relative = Path(record['path'])
        path = (root / relative).resolve()
        require(not relative.is_absolute() and '..' not in relative.parts and path.is_relative_to(root), 'Unsafe input path')
        require(path.is_file() and path.stat().st_size == record['size_bytes'] and sha(path) == record['sha256'], f'Frozen input changed: {relative}')
    return {'valid': True, 'checked_files': len(records), 'content_fingerprint': fingerprint}


def gpr_active(rule, knocked):
    if not rule:
        return True
    def evaluate(node):
        if isinstance(node, ast.Name):
            return node.id not in knocked
        if isinstance(node, ast.BoolOp):
            values = [evaluate(value) for value in node.values]
            if isinstance(node.op, ast.And):
                return all(values)
            if isinstance(node.op, ast.Or):
                return any(values)
        raise ValueError('Unsupported gene-rule expression')
    return evaluate(ast.parse(rule, mode='eval').body)


def medium_bounds(record):
    components = record['bigg_ids'].split(';')
    if record['aerobic'].lower() in ('yes', 'true', '1') and 'o2' not in components:
        components.append('o2')
    trace = set(record['trace_components'].split(';')) - {''}
    return {f'EX_{met}_e': -.001 if met in trace else -.1 if met in METALS else -1000.
            for met in components if met}


def metadata_conditions(root):
    media = {row['media']: row for row in table(root / 'data/reference/fitness_browser_media_bigg.tsv')}
    carbon = {row['condition']: row for row in table(root / 'data/reference/fitness_browser_carbon_sources_bigg.tsv')}
    with (root / 'data/fitness_browser/Putida/fit_logratios.tsv').open(newline='') as stream:
        header = next(csv.reader(stream, delimiter='\t'))
    columns = {column.split(' ')[0] for column in header if column not in {'orgId', 'locusId', 'sysName', 'geneName', 'desc'}}
    groups = {}
    for row in table(root / 'data/fitness_browser/Putida/experiments.tsv'):
        if row['expGroup'] != 'carbon source' or row['condition_2'] not in ('', 'Dimethyl Sulfoxide') or row['media'] not in media or row['expName'] not in columns:
            continue
        key = row['condition_1'] + ' | ' + row['media']
        if key not in groups:
            ids = [met for met in carbon.get(row['condition_1'], {}).get('bigg_ids', '').split(';') if met]
            groups[key] = {'key': key, 'name': row['condition_1'], 'media': row['media'],
                           'exchanges': [f'EX_{met}_e' for met in ids]}
    return [condition for condition in groups.values() if condition['exchanges']], media


def network_from_model(model):
    objective = {r.id: float(value) for r, value in cobra.util.solver.linear_reaction_coefficients(model).items()}
    require(objective == {'Growth': 1.} and model.objective_direction == 'max', 'Unexpected base objective')
    records = {r.id: {'stoich': {m.id: float(value) for m, value in r.metabolites.items()},
                      'bounds': tuple(r.bounds), 'gpr': r.gene_reaction_rule} for r in model.reactions}
    return records, {met.id: met.compartment for met in model.metabolites}, {r.id for r in model.exchanges}


def prepare_network(base, metabolites, exchanges, conditions, media, condition, recipe):
    """Build a numeric LP description, without modifying a COBRA model."""
    records = dict(base)
    metabolites, exchanges = dict(metabolites), set(exchanges)
    external = {metabolites[mid] for rid in exchanges for mid in records[rid]['stoich']}
    require(len(external) <= 1, 'Ambiguous parent external compartment')
    external_compartment = next(iter(external), 'e')
    added = []
    for name in sorted({row['media'] for row in conditions}):
        for ex in medium_bounds(media[name]):
            met = ex[3:-2]
            if ex in records or met + '_c' not in metabolites or met in EXCLUDED_COMPLETION:
                continue
            require('MEDt_' + met not in records, 'Unexpected medium-carrier collision')
            metabolites.setdefault(met + '_e', 'e')
            records[ex] = {'stoich': {met + '_e': -1.}, 'bounds': (0., 1000.), 'gpr': ''}
            records['MEDt_' + met] = {'stoich': {met + '_e': -1., met + '_c': 1.}, 'bounds': (0., 1000.), 'gpr': ''}
            # The inherited completion creates new metabolites in "e" while
            # the parent uses "C_e". COBRA excludes those new boundaries from
            # model.exchanges/context, but they remain present in the full LP.
            if metabolites[met + '_e'] == external_compartment:
                exchanges.add(ex)
            added.append(ex)
    bounds = {rid: (0., 1000.) if rid in exchanges else record['bounds'] for rid, record in records.items()}
    missing = []
    for ex, lower in medium_bounds(media[condition['media']]).items():
        if ex in bounds:
            bounds[ex] = (lower, bounds[ex][1])
        else:
            missing.append(ex)
    absent_carbon = []
    for ex in condition['exchanges']:
        if ex in bounds:
            bounds[ex] = (recipe['protocol']['carbon_uptake'], bounds[ex][1])
        else:
            absent_carbon.append(ex)
    context = {'condition': condition['key'], 'carbon_exchanges': condition['exchanges'],
               'medium_completion_added': sorted(added), 'missing_medium': missing, 'missing_carbon': absent_carbon,
               'exchange_bounds': {rid: list(bounds[rid]) for rid in sorted(exchanges)},
               'growth_bounds': list(bounds['Growth']), 'atpm_bounds': list(bounds['ATPM'])}
    return records, bounds, context


def restricted_bounds(records, bounds, specification, aliases):
    result = dict(bounds)
    knocked = {aliases[locus] for locus in specification.get('loci', [])}
    for rid, record in records.items():
        if not gpr_active(record['gpr'], knocked):
            result[rid] = (0., 0.)
    for rid in specification.get('close_reactions', []):
        require(rid in result, 'Missing declared reaction closure')
        result[rid] = (0., 0.)
    for field, rid in (('exchange_bounds', EX), ('growth_bounds', 'Growth')):
        if field in specification:
            value = specification[field]
            require(len(value) == 2 and all(finite(v) for v in value) and value[0] <= value[1], 'Invalid specified bounds')
            result[rid] = tuple(value)
    return result


def exact_pool(records):
    summed = {}
    for rid, record in records.items():
        value = sum((Fraction.from_float(coefficient) for mid, coefficient in record['stoich'].items() if mid in POOL), Fraction(0))
        if value:
            summed[rid] = value
    epsilon = -Fraction.from_float(records['Growth']['stoich']['q8h2_c'])
    require(summed == {'CHRPL': 1, EX: -1, COUM: -1, 'Growth': -epsilon, '4HBHYOX': -1, 'sink_2ohph_c': -1}, 'Exact precursor row combination differs')
    require(all(records[rid]['bounds'][0] == 0 for rid in ('CHRPL', '4HBHYOX', 'sink_2ohph_c')), 'Necessary supply bound lacks nonnegative source/disposal terms')
    return float(epsilon), summed


def primal_certificate(records, bounds, fluxes):
    require(set(fluxes) == set(records), 'Full primal reaction inventory differs')
    require(all(finite(value) for value in fluxes.values()), 'Nonfinite full primal flux')
    terms = {}
    for rid, record in records.items():
        for mid, coefficient in record['stoich'].items():
            terms.setdefault(mid, []).append(coefficient * fluxes[rid])
    residuals = {mid: math.fsum(values) for mid, values in terms.items()}
    # fsum and a sequential sparse sum can differ by bounded roundoff.
    roundoff_bound = max((2 * np.finfo(float).eps * len(values) * math.fsum(abs(v) for v in values)
                         for values in terms.values()), default=0.)
    violations = {rid: max(lower - fluxes[rid], fluxes[rid] - upper, 0.) for rid, (lower, upper) in bounds.items()}
    certificate = {'max_abs_S_residual': max(map(abs, residuals.values()), default=0.),
                   'n_S_rows_over_1e-6': sum(abs(value) > 1e-6 for value in residuals.values()),
                   'max_bound_violation': max(violations.values(), default=0.),
                   'n_bounds_over_1e-6': sum(value > 1e-6 for value in violations.values())}
    return certificate, residuals, roundoff_bound


def verify_case(row, expected, records, base_bounds, aliases, recipe, epsilon):
    require(row['specification'] == expected, f'Case specification differs: {expected}')
    bounds = restricted_bounds(records, base_bounds, expected, aliases)
    require([solve['solver'] for solve in row['solves']] == ['glpk', 'highs'], 'Solver inventory differs')
    objective_id = expected.get('objective', 'Growth')
    checks = []
    for solve in row['solves']:
        require(solve['status'].lower() == 'optimal' and finite(solve['objective']), 'Failed/nonfinite saved solve')
        require(solve['feasibility_tolerance'] == recipe['solver_tolerance'], 'Solver tolerance differs')
        fluxes = solve['full_primal_fluxes']
        actual, residuals, roundoff = primal_certificate(records, bounds, fluxes)
        reported = solve['primal_certificate']
        require(set(reported) == set(actual), 'Primal certificate inventory differs')
        for field, value in actual.items():
            require(finite(reported[field]) and reported[field] >= 0, 'Invalid reported primal certificate')
            if field.startswith('n_'):
                require(type(reported[field]) is int and reported[field] == value, f'Certificate count differs: {field}')
            else:
                allowance = max(1e-13, roundoff) if field == 'max_abs_S_residual' else 1e-13
                require(abs(reported[field] - value) <= allowance, f'Certificate scalar differs: {field}')
                require(value <= recipe['residual_tolerance'] and reported[field] <= recipe['residual_tolerance'], f'Primal feasibility fails: {field}')
        require(near(solve['objective'], fluxes[objective_id], recipe['objective_tolerance']), 'Primal objective differs')
        require(near(solve['growth_flux'], fluxes['Growth'], recipe['objective_tolerance']), 'Primal Growth differs')
        require(not any(rid.startswith('DIAG_') for rid in fluxes), 'Unexpected diagnostic reaction in growth LP')
        require(set(solve['quinone_balance_residuals']) == {'q8_c', 'q8h2_c'}, 'Quinone residual inventory differs')
        for mid, value in solve['quinone_balance_residuals'].items():
            require(finite(value) and abs(value - residuals[mid]) <= max(roundoff, 1e-13), 'Quinone residual differs')
        pool_residual = math.fsum((fluxes['CHRPL'], -fluxes[EX], -fluxes[COUM], -epsilon * fluxes['Growth'], -fluxes['4HBHYOX'], -fluxes['sink_2ohph_c']))
        require(abs(pool_residual) <= recipe['residual_tolerance'], 'Exact precursor balance fails')
        checks.append({'solver': solve['solver'], 'recomputed_primal_certificate': actual,
                       'precursor_balance_residual': pool_residual, 'certificate_roundoff_allowance': roundoff})
    first, second = row['solves']
    require(near(first['objective'], second['objective'], recipe['objective_tolerance']), 'Solver objectives disagree')
    require(near(row['objective'], first['objective']) and near(row['growth'], first['growth_flux']), 'Top-level scalars differ')
    expected_selected = {rid: first['full_primal_fluxes'][rid] for rid in recipe['witness_reactions']}
    require(row['selected_glpk_fluxes'] == expected_selected, 'Selected GLPK witness differs from full vector')
    require(finite(row['precursor_balance_residual']) and abs(row['precursor_balance_residual'] - checks[0]['precursor_balance_residual']) <= 1e-12, 'Saved precursor residual differs')
    return {'case': expected, 'n_reactions': len(records),
            'n_metabolite_rows': len({mid for record in records.values() for mid in record['stoich']}),
            'solvers': checks, 'objective_disagreement': abs(first['objective'] - second['objective'])}


def verify_results(root, study):
    manifest = read_json(study / 'manifest.json')
    frozen = verify_manifest(root, manifest)
    summary = read_json(study / 'summary.json')
    recipe = manifest['plan']['recipe']
    require(summary['study'] == recipe and summary['study_fingerprint'] == frozen['content_fingerprint'], 'Summary study identity differs')
    require(recipe['study_id'] == 'putida-precursor-supply-diagnostic-v1' and recipe['evaluation_role'] == 'development', 'Wrong study identity/role')
    require(recipe['protocol']['complete_medium_transport'] and not recipe['protocol'].get('knockout_genes'), 'Unsupported baseline preparation')
    require(recipe['growth_fractions'] == [.1, .5, .95] and recipe['dose_multipliers'] == [0, .1, .5, 1, 2], 'Unexpected diagnostic grid')
    base_path = root / recipe['base_model']
    require(sha(base_path) == recipe['base_model_sha256'], 'Parent model hash differs')
    model = cobra.io.read_sbml_model(str(base_path))
    base, metabolites, exchanges = network_from_model(model)
    epsilon, exact = exact_pool(base)
    require(summary['epsilon'] == epsilon and finite(summary['seconds']) and summary['seconds'] >= 0, 'Invalid summary epsilon/timing')
    mapping = read_json(root / recipe['parent_gene_map'])['model_to_browser']
    aliases = {}
    for locus in recipe['gene_loci']:
        hits = [gene for gene, name in mapping.items() if name == locus]
        require(len(hits) == 1 and hits[0] in model.genes and model.genes.get_by_id(hits[0]).reactions, 'Absent/ambiguous diagnostic locus')
        aliases[locus] = hits[0]
    require(len(set(aliases.values())) == len(aliases), 'Gene aliases collapse distinct loci')
    conditions, media = metadata_conditions(root)
    require(len(conditions) == 43, 'Mapped condition inventory changed')
    contexts = read_json(study / 'contexts.json')
    require(len(contexts) == len(conditions) and len(summary['condition_results']) == len(conditions), 'Condition summary/context inventory differs')
    with gzip.open(study / 'cases.jsonl.gz', 'rt') as stream:
        rows = [json.loads(line, object_pairs_hook=pairs_unique,
                           parse_constant=lambda value: (_ for _ in ()).throw(ValueError(f'Nonfinite JSON: {value}')))
                for line in stream]
    checks, minimum_checks, rebuilt_outcomes, prepared = [], [], [], {}
    position = 0
    def consume(expected, records, bounds):
        nonlocal position
        require(position < len(rows), 'Missing declared case')
        row = rows[position]
        checks.append(verify_case(row, expected, records, bounds, aliases, recipe, epsilon))
        position += 1
        return row
    for condition, saved_context in zip(conditions, contexts):
        key = condition['key']
        records, bounds, context = prepare_network(base, metabolites, exchanges, conditions, media, condition, recipe)
        require(context == saved_context, f'Medium context differs: {key}')
        require(exact_pool(records)[1] == exact, 'Medium completion changes precursor certificate')
        wt = consume({'id': 'wild_type', 'kind': 'growth', 'condition': key}, records, bounds)
        ko = consume({'id': 'pp5317', 'kind': 'growth', 'condition': key, 'loci': ['PP_5317']}, records, bounds)
        require(ko['growth'] <= wt['growth'] + recipe['objective_tolerance'], 'KO exceeds parent optimum')
        outcome = {'condition': key, 'wild_type': wt['growth'], 'pp5317': ko['growth'], 'minimum_uptake': []}
        if wt['growth'] < recipe['protocol']['growth_threshold']:
            outcome['minimum_uptake_omission'] = 'Wild-type model does not grow at the declared threshold'
        elif EX in condition['exchanges']:
            outcome['minimum_uptake_omission'] = '4HBZ is already a main carbon substrate; a trace-only cap would confound carbon starvation with precursor supply'
        else:
            for fraction in recipe['growth_fractions']:
                target = fraction * wt['growth']
                spec = {'id': f'minimum_{fraction}', 'kind': 'minimum_uptake', 'condition': key,
                        'loci': ['PP_5317'], 'objective': EX, 'fraction': fraction, 'target_growth': target,
                        'growth_bounds': [target, target], 'exchange_bounds': [-recipe['uptake_capacity'], 0]}
                row = consume(spec, records, bounds)
                uptake = max(0., -row['objective'])
                coumarate = bounds[COUM][0] < 0
                require(0 <= uptake <= recipe['uptake_capacity'] + recipe['objective_tolerance'], 'Minimum exceeds capacity')
                if not coumarate:
                    require(uptake + recipe['residual_tolerance'] >= epsilon * target, 'Minimum violates necessary precursor lower bound')
                minimum_checks.append({'condition': key, 'fraction': fraction, 'minimum': uptake,
                                       'epsilon_times_target': epsilon * target, 'coumarate_import_allowed': coumarate,
                                       'bound_applicable': not coumarate,
                                       'gap_above_bound': None if coumarate else uptake - epsilon * target})
                outcome['minimum_uptake'].append({'fraction': fraction, 'target_growth': target, 'minimum_uptake': uptake,
                    'epsilon_times_target': epsilon * target, 'coumarate_uptake_allowed': coumarate})
        rebuilt_outcomes.append(outcome)
        prepared[key] = records, bounds, wt, ko
    require(rebuilt_outcomes == summary['condition_results'], 'Summary condition/minimum values differ')
    key = recipe['diagnostic_condition']['name'] + ' | ' + recipe['diagnostic_condition']['media']
    records, bounds, wt, ko = prepared[key]
    require(wt['growth'] >= recipe['protocol']['growth_threshold'] and bounds[COUM][0] == 0, 'Invalid glucose dose context')
    doses = []
    for multiplier in recipe['dose_multipliers']:
        cap = multiplier * epsilon * wt['growth']
        row = consume({'id': f'dose_{multiplier}', 'kind': 'dose_growth', 'condition': key,
                       'loci': ['PP_5317'], 'dose_multiplier': multiplier, 'exchange_bounds': [-cap, 1000]}, records, bounds)
        require(epsilon * row['growth'] <= cap + recipe['residual_tolerance'], 'Dose response violates precursor upper bound')
        require(not doses or row['growth'] >= doses[-1]['growth'] - recipe['objective_tolerance'], 'Nested dose capacities reduce maximum growth')
        doses.append({'multiplier': multiplier, 'capacity': cap, 'growth': row['growth']})
    controls = {}
    for control in recipe['growth_controls']:
        spec = {'kind': 'growth_control', 'condition': key, 'exchange_bounds': [-recipe['uptake_capacity'], 1000], **control}
        controls[control['id']] = consume(spec, records, bounds)['growth']
    require(controls['supply_wild_type'] + recipe['objective_tolerance'] >= wt['growth'], 'Expanded WT supply lowers optimum')
    require(controls['supply_pp5317'] <= controls['supply_wild_type'] + recipe['objective_tolerance'], 'Supply KO exceeds supply WT')
    for control in recipe['growth_controls']:
        if control['id'] not in ('supply_wild_type', 'supply_pp5317'):
            require(controls[control['id']] <= controls['supply_pp5317'] + recipe['objective_tolerance'], 'Restricted control exceeds supply KO')
    donors = {}
    for control in recipe['donor_controls']:
        spec = {'kind': 'donor_secretion', 'condition': key, 'objective': EX,
                'growth_bounds': [recipe['donor_growth_fraction'] * wt['growth'], bounds['Growth'][1]],
                'exchange_bounds': [0, 1000], **control}
        donors[control['id']] = consume(spec, records, bounds)['objective']
    require(donors['donor_transport_closed'] <= donors['donor_wild_type'] + recipe['objective_tolerance'], 'Restricted donor exceeds WT donor')
    require(position == len(rows) == summary['n_cases'] and summary['n_solver_runs'] == 2 * position, 'Unexpected extra case or total')
    require([r['context'] for r in summary['energy_checks']] == ['parent', 'hypothetical_supply'], 'Energy context inventory differs')
    for energy in summary['energy_checks']:
        require(set(energy['values']) == {'atp', 'nadh', 'nadph', 'q8h2', 'h_p'}, 'Energy probe inventory differs')
        require(all(finite(value) and value >= 0 and value <= recipe['production_threshold'] for value in energy['values'].values()), 'Invalid saved energy result')
    require(summary['energy_checks'][0]['values'] == summary['energy_checks'][1]['values'], 'Closed-boundary energy repetitions differ')
    after = read_json(study / 'input_verification_after_run.json')
    require(after['valid'] and after['checked_files'] == frozen['checked_files'] and after['content_fingerprint'] == frozen['content_fingerprint'], 'Post-run input verification differs')
    verify_manifest(root, manifest)
    all_primal = [solve['recomputed_primal_certificate'] for check in checks for solve in check['solvers']]
    added_reactions = sorted(set(records) - set(base))
    added_boundaries = [rid for rid in added_reactions if len(records[rid]['stoich']) == 1]
    saved_exchange_ids = set(contexts[0]['exchange_bounds'])
    omitted_boundaries = [rid for rid in added_boundaries if rid not in saved_exchange_ids]
    return {'valid': True, 'n_conditions': len(conditions), 'n_cases': position, 'n_reconstructed_full_primal_vectors': 2 * position,
            'input_verification': frozen, 'epsilon': epsilon, 'gene_aliases': aliases,
            'medium_completion_inventory': {'added_reactions': added_reactions, 'added_boundaries': added_boundaries,
                'omitted_from_saved_exchange_inventory': omitted_boundaries,
                'note': 'Full-vector checks include every added exchange and carrier. Newly created extracellular metabolites use e, whereas the parent external compartment is C_e; COBRA omits these new EX reactions from model.exchanges and the serialized context inventory. Their complete LP bounds are independently reconstructed and checked here. Closed-boundary energy probes use reaction.boundary, which includes these reactions.'},
            'exact_precursor_coefficients': {rid: str(value) for rid, value in exact.items()},
            'max_recomputed_S_residual': max(r['max_abs_S_residual'] for r in all_primal),
            'max_recomputed_bound_violation': max(r['max_bound_violation'] for r in all_primal),
            'max_objective_disagreement': max(r['objective_disagreement'] for r in checks),
            'minimum_checks': minimum_checks, 'dose_results': doses, 'growth_controls': controls, 'donor_objectives': donors,
            'case_checks': checks,
            'scope': 'All declared full primal reaction vectors independently checked against reconstructed equations/bounds; no optimizations or phenotype-value reads. Status and solver agreement are not a dual certificate or exact optimality proof. Energy outputs checked as reported scalars only; their repeated closed-boundary contexts do not establish sharing. Feasible precursor supply/secretion is not measured cross-feeding, carryover, concentration or monoculture viability.'}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True, help='Fresh verification directory')
    args = parser.parse_args(argv)
    study, out = args.study.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    # Preserve the verifier before inspecting any numerical result.
    source = Path(__file__).resolve()
    shutil.copyfile(source, out / 'verifier_source.py')
    write_json(out / 'verifier_provenance.json', {'created_at': datetime.now(timezone.utc).isoformat(),
        'source_sha256': sha(source), 'python': sys.version, 'platform': platform.platform(),
        'versions': {name: version(name) for name in ('cobra', 'numpy')},
        'scope': 'Independent artifact verification; no solver calls.'})
    try:
        paths = [study / name for name in ('manifest.json', 'summary.json', 'contexts.json', 'cases.jsonl.gz', 'input_verification_after_run.json')]
        inputs = [{'path': str(path), 'sha256': sha(path), 'size_bytes': path.stat().st_size} for path in paths]
        write_json(out / 'result_inputs.json', inputs)
        report = verify_results(ROOT, study)
        require(all(sha(Path(row['path'])) == row['sha256'] for row in inputs) and sha(source) == sha(out / 'verifier_source.py'), 'Verification inputs changed')
        write_json(out / 'report.json', report)
        print(json.dumps({key: report[key] for key in ('valid', 'n_conditions', 'n_cases', 'n_reconstructed_full_primal_vectors', 'max_recomputed_S_residual', 'max_recomputed_bound_violation', 'max_objective_disagreement')}, indent=2))
        return 0
    except Exception as error:
        write_json(out / 'failure.json', {'type': type(error).__name__, 'message': str(error)})
        raise


if __name__ == '__main__':
    raise SystemExit(main())
