"""Quantify hypothetical precursor supply; do not patch or rescore any model.

The diagnostic uses only already represented transport and reactions. Feasible
secretion/import is not evidence that metabolites were shared in an experiment.
"""
from __future__ import annotations

import argparse
import gzip
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
from gembench.checks import energy_from_nothing
from gembench.fitness_browser import base_medium, carbon_source_conditions
from gembench.media import apply_medium
from gembench.protocols.carbon_fitness_generic import GenericParams, complete_medium_transport
from gembench.study import freeze_study, verify_study
from scripts.diagnose_quinone_producibility import matrix_model, W
from scripts.map_quinone_precursor import BASE, BASE_HASH, PRECURSOR, TAIL, COUMARATE, combined_rows
from scripts.run_quinone_repair import input_path, metadata_organism, read_json, write_json
from scripts.verify_quinone_repair_solvers import assert_ordinary_lp

EXCHANGE = 'EX_4hbz_e'
PLAN = ROOT / 'data/studies/precursor_supply_v1.json'


def evidence_hashes(record):
    """Accept the two explicit provenance schemas used by the source audits."""
    hashes = dict(record.get('input_sha256', {}))
    for item in record.get('provenance', {}).get('inputs', []):
        path, digest = item['path'], item['sha256']
        if path in hashes and hashes[path] != digest:
            raise ValueError(f'Conflicting evidence hashes: {path}')
        hashes[path] = digest
    if not hashes:
        raise ValueError('Evidence record has no source hashes')
    return hashes


def freeze(root, recipe):
    old = read_json(root / recipe['parent_manifest'])
    registry = read_json(root / 'data/studies/exposure_registry_v1.json')
    paths = {r['path'] for r in old['files']}
    paths.update([BASE, recipe['parent_manifest'], recipe['parent_gene_map'],
                  'scripts/run_precursor_supply.py', 'scripts/map_quinone_precursor.py',
                  'tests/test_precursor_supply.py', 'data/studies/precursor_supply_v1.json',
                  'docs/studies/precursor-supply-v1.md'])
    paths.update(str(p.relative_to(root)) for p in (root / recipe['evidence_directory']).rglob('*') if p.is_file())
    for evidence in recipe['evidence_inputs']:
        record = read_json(root / evidence)
        paths.update(evidence_hashes(record))
    return freeze_study(root, {'schema_version': 1, 'study_id': recipe['study_id'],
        'evaluation_role': 'development', 'development_organisms': registry['development_organisms'],
        'quarantined_organisms': registry['quarantined_organisms'], 'paths': sorted(paths),
        'recipe': recipe, 'runtime': {**software_versions(), 'glpk': swiglpk.glp_version(), 'highs': highspy.Highs().version()},
        'note': 'Frozen before diagnostic solves; question follows exposed development errors. No independent validation claimed.'})


def verify(root, manifest):
    record = verify_study(root, manifest)
    if not record['valid']:
        raise ValueError(f'Frozen inputs changed: {record}')
    return record


def prepare(base, condition, all_conditions, recipe):
    model = base.copy()
    params = GenericParams(**recipe['protocol'])
    model.solver = 'glpk'
    model.solver.configuration.tolerances.feasibility = recipe['solver_tolerance']
    for r in model.exchanges:
        r.bounds = (0, 1000)
    added = complete_medium_transport(model, sorted({c.media for c in all_conditions if c.bigg_ids}), params.medium_completion_exclude)
    missing = apply_medium(model, base_medium(condition.media), close_all=True)
    absent_carbon = []
    for rid in condition.exchanges:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).lower_bound = params.carbon_uptake
        else:
            absent_carbon.append(rid)
    context = {'condition': condition.key, 'carbon_exchanges': condition.exchanges,
        'medium_completion_added': added, 'missing_medium': missing, 'missing_carbon': absent_carbon,
        'exchange_bounds': {r.id: list(r.bounds) for r in sorted(model.exchanges, key=lambda r: r.id)},
        'growth_bounds': list(model.reactions.Growth.bounds), 'atpm_bounds': list(model.reactions.ATPM.bounds)}
    assert_ordinary_lp(model)
    return model, context


def apply_case(model, specification, aliases):
    """Apply only declared bounds, deletions and a single-reaction objective."""
    allowed = {'id', 'kind', 'loci', 'close_reactions', 'exchange_bounds', 'growth_bounds',
               'objective', 'fraction', 'dose_multiplier', 'condition', 'target_growth'}
    if set(specification) - allowed:
        raise ValueError('Unknown case specification field')
    for locus in specification.get('loci', []):
        gid = aliases[locus]
        if gid not in model.genes or not model.genes.get_by_id(gid).reactions:
            raise ValueError(f'Unrepresented diagnostic target: {locus}')
        model.genes.get_by_id(gid).knock_out()
    for rid in specification.get('close_reactions', []):
        if rid not in model.reactions:
            raise ValueError(f'Unrepresented diagnostic reaction: {rid}')
        model.reactions.get_by_id(rid).bounds = (0, 0)
    if 'exchange_bounds' in specification:
        bounds = specification['exchange_bounds']
        if len(bounds) != 2 or not all(math.isfinite(x) for x in bounds) or bounds[0] > bounds[1]:
            raise ValueError('Invalid diagnostic exchange bounds')
        model.reactions.get_by_id(EXCHANGE).bounds = tuple(bounds)
    if 'growth_bounds' in specification:
        model.reactions.Growth.bounds = tuple(specification['growth_bounds'])
    model.objective = specification.get('objective', 'Growth')
    model.objective_direction = 'max'
    assert_ordinary_lp(model)


def solve_case(model, solver, tolerance):
    """Retain both complete primal witnesses for independent artifact checks."""
    matrix = matrix_model(model)
    if solver == 'glpk':
        model.solver = 'glpk'
        model.solver.configuration.tolerances.feasibility = tolerance
        solution = model.optimize()
        status = model.solver.status
        values = solution.fluxes.reindex(matrix.rxns).to_numpy() if status == 'optimal' else None
        objective = solution.objective_value if status == 'optimal' else None
    elif solver == 'highs':
        answer = W.solve_highs(matrix, method='simplex', threads=0, feas_tol=tolerance,
                               opt_tol=tolerance, time_limit=60.0, extra_options={'parallel': 'off'})
        status = answer.status
        values = answer.x if status == 'Optimal' else None
        objective = answer.objective if status == 'Optimal' else None
    else:
        raise ValueError(f'Unsupported diagnostic solver: {solver}')
    record = {'solver': solver, 'feasibility_tolerance': tolerance,
              'status': status, 'objective': objective}
    if values is not None:
        fluxes = {str(rid): float(value) for rid, value in zip(matrix.rxns, values)}
        if not all(math.isfinite(value) for value in fluxes.values()):
            raise RuntimeError('Nonfinite full primal witness')
        record.update(primal_certificate=W.certify(matrix, values),
                      growth_flux=fluxes['Growth'], full_primal_fluxes=fluxes)
        residual = matrix.S @ values
        record['quinone_balance_residuals'] = {
            mid: float(residual[list(matrix.mets).index(mid)]) for mid in ('q8_c', 'q8h2_c')}
    return record


def checked_solve(model, recipe):
    answers, selected = [], None
    for solver in ('glpk', 'highs'):
        answer = solve_case(model, solver, recipe['solver_tolerance'])
        if answer['status'].lower() != 'optimal' or answer['objective'] is None or not math.isfinite(answer['objective']):
            raise RuntimeError(f'Diagnostic LP did not solve optimally: {solver}/{answer["status"]}')
        certificate = answer.get('primal_certificate', {})
        for field in ('max_abs_S_residual', 'max_bound_violation'):
            value = certificate.get(field)
            if value is None or not math.isfinite(value) or value > recipe['residual_tolerance']:
                raise RuntimeError(f'Failed numerical certificate: {solver}/{field}')
        if solver == 'glpk':
            selected = {rid: float(model.reactions.get_by_id(rid).flux) for rid in recipe['witness_reactions']}
            if not all(math.isfinite(x) for x in selected.values()):
                raise RuntimeError('Nonfinite selected flux witness')
        answers.append(answer)
    if not math.isclose(answers[0]['objective'], answers[1]['objective'], abs_tol=recipe['objective_tolerance'], rel_tol=recipe['objective_tolerance']):
        raise RuntimeError('Independent solver objectives disagree')
    epsilon = -model.reactions.Growth.get_coefficient('q8h2_c')
    balance = selected['CHRPL'] - selected[EXCHANGE] - selected['EX_T4hcinnm_e'] - epsilon * selected['Growth'] - selected['4HBHYOX'] - selected['sink_2ohph_c']
    if abs(balance) > recipe['residual_tolerance']:
        raise RuntimeError('Selected witness violates exact precursor balance')
    return {'objective': answers[0]['objective'], 'growth': answers[0]['growth_flux'],
        'solves': answers, 'selected_glpk_fluxes': selected, 'precursor_balance_residual': balance,
        'witness_scope': 'One optimal GLPK solution; individual fluxes need not equal another optimum.'}


def minimum_spec(condition, fraction, wild_type, capacity):
    target = fraction * wild_type
    if not 0 < fraction < 1 or not math.isfinite(target) or target <= 0 or not math.isfinite(capacity) or capacity <= 0:
        raise ValueError('Require a finite positive target below the wild type and positive capacity')
    return {'id': f'minimum_{fraction}', 'kind': 'minimum_uptake', 'condition': condition,
        'loci': ['PP_5317'], 'objective': EXCHANGE, 'fraction': fraction, 'target_growth': target,
        'growth_bounds': [target, target], 'exchange_bounds': [-capacity, 0]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--freeze-only', action='store_true')
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('Preserve all attempts; select a new output directory')
    args.out.mkdir(parents=True)
    started = time.time()
    try:
        recipe = read_json(PLAN)
        if recipe['evaluation_role'] != 'development' or recipe['base_model'] != BASE or recipe['base_model_sha256'] != BASE_HASH:
            raise ValueError('Unexpected study identity or parent')
        verify(ROOT, read_json(ROOT / recipe['parent_manifest']))
        for filename in recipe['evidence_inputs']:
            evidence = read_json(ROOT / filename)
            for path, digest in evidence_hashes(evidence).items():
                if sha256_of(input_path(ROOT, path)) != digest:
                    raise ValueError(f'Evidence source changed: {path}')
        manifest = freeze(ROOT, recipe)
        write_json(args.out / 'manifest.json', manifest)
        if args.freeze_only:
            return
        if sha256_of(ROOT / BASE) != BASE_HASH:
            raise ValueError('Parent model hash changed')
        base = cobra.io.read_sbml_model(str(ROOT / BASE))
        context_evidence = read_json(ROOT / recipe['model_context'])
        rows = combined_rows(base, PRECURSOR + TAIL + COUMARATE)
        if rows != context_evidence['row_combinations']['total']['nonzero_rows']:
            raise ValueError('Exact precursor certificate changed')
        mapping = read_json(ROOT / recipe['parent_gene_map'])['model_to_browser']
        aliases = {}
        for locus in recipe['gene_loci']:
            hits = [g for g, value in mapping.items() if value == locus]
            if len(hits) != 1:
                raise ValueError(f'Ambiguous/absent target: {locus}')
            aliases[locus] = hits[0]
        metadata = metadata_organism(ROOT)
        all_conditions = carbon_source_conditions(metadata)
        conditions = [c for c in all_conditions if c.bigg_ids]
        cases, contexts, outcomes, energies = [], [], [], []
        epsilon = -base.reactions.Growth.get_coefficient('q8h2_c')
        selected_model = None
        selected_wt = None
        with gzip.open(args.out / 'cases.jsonl.gz', 'xt') as stream:
            def emit(model, specification):
                verify(ROOT, manifest)
                with model:
                    apply_case(model, specification, aliases)
                    result = checked_solve(model, recipe)
                    row = {'specification': specification, **result}
                    cases.append(row)
                    stream.write(json.dumps(row, allow_nan=False) + '\n')
                    stream.flush()
                    return row

            for condition in conditions:
                model, context = prepare(base, condition, all_conditions, recipe)
                contexts.append(context)
                wt = emit(model, {'id': 'wild_type', 'kind': 'growth', 'condition': condition.key})
                ko = emit(model, {'id': 'pp5317', 'kind': 'growth', 'condition': condition.key, 'loci': ['PP_5317']})
                result = {'condition': condition.key, 'wild_type': wt['growth'], 'pp5317': ko['growth'], 'minimum_uptake': []}
                if wt['growth'] < recipe['protocol']['growth_threshold']:
                    result['minimum_uptake_omission'] = 'Wild-type model does not grow at the declared threshold'
                elif EXCHANGE in condition.exchanges:
                    result['minimum_uptake_omission'] = '4HBZ is already a main carbon substrate; a trace-only cap would confound carbon starvation with precursor supply'
                else:
                    for fraction in recipe['growth_fractions']:
                        spec = minimum_spec(condition.key, fraction, wt['growth'], recipe['uptake_capacity'])
                        answer = emit(model, spec)
                        if answer['objective'] > recipe['objective_tolerance']:
                            raise RuntimeError('Uptake-only minimum returned secretion')
                        result['minimum_uptake'].append({'fraction': fraction, 'target_growth': spec['target_growth'],
                            'minimum_uptake': max(0, -answer['objective']), 'epsilon_times_target': epsilon * spec['target_growth'],
                            'coumarate_uptake_allowed': model.reactions.EX_T4hcinnm_e.lower_bound < 0})
                outcomes.append(result)
                if condition.name == recipe['diagnostic_condition']['name'] and condition.media == recipe['diagnostic_condition']['media']:
                    selected_model, selected_wt = model, wt['growth']
                print('Condition diagnostics complete:', condition.key, flush=True)
            if selected_model is None or selected_wt < recipe['protocol']['growth_threshold']:
                raise ValueError('Declared glucose diagnostic condition is not growing or absent')
            key = ' | '.join(recipe['diagnostic_condition'][k] for k in ('name', 'media'))
            for multiplier in recipe['dose_multipliers']:
                emit(selected_model, {'id': f'dose_{multiplier}', 'kind': 'dose_growth', 'condition': key,
                    'loci': ['PP_5317'], 'dose_multiplier': multiplier,
                    'exchange_bounds': [-multiplier * epsilon * selected_wt, 1000]})
            for control in recipe['growth_controls']:
                emit(selected_model, {'kind': 'growth_control', 'condition': key,
                    'exchange_bounds': [-recipe['uptake_capacity'], 1000], **control})
            for control in recipe['donor_controls']:
                emit(selected_model, {'kind': 'donor_secretion', 'condition': key, 'objective': EXCHANGE,
                    'growth_bounds': [recipe['donor_growth_fraction'] * selected_wt, selected_model.reactions.Growth.upper_bound],
                    'exchange_bounds': [0, 1000], **control})
            for label, bound in [('parent', 0), ('hypothetical_supply', -recipe['uptake_capacity'])]:
                with selected_model:
                    selected_model.reactions.get_by_id(EXCHANGE).lower_bound = bound
                    values = energy_from_nothing(selected_model)
                    if any(v > recipe['production_threshold'] for v in values.values()):
                        raise ValueError('Energy-from-nothing diagnostic failed')
                    energies.append({'context': label, 'values': values,
                        'note': 'The energy test itself closes all boundaries, including the hypothetical supply.'})
        write_json(args.out / 'contexts.json', contexts)
        write_json(args.out / 'input_verification_after_run.json', verify(ROOT, manifest))
        write_json(args.out / 'summary.json', {'study': recipe, 'study_fingerprint': manifest['content_fingerprint'],
            'condition_results': outcomes, 'energy_checks': energies, 'n_cases': len(cases),
            'n_solver_runs': 2 * len(cases), 'epsilon': epsilon, 'seconds': time.time() - started,
            'scope': 'Hypothetical supply and secretion feasibility only. No model patch, new reaction, fitness rescore, concentration estimate or evidence of observed cross-feeding.'})
    except Exception as error:
        write_json(args.out / 'failure.json', {'type': type(error).__name__, 'message': str(error),
            'seconds': time.time() - started, 'note': 'Partial outputs preserved; any reviewed retry requires a fresh directory.'})
        raise


if __name__ == '__main__':
    main()
