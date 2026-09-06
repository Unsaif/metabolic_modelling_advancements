"""Independently check the fixed 43-condition numerical comparison artifacts.

No solver, new preparation helper, runner or shared certificate function is
imported. Reconstructs the historical LP and explicit medium oracle, both sparse
and compensated mass balances, objective consistency and the declared decisions.
Rejected and unresolved numerical results remain visible rather than aborting
verification merely for failing their declared acceptance gate.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import gzip
from itertools import combinations
import json
import logging
import math
from pathlib import Path
import sys
import time
import warnings

import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import verify_medium_preparation as V
ROOT = V.ROOT
METHODS = ('glpk_native', 'highs_primal_simplex', 'highs_simplex_tighter')


def close(a, b, tolerance):
    return V.finite(a) and V.finite(b) and math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)


def read_gzip(path):
    with gzip.open(path, 'rt') as handle:
        return V.loads(handle.read())


def assert_number(observed, expected, label):
    V.require(close(observed, expected, 1e-12), f'Numeric artifact mismatch: {label}')


def ranked_subset(observed, candidates, identity, score, limit, label):
    """Allow equal-score ties without trusting set iteration or ordering."""
    V.require(len(observed) == min(limit, len(candidates)), f'Wrong diagnostic subset length: {label}')
    selected = [row[identity] for row in observed]
    V.require(len(set(selected)) == len(selected) and set(selected) <= set(candidates), f'Invalid diagnostic identities: {label}')
    weights = [score(candidates[key]) for key in selected]
    V.require(weights == sorted(weights, reverse=True), f'Diagnostic ranking is not descending: {label}')
    if selected:
        V.require(min(weights) >= max((score(row) for key, row in candidates.items() if key not in selected), default=-math.inf),
                  f'Diagnostic ranking omits a larger value: {label}')


def audit_vector(signature, record, recipe):
    rxns, metabolites = signature['reactions'], signature['metabolites']
    fluxes = record.get('full_primal_fluxes')
    complete = isinstance(fluxes, dict) and set(fluxes) == set(rxns) and all(V.finite(x) for x in fluxes.values())
    if not complete:
        V.require(record.get('finite_primal_vector', False) is False, 'Incomplete/nonfinite vector reported as finite')
        return None
    V.require(record.get('finite_primal_vector') is True, 'Complete finite vector incorrectly reported')
    rids, mids = list(rxns), list(metabolites)
    mi = {mid: index for index, mid in enumerate(mids)}
    rows, cols, values = [], [], []
    terms = {mid: {} for mid in mids}
    for column, rid in enumerate(rids):
        for mid, coefficient in rxns[rid]['stoichiometry'].items():
            V.require(mid in mi and V.finite(coefficient), 'Invalid reconstructed stoichiometry')
            rows.append(mi[mid]); cols.append(column); values.append(coefficient)
            terms[mid][rid] = {'reaction': rid, 'coefficient': coefficient,
                              'flux': fluxes[rid], 'term': coefficient * fluxes[rid]}
    # Independent matrix assembly; no runner/certificate translator is reused.
    matrix = csr_matrix((values, (rows, cols)), shape=(len(mids), len(rids)))
    sparse = np.asarray(matrix @ np.asarray([fluxes[rid] for rid in rids], dtype=float)).reshape(-1)
    residuals = {}
    row_records = {}
    for index, mid in enumerate(mids):
        compensated = math.fsum(row['term'] for row in terms[mid].values())
        residuals[mid] = {'sparse': float(sparse[index]), 'fsum': compensated}
        row_records[mid] = {'metabolite': mid, 'sparse_residual': float(sparse[index]),
            'fsum_residual': compensated, 'sum_abs_terms': math.fsum(abs(row['term']) for row in terms[mid].values()),
            'n_terms': len(terms[mid])}
    observed_residuals = record['complete_metabolite_residuals']
    V.require(set(observed_residuals) == set(mids), 'Incomplete metabolite residual vector')
    for mid, expected in residuals.items():
        V.require(set(observed_residuals[mid]) == {'sparse', 'fsum'}, 'Unexpected residual fields')
        for kind, value in expected.items():
            assert_number(observed_residuals[mid][kind], value, f'{mid}/{kind}')
    bounds = {rid: {'reaction': rid, 'flux': fluxes[rid], 'bounds': row['bounds'],
                   'violation': max(0., row['bounds'][0] - fluxes[rid], fluxes[rid] - row['bounds'][1])}
              for rid, row in rxns.items()}
    sparse_max = max(abs(row['sparse']) for row in residuals.values())
    compensated_max = max(abs(row['fsum']) for row in residuals.values())
    bound_max = max(row['violation'] for row in bounds.values())
    certificate = {'max_abs_S_residual': sparse_max, 'max_bound_violation': bound_max,
                   'n_S_rows_over_1e-6': sum(abs(row['sparse']) > 1e-6 for row in residuals.values()),
                   'n_bounds_over_1e-6': sum(row['violation'] > 1e-6 for row in bounds.values())}
    V.require(set(record['primal_certificate']) == set(certificate), 'Unexpected certificate fields')
    for key, value in certificate.items():
        assert_number(record['primal_certificate'][key], value, key)
    objective = math.fsum(coefficient * fluxes[rid] for rid, coefficient in signature['objective'].items())
    assert_number(record['recomputed_objective'], objective, 'recomputed objective')
    assert_number(record['max_fsum_residual'], compensated_max, 'max fsum residual')
    threshold = recipe['residual_tolerance']
    V.require(record['n_sparse_rows_above_original_tolerance'] == sum(abs(row['sparse']) > threshold for row in residuals.values()),
              'Wrong sparse residual exceedance count')
    V.require(record['n_fsum_rows_above_original_tolerance'] == sum(abs(row['fsum']) > threshold for row in residuals.values()),
              'Wrong compensated residual exceedance count')
    for key, score, limit in [('worst_bounds', lambda row: row['violation'], 20),
                              ('largest_absolute_fluxes', lambda row: abs(row['flux']), 30)]:
        ranked_subset(record[key], bounds, 'reaction', score, limit, key)
        V.require(all(row == bounds[row['reaction']] for row in record[key]), f'Incorrect diagnostic bound records: {key}')
    ranked_subset(record['worst_rows'], row_records, 'metabolite',
                  lambda row: max(abs(row['sparse_residual']), abs(row['fsum_residual'])), 20, 'worst_rows')
    for row in record['worst_rows']:
        expected = row_records[row['metabolite']]
        V.require(set(row) == set(expected) | {'largest_terms'}, 'Unexpected worst-row fields')
        for key, value in expected.items():
            if key in {'metabolite', 'n_terms'}:
                V.require(row[key] == value, 'Wrong row identity/term count')
            else:
                assert_number(row[key], value, key)
        possible = terms[row['metabolite']]
        ranked_subset(row['largest_terms'], possible, 'reaction', lambda item: abs(item['term']), 20, 'largest_terms')
        V.require(all(item == possible[item['reaction']] for item in row['largest_terms']), 'Incorrect diagnostic row terms')
    gate = sparse_max <= threshold and bound_max <= threshold
    V.require(record['passes_original_primal_gate'] == gate, 'Original primal gate flag differs')
    optimal = str(record.get('status', '')).lower() == 'optimal'
    V.require(record['passes_original_case_gate'] == (optimal and V.finite(record.get('objective')) and gate),
              'Original case gate flag differs')
    return {'certificate': certificate, 'max_fsum_residual': compensated_max, 'objective': objective,
            'n_fluxes': len(rids), 'n_balance_rows': len(mids)}


def method_assessment(record, audit, recipe):
    complete = audit is not None
    original = complete and all(audit['certificate'][key] <= recipe['residual_tolerance']
                               for key in ('max_abs_S_residual', 'max_bound_violation'))
    compensated = complete and audit['max_fsum_residual'] <= recipe['residual_tolerance']
    consistent = complete and close(record.get('objective'), audit['objective'], recipe['objective_tolerance'])
    checks = [('no_capture_exception', record.get('exception') is None),
              ('optimal_status', str(record.get('status', '')).lower() == 'optimal'),
              ('complete_finite_vector', complete), ('original_primal_gate', original),
              ('compensated_residual_gate', compensated), ('objective_matches_vector', consistent)]
    reasons = [name for name, valid in checks if not valid]
    return {'original_primal_gate_recomputed': original, 'complete_finite_vector': complete,
            'compensated_residual_gate': compensated, 'objective_matches_vector': consistent,
            'accepted_for_comparison': not reasons, 'failed_requirements': reasons}


def case_assessment(records, audits, recipe):
    methods = {record['method']['id']: method_assessment(record, audit, recipe) for record, audit in zip(records, audits)}
    optimal = [record for record in records if str(record.get('status', '')).lower() == 'optimal']
    pairs = [{'first': a['method']['id'], 'second': b['method']['id'], 'first_objective': a.get('objective'),
              'second_objective': b.get('objective'), 'agree': close(a.get('objective'), b.get('objective'), recipe['objective_tolerance'])}
             for a, b in combinations(optimal, 2)]
    finite = all(V.finite(record.get('objective')) for record in optimal)
    agree = len(optimal) >= 2 and finite and all(pair['agree'] for pair in pairs)
    statuses = {str(record.get('status', '')).lower() for record in records}
    feasibility_conflict = 'optimal' in statuses and any('infeasible' in status for status in statuses)
    conflict = 'optimal' in statuses and any('infeasible' in status or 'unbounded' in status for status in statuses)
    cross = methods[METHODS[0]]['accepted_for_comparison'] and any(methods[key]['accepted_for_comparison'] for key in METHODS[1:])
    return {'methods': methods, 'n_reported_optimal': len(optimal), 'all_optimal_objectives_finite': finite,
            'pairwise_reported_optimal_objectives': pairs, 'all_reported_optimal_objectives_agree': agree,
            'reported_feasibility_conflict': feasibility_conflict, 'reported_status_conflict': conflict,
            'glpk_and_at_least_one_highs_accepted': cross,
            'accepted_cross_solver_evidence': cross and agree and not conflict,
            'all_three_methods_accepted': all(item['accepted_for_comparison'] for item in methods.values()) and agree and not conflict}


def infeasibility_review(records, audits, recipe):
    """A solver status is distinct from a primal witness and from zero growth."""
    statuses = [str(record.get('status', '')).lower() for record in records]
    conflicts = []
    for record, audit, status in zip(records, audits, statuses):
        if status != 'infeasible' or audit is None:
            continue
        within = (audit['certificate']['max_abs_S_residual'] <= recipe['residual_tolerance'] and
                  audit['certificate']['max_bound_violation'] <= recipe['residual_tolerance'] and
                  audit['max_fsum_residual'] <= recipe['residual_tolerance'])
        if within:
            conflicts.append({'method': record['method']['id'], 'status': record['status'],
                              'certificate': audit['certificate'], 'max_fsum_residual': audit['max_fsum_residual'],
                              'interpretation': 'A tolerance-feasible complete witness conflicts with this infeasible status; this is not an exact-arithmetic feasibility proof.'})
    return all(status == 'infeasible' for status in statuses), conflicts


def verify(run):
    required = [run / name for name in ('manifest.json', 'summary.json', 'cases.jsonl.gz', 'source_context.json')]
    V.require(all(path.is_file() for path in required) and not (run / 'failure.json').exists(), 'Completed numerical artifacts required')
    direct = {str(path.relative_to(ROOT)): V.sha(path) for path in required}
    manifest, summary, context = V.read(required[0]), V.read(required[1]), V.read(required[3])
    frozen = V.verify_manifest(ROOT, manifest)
    recipe = manifest['plan']['recipe']
    V.require(summary['study'] == recipe and summary['study_fingerprint'] == manifest['content_fingerprint']
              and summary['completed'] is True, 'Numerical summary declaration differs')
    V.require(recipe['model_label'] == 'Putida_iJN1463' and recipe['condition_count'] == 43 and recipe['solve_count'] == 129
              and [m['id'] for m in recipe['methods']] == list(METHODS), 'Unexpected numerical panel')
    V.require(recipe['residual_tolerance'] == recipe['objective_tolerance'] == 1e-8, 'Original acceptance threshold changed')
    source = V.load_source()
    base, provenance = source.build_parent('Putida_iJN1463')
    V.require(context['provenance'] == provenance and all(frozen.get(path) == digest for path, digest in provenance['input_sha256'].items()),
              'Parent source reconstruction differs or has unfrozen dependencies')
    conditions = source.metadata_conditions('Putida')
    declared = [{'key': c.key, 'name': c.name, 'media': c.media, 'exchanges': c.exchanges} for c in conditions]
    V.require(context['conditions'] == declared and len(conditions) == 43, 'Condition metadata differs')
    V.require(context['effective_highs_options'] == manifest['plan']['effective_highs_options'], 'Effective option records differ')
    media = sorted({condition.media for condition in conditions})
    policy = V.policy_record(recipe['compartments'])
    V.require(policy['exchange_metabolites'] == {'EX_AEP_e': '2ameph_e'}, 'Source exchange binding differs')
    legacy = base.copy()
    V.complete_medium_transport(legacy, media, policy['completion_exclude'])
    expected_raw = set()
    outcomes, methods = [], {key: {'statuses': Counter(), 'original': 0, 'accepted': 0} for key in METHODS}
    finite_vectors, flux_count, row_count, boundary_count = 0, 0, 0, 0
    maxima = {'sparse_residual': 0., 'fsum_residual': 0., 'bound_violation': 0., 'objective_spread': 0.}
    per_method_maxima = {key: {'sparse_residual': 0., 'fsum_residual': 0., 'bound_violation': 0.} for key in METHODS}
    rejected_records = []
    unanimous_infeasible, status_witness_conflicts = [], []
    with gzip.open(required[2], 'rt') as handle:
        for index, line in enumerate(handle):
            V.require(index < len(conditions), 'Unexpected extra numerical case')
            case = V.loads(line)
            condition = conditions[index]
            V.require(case['index'] == index and case['model'] == recipe['model_label'] and case['condition'] == condition.key,
                      'Numerical cases are not the ordered declared inventory')
            medium = V.medium_for(condition)
            old = legacy.copy(); V.apply_medium(old, medium, close_all=True)
            prepared, report = V.independent_target(base, medium, media, policy)
            signature = V.canonical(prepared)
            V.require(case['prepared_signature'] == signature, 'Full saved prepared LP differs from independent oracle')
            V.verify_signatures(case, V.canonical(old), signature)
            V.check_report(case['medium_report'], report)
            boundary_count += len(report['boundary_reactions'])
            directory = run / 'raw_methods' / f'{index:02d}'
            raw_signature = directory / 'prepared_signature.json.gz'
            expected_raw.add(raw_signature)
            direct[str(raw_signature.relative_to(ROOT))] = V.sha(raw_signature)
            V.require(read_gzip(raw_signature) == case['prepared_signature'], 'Raw prepared signature differs')
            V.require([r['method'] for r in case['solves']] == recipe['methods'], 'Numerical methods/settings differ')
            audits = []
            for record in case['solves']:
                method = record['method']['id']
                raw = directory / f'{method}.json.gz'
                expected_raw.add(raw); direct[str(raw.relative_to(ROOT))] = V.sha(raw)
                V.require(read_gzip(raw) == record, 'Raw method result differs from case record')
                V.require(record['algebra_and_chemistry_sha256'] == case['algebra_and_chemistry_sha256'], 'Method LP hash differs')
                if record['method']['solver'] == 'glpk' and 'effective_settings' in record:
                    V.require(record['effective_settings']['feasibility_tolerance'] == record['method']['solver_tolerance']
                              == record['effective_settings']['simplex_tol_bnd'], 'Effective GLPK feasibility setting differs')
                if record['method']['solver'] == 'highs' and 'solver_info' in record:
                    V.require(float(record['solver_info']['feas_tol']) == record['method']['solver_tolerance'] and
                              record['solver_info']['presolve'] == 'choose', 'Effective HiGHS settings differ')
                # Values are independently verified above; use recorded reaction order
                # to reproduce the sparse arithmetic, alongside order-robust fsum.
                audit = audit_vector(case['prepared_signature'], record, recipe)
                audits.append(audit)
                if audit:
                    finite_vectors += 1; flux_count += audit['n_fluxes']; row_count += audit['n_balance_rows']
                    maxima['sparse_residual'] = max(maxima['sparse_residual'], audit['certificate']['max_abs_S_residual'])
                    maxima['fsum_residual'] = max(maxima['fsum_residual'], audit['max_fsum_residual'])
                    maxima['bound_violation'] = max(maxima['bound_violation'], audit['certificate']['max_bound_violation'])
                    for metric, value in [('sparse_residual', audit['certificate']['max_abs_S_residual']),
                                          ('fsum_residual', audit['max_fsum_residual']),
                                          ('bound_violation', audit['certificate']['max_bound_violation'])]:
                        per_method_maxima[method][metric] = max(per_method_maxima[method][metric], value)
                methods[method]['statuses'][str(record.get('status', 'missing'))] += 1
            assessment = case_assessment(case['solves'], audits, recipe)
            V.require(case['assessment'] == assessment, 'Independent case acceptance differs')
            unanimous, conflicts = infeasibility_review(case['solves'], audits, recipe)
            if unanimous:
                unanimous_infeasible.append(condition.key)
            status_witness_conflicts.extend({'condition': condition.key, **conflict} for conflict in conflicts)
            outcomes.append({'condition': condition.key, **assessment})
            finite_objectives = [r['objective'] for r in case['solves'] if str(r.get('status', '')).lower() == 'optimal' and V.finite(r.get('objective'))]
            if finite_objectives:
                maxima['objective_spread'] = max(maxima['objective_spread'], max(finite_objectives) - min(finite_objectives))
            for method, row in assessment['methods'].items():
                methods[method]['original'] += int(row['original_primal_gate_recomputed'])
                methods[method]['accepted'] += int(row['accepted_for_comparison'])
                if not row['accepted_for_comparison']:
                    rejected_records.append({'condition': condition.key, 'method': method,
                                             'failed_requirements': row['failed_requirements']})
    V.require(len(outcomes) == 43 and set((run / 'raw_methods').rglob('*.json.gz')) == expected_raw, 'Missing or extra raw case inventory')
    counts = {method: {'status_counts': dict(row['statuses']), 'original_primal_gate_passes': row['original'],
                       'accepted_for_comparison': row['accepted'], 'not_accepted_for_comparison': 43 - row['accepted']}
              for method, row in methods.items()}
    cross = sum(row['accepted_cross_solver_evidence'] for row in outcomes)
    all_three = sum(row['all_three_methods_accepted'] for row in outcomes)
    V.require(summary['n_conditions'] == 43 and summary['n_solves'] == 129 and
              summary['n_complete_finite_vectors'] == finite_vectors and summary['methods'] == counts and
              summary['conditions'] == outcomes and summary['n_conditions_with_cross_solver_evidence'] == cross and
              summary['n_conditions_all_three_methods_accepted'] == all_three, 'Summary counts or decisions differ')
    V.require(summary['input_verification_after']['valid'], 'Comparison did not verify final inputs')
    V.require(frozen == V.verify_manifest(ROOT, manifest), 'Frozen inputs changed during independent verification')
    V.require(direct == {path: V.sha(V.resolve_input(ROOT, path)) for path in direct}, 'Artifacts changed during independent verification')
    return {'status': 'passed_with_status_witness_conflicts' if status_witness_conflicts else 'passed',
            'cases_verified': 43, 'method_records_verified': 129,
            'complete_finite_vectors_verified': finite_vectors, 'reaction_flux_values_verified': flux_count,
            'mass_balance_rows_verified': row_count, 'boundary_records_verified': boundary_count,
            'raw_files_verified': len(expected_raw), 'frozen_input_files_verified': len(frozen),
            'methods': counts, 'conditions_with_cross_solver_evidence': cross,
            'conditions_all_three_methods_accepted': all_three, 'maxima': maxima,
            'per_method_maxima': per_method_maxima, 'rejected_method_records': rejected_records,
            'unanimous_infeasible_conditions': unanimous_infeasible,
            'n_unanimous_infeasible_conditions': len(unanimous_infeasible),
            'infeasible_status_with_tolerance_feasible_witness': status_witness_conflicts,
            'frozen_input_sha256': frozen, 'artifact_sha256': direct, 'study_fingerprint': manifest['content_fingerprint'],
            'scope': 'Independent artifact reconstruction and numerical-decision verification; no optimizer called. Unanimous infeasible statuses lack independent ray certificates and are separate from accepted optimal evidence. Stale or zero vectors are not growth evidence. A passed audit does not make rejected vectors acceptable, prove dual optimality, validate biology, or rewrite failed primary studies.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    run = args.run.resolve()
    V.require(run.is_relative_to(ROOT), 'Numerical artifacts must be in the repository')
    out = args.out or run / 'independent_verification'
    out.mkdir(parents=True, exist_ok=False)
    sources = {str(path.relative_to(ROOT)): V.sha(path) for path in (Path(__file__).resolve(), Path(V.__file__).resolve())}
    for path, name in [(Path(__file__).resolve(), 'checker_source.py'), (Path(V.__file__).resolve(), 'oracle_source.py')]:
        (out / name).write_bytes(path.read_bytes())
    (out / 'attempt.json').write_text(json.dumps({'started_utc': datetime.now(timezone.utc).isoformat(), 'source_sha256': sources,
        'scope': 'Independent checker sources snapshotted before reading numerical results.'}, indent=2) + '\n')
    started = time.time()
    logging.getLogger('cobra').setLevel(logging.ERROR)
    warnings.filterwarnings('ignore', category=FutureWarning, module=r'cobra\.medium\.boundary_types')
    def forbidden(*args, **kwargs):
        raise AssertionError('Independent checker must not optimize')
    V.cobra.Model.optimize = forbidden; V.cobra.Model.slim_optimize = forbidden
    try:
        result = verify(run)
        V.require(sources == {path: V.sha(V.resolve_input(ROOT, path)) for path in sources}, 'Checker source changed during verification')
        result.update(source_sha256=sources, elapsed_seconds=time.time() - started)
        (out / 'result.json').write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
        print(json.dumps({key: result[key] for key in ('status', 'cases_verified', 'method_records_verified',
            'complete_finite_vectors_verified', 'conditions_with_cross_solver_evidence', 'maxima')}))
    except Exception as error:
        (out / 'failure.json').write_text(json.dumps({'type': type(error).__name__, 'message': str(error),
            'elapsed_seconds': time.time() - started, 'scope': 'Checker attempt retained; primary outputs unchanged.'}, indent=2) + '\n')
        raise


if __name__ == '__main__':
    main()
