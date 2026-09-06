"""Fixed post-failure three-method comparison on all 43 curated conditions.

Every returned vector is persisted before LP-integrity or acceptance checks.
Numerical failures remain results; model algebra and acceptance rules never adapt.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
from itertools import combinations
import json
import math
from pathlib import Path
import sys
import time

import highspy
import swiglpk

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import diagnose_medium_numerics as D
from scripts import run_medium_preparation_v2 as R
from gembench.study import freeze_study, verify_study

PLAN = ROOT / 'data/studies/medium_curated_numerics_v1.json'
METHOD_IDS = ('glpk_native', 'highs_primal_simplex', 'highs_simplex_tighter')


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def close(a, b, tolerance):
    return finite(a) and finite(b) and math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)


def assess_record(record, recipe, reaction_ids):
    tolerance = recipe['residual_tolerance']
    certificate = record.get('primal_certificate')
    if not isinstance(certificate, dict):
        certificate = {}
    original = all(finite(certificate.get(key)) and 0 <= certificate[key] <= tolerance
                   for key in ('max_abs_S_residual', 'max_bound_violation'))
    fluxes = record.get('full_primal_fluxes')
    full_vector = isinstance(fluxes, dict) and set(fluxes) == set(reaction_ids) and all(finite(v) for v in fluxes.values())
    compensated = finite(record.get('max_fsum_residual')) and 0 <= record['max_fsum_residual'] <= tolerance
    objective_consistent = close(record.get('objective'), record.get('recomputed_objective'), recipe['objective_tolerance'])
    optimal = str(record.get('status', '')).lower() == 'optimal'
    exception_free = record.get('exception') is None
    reasons = [label for label, passed in [('no_capture_exception', exception_free), ('optimal_status', optimal), ('complete_finite_vector', full_vector),
              ('original_primal_gate', original), ('compensated_residual_gate', compensated),
              ('objective_matches_vector', objective_consistent)] if not passed]
    return {'original_primal_gate_recomputed': original, 'complete_finite_vector': full_vector,
            'compensated_residual_gate': compensated, 'objective_matches_vector': objective_consistent,
            'accepted_for_comparison': not reasons, 'failed_requirements': reasons}


def assess_case(records, recipe, reaction_ids):
    D.require([record['method']['id'] for record in records] == list(METHOD_IDS), 'Method inventory differs')
    assessments = {record['method']['id']: assess_record(record, recipe, reaction_ids) for record in records}
    optimal = [record for record in records if str(record.get('status', '')).lower() == 'optimal']
    finite_optimal = all(finite(record.get('objective')) for record in optimal)
    pairwise = [{'first': a['method']['id'], 'second': b['method']['id'],
                 'first_objective': a.get('objective'), 'second_objective': b.get('objective'),
                 'agree': close(a.get('objective'), b.get('objective'), recipe['objective_tolerance'])}
                for a, b in combinations(optimal, 2)]
    agreement = len(optimal) >= 2 and finite_optimal and all(pair['agree'] for pair in pairwise)
    statuses = {str(record.get('status', '')).lower() for record in records}
    feasibility_conflict = 'optimal' in statuses and any('infeasible' in status for status in statuses)
    status_conflict = 'optimal' in statuses and any('infeasible' in status or 'unbounded' in status for status in statuses)
    glpk_and_highs = assessments['glpk_native']['accepted_for_comparison'] and any(
        assessments[key]['accepted_for_comparison'] for key in METHOD_IDS[1:])
    return {'methods': assessments, 'n_reported_optimal': len(optimal), 'all_optimal_objectives_finite': finite_optimal,
            'pairwise_reported_optimal_objectives': pairwise, 'all_reported_optimal_objectives_agree': agreement,
            'reported_feasibility_conflict': feasibility_conflict, 'reported_status_conflict': status_conflict,
            'glpk_and_at_least_one_highs_accepted': glpk_and_highs,
            'accepted_cross_solver_evidence': glpk_and_highs and agreement and not status_conflict,
            'all_three_methods_accepted': all(item['accepted_for_comparison'] for item in assessments.values())
                                          and agreement and not status_conflict}


def write_gzip(path, value):
    with gzip.open(path, 'xt') as stream:
        json.dump(D.safe(value), stream, allow_nan=False)
        stream.write('\n')


def capture_and_store(model, method, recipe, signature, path):
    settings = {'solver_tolerance': method['solver_tolerance'], 'residual_tolerance': recipe['residual_tolerance'],
                'objective_tolerance': recipe['objective_tolerance']}
    try:
        record = D.run_method(model, method, settings)
    except Exception as error:
        record = {'method': method, 'status': 'exception', 'objective': None, 'full_primal_fluxes': None,
                  'exception': {'type': type(error).__name__, 'message': str(error)}}
    record['algebra_and_chemistry_sha256'] = R.digest(signature)
    write_gzip(path, record)  # Always precedes integrity and acceptance checks.
    D.require(R.signature(model) == signature, 'Solver changed the declared LP; raw result retained')
    return record


def supported_options(methods):
    records = {}
    for method in methods:
        if method['solver'] != 'highs':
            continue
        engine = highspy.Highs()
        options = {'solver': method['method'], 'presolve': 'choose', 'threads': method['threads'],
                   'time_limit': method['time_limit_seconds'], 'primal_feasibility_tolerance': method['solver_tolerance'],
                   'dual_feasibility_tolerance': method['solver_tolerance'], **method['extra_options']}
        effective = {}
        for key, value in options.items():
            status = engine.setOptionValue(key, value)
            get_status, actual = engine.getOptionValue(key)
            D.require(status == get_status == highspy.HighsStatus.kOk and actual == value, f'Unsupported numerical option: {key}')
            effective[key] = actual
        for key in ('simplex_strategy', 'simplex_scale_strategy'):
            status, value = engine.getOptionValue(key)
            D.require(status == highspy.HighsStatus.kOk, f'Missing numerical default: {key}')
            effective[key] = value
        records[method['id']] = effective
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    try:
        recipe = R.read(PLAN)
        D.require(recipe['condition_count'] == 43 and recipe['solve_count'] == 129, 'Wrong declared panel size')
        D.require([method['id'] for method in recipe['methods']] == list(METHOD_IDS), 'Wrong declared method inventory')
        parent = R.read(ROOT / recipe['parent_manifest'])
        D.require(verify_study(ROOT, parent)['valid'], 'Failed primary-study inputs changed')
        inventory = R.inventory_module()
        evidence = R.read(ROOT / recipe['inventory_record'])
        declared = next(item for item in evidence['models'] if item['label'] == recipe['model_label'])
        conditions = inventory.metadata_conditions('Putida')
        actual = [{'key': c.key, 'name': c.name, 'media': c.media, 'exchanges': c.exchanges} for c in conditions]
        D.require(actual == declared['conditions'] and len(actual) == 43 and len({c['key'] for c in actual}) == 43,
                  'Condition inventory differs from frozen source metadata')
        effective = supported_options(recipe['methods'])
        exposure = R.read(ROOT / 'data/studies/exposure_registry_v1.json')
        paths = {record['path'] for record in parent['files']} | set(recipe['input_paths']) | set(recipe['diagnostic_inputs'])
        manifest = freeze_study(ROOT, {'schema_version': 1, 'study_id': recipe['study_id'], 'evaluation_role': 'development',
            'development_organisms': exposure['development_organisms'], 'quarantined_organisms': exposure['quarantined_organisms'],
            'paths': sorted(paths), 'recipe': recipe, 'effective_highs_options': effective,
            'runtime': {**R.software_versions(), 'glpk': swiglpk.glp_version(), 'highs': highspy.Highs().version()}})
        D.require(manifest['repository']['available'] and manifest['repository']['head'] and not manifest['repository']['dirty'],
                  'Commit the complete proposal and inputs before running from a clean checkout')
        R.write(out / 'manifest.json', manifest)
        base, provenance = inventory.build_parent(recipe['model_label'])
        policy = R.CompartmentPolicy(**recipe['compartments'])
        media_names = sorted({c.media for c in conditions})
        completion = [R.base_medium(name) for name in media_names]
        legacy = base.copy()
        R.complete_medium_transport(legacy, media_names, list(policy.completion_exclude))
        R.write(out / 'source_context.json', {'provenance': provenance, 'conditions': actual, 'effective_highs_options': effective})
        outcomes, statuses = [], {method: Counter() for method in METHOD_IDS}
        accepted = Counter()
        original_passes = Counter()
        finite_vectors = 0
        with gzip.open(out / 'cases.jsonl.gz', 'xt') as stream:
            for index, condition in enumerate(conditions):
                medium = R.medium_for(condition)
                prepared = R.prepare_medium(base, medium, policy=policy, completion_media=completion, missing_policy='report')
                old = legacy.copy()
                R.apply_medium(old, medium, close_all=True)
                signature = R.signature(prepared.model)
                D.require(signature == R.signature(old), 'Fresh legacy/v2 LP or chemical metadata differs')
                directory = out / 'raw_methods' / f'{index:02d}'
                directory.mkdir(parents=True)
                write_gzip(directory / 'prepared_signature.json.gz', signature)
                models = [prepared.model] + [prepared.model.copy() for _ in recipe['methods'][1:]]
                records = []
                for method, model in zip(recipe['methods'], models):
                    D.require(R.signature(model) == signature, 'Method starts with a different LP')
                    records.append(capture_and_store(model, method, recipe, signature, directory / f'{method["id"]}.json.gz'))
                # Raw method vectors have already been saved, including failures.
                assessment = assess_case(records, recipe, signature['reactions'])
                row = {'index': index, 'model': recipe['model_label'], 'condition': condition.key,
                       'medium_report': prepared.report, 'prepared_signature': signature,
                       'algebra_and_chemistry_sha256': R.digest(signature), 'legacy_signature_sha256': R.digest(R.signature(old)),
                       'solves': records, 'assessment': assessment}
                stream.write(json.dumps(D.safe(row), allow_nan=False) + '\n')
                stream.flush()
                for record in records:
                    key = record['method']['id']
                    statuses[key][str(record.get('status', 'missing'))] += 1
                    accepted[key] += int(assessment['methods'][key]['accepted_for_comparison'])
                    original_passes[key] += int(assessment['methods'][key]['original_primal_gate_recomputed'])
                    finite_vectors += int(assessment['methods'][key]['complete_finite_vector'])
                outcomes.append({'condition': condition.key, **assessment})
                print('Completed:', index + 1, '/', len(conditions), condition.key,
                      'cross-solver evidence:', assessment['accepted_cross_solver_evidence'], flush=True)
        D.require(len(outcomes) == recipe['condition_count'], 'Incomplete numerical panel')
        verified = verify_study(ROOT, manifest)
        D.require(verified['valid'], 'Frozen inputs changed during comparison')
        R.write(out / 'summary.json', {'completed': True, 'study': recipe, 'study_fingerprint': manifest['content_fingerprint'],
            'n_conditions': len(outcomes), 'n_solves': len(outcomes) * len(METHOD_IDS), 'n_complete_finite_vectors': finite_vectors,
            'methods': {key: {'status_counts': dict(statuses[key]), 'original_primal_gate_passes': original_passes[key],
                             'accepted_for_comparison': accepted[key], 'not_accepted_for_comparison': len(outcomes) - accepted[key]}
                        for key in METHOD_IDS}, 'conditions': outcomes,
            'n_conditions_with_cross_solver_evidence': sum(item['accepted_cross_solver_evidence'] for item in outcomes),
            'n_conditions_all_three_methods_accepted': sum(item['all_three_methods_accepted'] for item in outcomes),
            'input_verification_after': verified, 'seconds': time.time() - started,
            'scope': 'Fixed post-failure numerical comparison. Primal feasibility and cross-method objective agreement are not a mathematical optimality proof, biochemical validation, phenotype rescore or rewrite of failed primary runs.'})
    except Exception as error:
        R.write(out / 'failure.json', {'type': type(error).__name__, 'message': str(error), 'seconds': time.time() - started})
        raise


if __name__ == '__main__':
    main()
