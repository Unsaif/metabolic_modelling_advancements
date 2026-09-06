"""One-condition numerical diagnosis, preserving all attempts and primal vectors.

No LP coefficients, bounds or acceptance tolerances are changed. This is a
post-failure solver-method diagnosis, not a replacement for the failed study.
All four methods, sources and input hashes are recorded before any solve.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
from importlib.metadata import version
import inspect
import json
import math
from pathlib import Path
import platform
import shutil
import sys
import time

import cobra
import highspy
import numpy as np
from optlang import glpk_interface
import swiglpk

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import run_medium_preparation_v2 as R
from gembench.study import verify_study


METHODS = [
    {'id': 'glpk_native', 'solver': 'glpk', 'basis': 'native',
     'scaling': 'GLP_SF_AUTO already performed by optlang.Model._optimize'},
    {'id': 'glpk_advanced_basis', 'solver': 'glpk', 'basis': 'glp_adv_basis(problem,0)',
     'scaling': 'GLP_SF_AUTO already performed by optlang.Model._optimize'},
    {'id': 'highs_simplex', 'solver': 'highs', 'method': 'simplex',
     'extra_options': {'parallel': 'off'}, 'time_limit_seconds': 60., 'threads': 0},
    {'id': 'highs_ipm', 'solver': 'highs', 'method': 'ipm',
     'extra_options': {'parallel': 'off', 'ipm_optimality_tolerance': 1e-9, 'run_crossover': 'on'},
     'time_limit_seconds': 60., 'threads': 0},
]


def require(test, message):
    if not test:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def safe(value):
    if isinstance(value, np.generic):
        return safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return repr(value)
    if isinstance(value, dict):
        return {str(key): safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe(item) for item in value]
    return value


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(safe(value), stream, indent=2, allow_nan=False)
        stream.write('\n')


def file_record(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': sha(path), 'size_bytes': path.stat().st_size}


def unchanged(records):
    require(all(file_record(record['path']) == record for record in records), 'Diagnostic input changed')


def vector_audit(model, matrix, values, tolerance):
    require(len(values) == len(matrix.rxns), 'Incomplete primal vector')
    fluxes = {str(rid): float(value) for rid, value in zip(matrix.rxns, values)}
    if not all(math.isfinite(value) for value in fluxes.values()):
        return {'full_primal_fluxes': fluxes, 'finite_primal_vector': False, 'passes_original_primal_gate': False}
    sparse = np.asarray(matrix.S @ values).reshape(-1)
    certificate = R.W.certify(matrix, values)
    rows, complete_residuals = [], {}
    for index, metabolite in enumerate(model.metabolites):
        terms = [{'reaction': reaction.id, 'coefficient': float(reaction.metabolites[metabolite]),
                  'flux': fluxes[reaction.id], 'term': float(reaction.metabolites[metabolite]) * fluxes[reaction.id]}
                 for reaction in metabolite.reactions]
        residual = math.fsum(term['term'] for term in terms)
        absolute_sum = math.fsum(abs(term['term']) for term in terms)
        complete_residuals[metabolite.id] = {'sparse': float(sparse[index]), 'fsum': residual}
        rows.append({'metabolite': metabolite.id, 'sparse_residual': float(sparse[index]),
                     'fsum_residual': residual, 'sum_abs_terms': absolute_sum, 'n_terms': len(terms),
                     'largest_terms': sorted(terms, key=lambda term: abs(term['term']), reverse=True)[:20]})
    violations = [{'reaction': reaction.id, 'flux': fluxes[reaction.id], 'bounds': list(reaction.bounds),
                   'violation': max(reaction.lower_bound - fluxes[reaction.id], fluxes[reaction.id] - reaction.upper_bound, 0.)}
                  for reaction in model.reactions]
    objective = math.fsum(float(coefficient) * fluxes[reaction.id]
                         for reaction, coefficient in cobra.util.solver.linear_reaction_coefficients(model).items())
    return {'full_primal_fluxes': fluxes, 'finite_primal_vector': True,
            'primal_certificate': certificate, 'recomputed_objective': objective,
            'complete_metabolite_residuals': complete_residuals,
            'max_fsum_residual': max(abs(row['fsum_residual']) for row in rows),
            'n_sparse_rows_above_original_tolerance': sum(abs(row['sparse_residual']) > tolerance for row in rows),
            'n_fsum_rows_above_original_tolerance': sum(abs(row['fsum_residual']) > tolerance for row in rows),
            'worst_rows': sorted(rows, key=lambda row: max(abs(row['sparse_residual']), abs(row['fsum_residual'])), reverse=True)[:20],
            'worst_bounds': sorted(violations, key=lambda row: row['violation'], reverse=True)[:20],
            'largest_absolute_fluxes': sorted(violations, key=lambda row: abs(row['flux']), reverse=True)[:30],
            'passes_original_primal_gate': all(math.isfinite(certificate[key]) and 0 <= certificate[key] <= tolerance
                 for key in ('max_abs_S_residual', 'max_bound_violation'))}


def run_method(model, method, settings):
    R.assert_ordinary_lp(model)
    matrix = R.matrix_model(model)
    record = {'method': method, 'started_at': datetime.now(timezone.utc).isoformat()}
    started = time.time()
    try:
        values = None
        if method['solver'] == 'glpk':
            model.solver = 'glpk'
            model.solver.configuration.tolerances.feasibility = settings['solver_tolerance']
            if method['basis'] != 'native':
                swiglpk.glp_adv_basis(model.solver.problem, 0)
            record['effective_settings'] = {'feasibility_tolerance': model.solver.configuration.tolerances.feasibility,
                'presolve': model.solver.configuration.presolve, 'timeout': model.solver.configuration.timeout,
                'simplex_tol_bnd': model.solver.configuration._smcp.tol_bnd,
                'simplex_tol_dj': model.solver.configuration._smcp.tol_dj}
            answer = model.optimize()
            record.update(status=model.solver.status, objective=answer.objective_value)
            if answer.fluxes is not None:
                values = answer.fluxes.reindex(matrix.rxns).to_numpy()
        else:
            answer = R.W.solve_highs(matrix, method=method['method'], threads=method['threads'],
                feas_tol=settings['solver_tolerance'], opt_tol=settings['solver_tolerance'],
                time_limit=method['time_limit_seconds'], extra_options=method['extra_options'])
            record.update(status=answer.status, objective=answer.objective, solver_info=answer.info)
            values = answer.x
        if values is not None:
            # Preserve the vector even if a subsequent diagnostic calculation
            # raises (for example, overflow while accumulating residuals).
            record['full_primal_fluxes'] = {str(rid): float(value) for rid, value in zip(matrix.rxns, values)}
            record.update(vector_audit(model, matrix, values, settings['residual_tolerance']))
        else:
            record.update(full_primal_fluxes=None, finite_primal_vector=False, passes_original_primal_gate=False)
        record['passes_original_case_gate'] = (record['status'].lower() == 'optimal'
            and isinstance(record['objective'], (int, float)) and math.isfinite(record['objective'])
            and record['passes_original_primal_gate'])
    except Exception as error:
        record.update(exception={'type': type(error).__name__, 'message': str(error)}, passes_original_case_gate=False)
    record['seconds'] = time.time() - started
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve()
    shutil.copyfile(source, out / 'diagnostic_source.py')
    try:
        old_directory = ROOT / 'results/medium_preparation_2026_09_06/runs/main_v2'
        manifest = R.read(old_directory / 'manifest.json')
        verification = verify_study(ROOT, manifest)
        require(verification['valid'], 'Frozen main_v2 input changed')
        recipe = R.read(R.PLAN)
        require(recipe == manifest['plan']['recipe'], 'Diagnostic recipe differs from failed run')
        runtime_source = Path(inspect.getfile(glpk_interface))
        inputs = [file_record(ROOT / record['path']) for record in manifest['files']]
        inputs += [file_record(old_directory / name) for name in ('manifest.json', 'failure.json', 'cases.jsonl.gz')]
        inputs += [file_record(source), file_record(runtime_source)]
        shutil.copyfile(runtime_source, out / 'optlang_glpk_interface_source.py')
        write(out / 'input_hashes.json', inputs)
        settings = {'study': 'post_failure_single_curated_condition_numerical_diagnosis',
            'model_label': 'Putida_iJN1463', 'condition': 'L-Arginine | MOPS minimal media_noCarbon',
            'solver_tolerance': recipe['solver_tolerance'], 'residual_tolerance': recipe['residual_tolerance'],
            'objective_tolerance': recipe['objective_tolerance'], 'methods': METHODS,
            'scope': 'Exactly one already failed model/condition. Methods predeclared before these solves; no changed coefficients, bounds, medium, GPRs or acceptance thresholds. Source metadata/fitness headers only; no numeric fitness access. Baseline failures remain failures.'}
        write(out / 'presolve_settings.json', settings)
        write(out / 'runtime.json', {'created_at': datetime.now(timezone.utc).isoformat(), 'python': sys.version,
            'platform': platform.platform(), 'versions': {name: version(name) for name in ('cobra', 'optlang', 'numpy', 'scipy', 'highspy', 'swiglpk')},
            'glpk': swiglpk.glp_version(), 'highs': highspy.Highs().version(), 'script_sha256': sha(source),
            'frozen_input_verification': verification,
            'native_glpk_auto_scaling_source': inspect.getsource(glpk_interface.Model._optimize)})
        inventory = R.inventory_module()
        base, provenance = inventory.build_parent(settings['model_label'])
        conditions = inventory.metadata_conditions('Putida')
        condition = conditions[0]
        require(condition.key == settings['condition'], 'First curated condition differs')
        entry = next(item for item in recipe['models'] if item['label'] == settings['model_label'])
        policy = R.CompartmentPolicy(**entry['compartments'])
        medium = R.medium_for(condition)
        media_names = sorted({item.media for item in conditions})
        completion_media = [R.base_medium(name) for name in media_names]
        prepared = R.prepare_medium(base, medium, policy=policy, completion_media=completion_media, missing_policy='report')
        legacy = base.copy()
        R.complete_medium_transport(legacy, media_names, list(policy.completion_exclude))
        R.apply_medium(legacy, medium, close_all=True)
        signature = R.signature(prepared.model)
        require(signature == R.signature(legacy), 'Prepared LP differs from legacy or frozen policy')
        # Match the runner's pre-solve idempotence check, which itself never optimizes.
        again = R.prepare_medium(prepared.model, medium, policy=policy, completion_media=completion_media, missing_policy='report')
        require(R.signature(again.model) == signature, 'Preparation not idempotent')
        cobra.io.write_sbml_model(prepared.model, str(out / 'prepared_model.xml.gz'))
        with gzip.open(out / 'prepared_signature.json.gz', 'xt') as stream:
            json.dump(signature, stream, allow_nan=False)
        write(out / 'prepared_context.json', {'condition': condition.key, 'medium_report': prepared.report,
            'provenance': provenance, 'algebra_and_chemistry_sha256': R.digest(signature),
            'n_reactions': len(prepared.model.reactions), 'n_metabolites': len(prepared.model.metabolites),
            'bound_max_abs': max(abs(bound) for reaction in prepared.model.reactions for bound in reaction.bounds),
            'coefficient_min_abs': min(abs(value) for reaction in prepared.model.reactions for value in reaction.metabolites.values() if value),
            'coefficient_max_abs': max(abs(value) for reaction in prepared.model.reactions for value in reaction.metabolites.values())})
        # Freeze the complete preparation artifacts and method settings before solving.
        artifacts = [file_record(path) for path in sorted(out.iterdir()) if path.is_file()]
        write(out / 'presolve_artifact_hashes.json', artifacts)
        unchanged(inputs)
        unchanged(artifacts)
        models = [prepared.model] + [prepared.model.copy() for _ in METHODS[1:]]
        results = []
        for method, model in zip(METHODS, models):
            require(R.signature(model) == signature, 'Method changed the LP before solving')
            result = run_method(model, method, settings)
            require(R.signature(model) == signature, 'Numerical method changed LP algebra/bounds')
            write(out / f'{method["id"]}.json', result)
            results.append(result)
            print(json.dumps({'method': method['id'], 'status': result.get('status'), 'objective': safe(result.get('objective')),
                'certificate': safe(result.get('primal_certificate')), 'passed': result['passes_original_case_gate']}), flush=True)
        unchanged(inputs)
        unchanged(artifacts)
        verified_after = verify_study(ROOT, manifest)
        require(verified_after['valid'], 'Frozen source changed during diagnosis')
        compact = []
        for result in results:
            compact.append({key: result.get(key) for key in ('method', 'status', 'objective', 'seconds', 'primal_certificate',
                'max_fsum_residual', 'n_sparse_rows_above_original_tolerance', 'n_fsum_rows_above_original_tolerance',
                'passes_original_case_gate', 'exception')})
        objectives = [float(result['objective']) for result in results if result.get('status', '').lower() == 'optimal'
                      and isinstance(result.get('objective'), (int, float)) and math.isfinite(result['objective'])]
        write(out / 'summary.json', {'settings': settings, 'methods': compact,
            'max_optimal_objective_spread': max(objectives) - min(objectives) if objectives else None,
            'input_verification_after': verified_after, 'source_and_pre_solve_artifacts_unchanged': True,
            'scope': 'Post-failure numerical diagnosis for one condition; full primal vectors retained whether accepted or rejected. Solver-method success does not rewrite the failed primary run, demonstrate all-panel reliability, or validate model biology.'})
    except Exception as error:
        write(out / 'failure.json', {'type': type(error).__name__, 'message': str(error)})
        raise


if __name__ == '__main__':
    main()
