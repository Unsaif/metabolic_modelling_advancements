"""Four predeclared HiGHS follow-ups on the exact saved one-condition LP.

The first diagnostic is immutable. This attempt changes solver options only,
including one tighter solver tolerance; the acceptance gate remains 1e-8.
All supported/effective option values and input hashes precede any solve.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
from importlib.metadata import version
import json
import math
from pathlib import Path
import platform
import shutil
import sys
import time

import cobra
import highspy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import diagnose_medium_numerics as D
from scripts import run_medium_preparation_v2 as R
from gembench.study import verify_study


METHODS = [
    {'id': 'highs_primal_simplex', 'tolerance': 1e-9, 'presolve': 'choose', 'extra_options': {'simplex_strategy': 4}},
    {'id': 'highs_dual_no_presolve', 'tolerance': 1e-9, 'presolve': 'off', 'extra_options': {'simplex_strategy': 1}},
    {'id': 'highs_simplex_tighter', 'tolerance': 1e-10, 'presolve': 'choose', 'extra_options': {}},
    {'id': 'highs_simplex_no_scaling', 'tolerance': 1e-9, 'presolve': 'choose', 'extra_options': {'simplex_scale_strategy': 0}},
]


def validate_options(method):
    engine = highspy.Highs()
    options = {'output_flag': False, 'solver': 'simplex', 'presolve': method['presolve'], 'parallel': 'off',
               'primal_feasibility_tolerance': method['tolerance'], 'dual_feasibility_tolerance': method['tolerance'],
               'time_limit': 60., 'threads': 0, **method['extra_options']}
    records = {}
    for key, value in options.items():
        status, kind = engine.getOptionType(key)
        D.require(status == highspy.HighsStatus.kOk, f'Unsupported option: {key}')
        set_status = engine.setOptionValue(key, value)
        get_status, effective = engine.getOptionValue(key)
        D.require(set_status == get_status == highspy.HighsStatus.kOk and effective == value,
                  f'Option not applied as requested: {key}')
        records[key] = {'type': str(kind), 'requested': value, 'effective': effective}
    for key in ('simplex_strategy', 'simplex_scale_strategy', 'primal_residual_tolerance', 'dual_residual_tolerance'):
        if key not in records:
            status, value = engine.getOptionValue(key)
            D.require(status == highspy.HighsStatus.kOk, f'Missing relevant default: {key}')
            records[key] = {'effective': value, 'source': 'unchanged default'}
    return records


def solve(model, method, acceptance):
    R.assert_ordinary_lp(model)
    matrix = R.matrix_model(model)
    started = time.time()
    record = {'method': method, 'started_at': datetime.now(timezone.utc).isoformat()}
    try:
        answer = R.W.solve_highs(matrix, method='simplex', presolve=method['presolve'],
            threads=0, feas_tol=method['tolerance'], opt_tol=method['tolerance'], time_limit=60.,
            extra_options={'parallel': 'off', **method['extra_options']})
        record.update(status=answer.status, objective=answer.objective, solver_info=answer.info)
        if answer.x is not None:
            record['full_primal_fluxes'] = {str(rid): float(value) for rid, value in zip(matrix.rxns, answer.x)}
            record.update(D.vector_audit(model, matrix, answer.x, acceptance))
        else:
            record.update(full_primal_fluxes=None, finite_primal_vector=False, passes_original_primal_gate=False)
        record['passes_original_case_gate'] = (answer.status.lower() == 'optimal' and math.isfinite(answer.objective)
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
    shutil.copyfile(Path(D.__file__).resolve(), out / 'attempt01_helpers_source.py')
    try:
        parent = ROOT / 'results/medium_preparation_2026_09_06/numerical_diagnostic'
        prior_inputs = R.read(parent / 'input_hashes.json')
        prior_artifacts = R.read(parent / 'presolve_artifact_hashes.json')
        D.unchanged(prior_inputs)
        D.unchanged(prior_artifacts)
        manifest = R.read(ROOT / 'results/medium_preparation_2026_09_06/runs/main_v2/manifest.json')
        D.require(verify_study(ROOT, manifest)['valid'], 'Frozen study changed')
        inputs = [*prior_inputs, D.file_record(source)]
        inputs.extend(D.file_record(path) for path in sorted(parent.iterdir()) if path.is_file())
        D.write(out / 'input_hashes.json', inputs)
        first_settings = R.read(parent / 'presolve_settings.json')
        settings = {'scope': 'Second post-failure numerical diagnosis, limited to four specified options on the exact saved LP.',
            'condition': first_settings['condition'], 'model_label': first_settings['model_label'],
            'acceptance_residual_tolerance': first_settings['residual_tolerance'],
            'objective_tolerance': first_settings['objective_tolerance'], 'methods': METHODS,
            'option_documentation': 'https://ergo-code.github.io/HiGHS/dev/options/definitions/',
            'no_changes': ['stoichiometry', 'reaction bounds', 'objective coefficients', 'GPRs', 'medium', '1e-8 acceptance gate'],
            'changes': 'Solver algorithm options only; one declared arm tightens solver tolerances from1e-9 to1e-10.'}
        D.require(settings['acceptance_residual_tolerance'] == 1e-8, 'Acceptance gate changed')
        D.write(out / 'presolve_settings.json', settings)
        D.write(out / 'supported_effective_options.json', {method['id']: validate_options(method) for method in METHODS})
        D.write(out / 'runtime.json', {'created_at': datetime.now(timezone.utc).isoformat(), 'python': sys.version,
            'platform': platform.platform(), 'highs': highspy.Highs().version(),
            'versions': {name: version(name) for name in ('cobra', 'numpy', 'scipy', 'highspy')},
            'script_sha256': D.sha(source), 'attempt01_helper_sha256': D.sha(D.__file__)})
        model = cobra.io.read_sbml_model(str(parent / 'prepared_model.xml.gz'))
        with gzip.open(parent / 'prepared_signature.json.gz', 'rt') as stream:
            signature = json.load(stream)
        D.require(R.signature(model) == signature, 'Saved SBML does not reproduce exact original LP and chemical metadata')
        D.write(out / 'prepared_context.json', {'source_model': D.file_record(parent / 'prepared_model.xml.gz'),
            'algebra_and_chemistry_sha256': R.digest(signature), 'n_reactions': len(model.reactions),
            'n_metabolites': len(model.metabolites), 'condition': settings['condition'],
            'exact_signature_matches_attempt01': True})
        artifacts = [D.file_record(path) for path in sorted(out.iterdir()) if path.is_file()]
        D.write(out / 'presolve_artifact_hashes.json', artifacts)
        D.unchanged(inputs)
        D.unchanged(artifacts)
        results = []
        for method in METHODS:
            result = solve(model, method, settings['acceptance_residual_tolerance'])
            D.require(R.signature(model) == signature, 'Solver method changed the LP')
            D.write(out / f'{method["id"]}.json', result)
            results.append(result)
            print(json.dumps(D.safe({'method': method['id'], 'status': result.get('status'), 'objective': result.get('objective'),
                'certificate': result.get('primal_certificate'), 'passed': result['passes_original_case_gate']})), flush=True)
        D.unchanged(inputs)
        D.unchanged(artifacts)
        after = verify_study(ROOT, manifest)
        D.require(after['valid'], 'Frozen study input changed')
        native = R.read(parent / 'glpk_native.json')['objective']
        compact = [{key: result.get(key) for key in ('method', 'status', 'objective', 'seconds', 'primal_certificate',
            'max_fsum_residual', 'n_sparse_rows_above_original_tolerance', 'passes_original_case_gate', 'exception')}
                   for result in results]
        for record in compact:
            record['absolute_difference_from_prior_native_glpk_objective'] = abs(record['objective'] - native) \
                if isinstance(record['objective'], (int, float)) and math.isfinite(record['objective']) else None
        D.write(out / 'summary.json', {'settings': settings, 'methods': compact,
            'input_verification_after': after, 'source_and_pre_solve_artifacts_unchanged': True,
            'scope': 'One-condition post-failure method diagnostics. No original failures overwritten, no all-panel conclusion, no changed biochemical or acceptance assumptions.'})
    except Exception as error:
        D.write(out / 'failure.json', {'type': type(error).__name__, 'message': str(error)})
        raise


if __name__ == '__main__':
    main()
