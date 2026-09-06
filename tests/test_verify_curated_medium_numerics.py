"""Independent numerical-checker fixtures; no solver or phenotype data."""
from copy import deepcopy
import importlib.util
from pathlib import Path
import sys

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
PATH = REPO / 'scripts/verify_curated_medium_numerics.py'
SPEC = importlib.util.spec_from_file_location('independent_curated_checker', PATH)
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)
RECIPE = {'residual_tolerance': 1e-8, 'objective_tolerance': 1e-8}


def fixture(growth=5.):
    structure = {'reactions': {'UPT': {'stoichiometry': {'a_c': 1.}, 'bounds': [0., 10.], 'gpr': ''},
                               'Growth': {'stoichiometry': {'a_c': -1.}, 'bounds': [0., 1000.], 'gpr': 'g'}},
                 'metabolites': {'a_c': {'formula': 'C', 'charge': 0}}, 'objective': {'Growth': 1.}, 'direction': 'max'}
    residue = 5. - growth
    terms = [{'reaction': 'UPT', 'coefficient': 1., 'flux': 5., 'term': 5.},
             {'reaction': 'Growth', 'coefficient': -1., 'flux': growth, 'term': -growth}]
    bounds = [{'reaction': 'UPT', 'flux': 5., 'bounds': [0., 10.], 'violation': 0.},
              {'reaction': 'Growth', 'flux': growth, 'bounds': [0., 1000.], 'violation': 0.}]
    gate = abs(residue) <= RECIPE['residual_tolerance']
    record = {'method': {'id': 'glpk_native'}, 'status': 'optimal', 'objective': growth,
        'full_primal_fluxes': {'UPT': 5., 'Growth': growth}, 'finite_primal_vector': True,
        'complete_metabolite_residuals': {'a_c': {'sparse': residue, 'fsum': residue}},
        'primal_certificate': {'max_abs_S_residual': abs(residue), 'max_bound_violation': 0.,
                              'n_S_rows_over_1e-6': int(abs(residue) > 1e-6), 'n_bounds_over_1e-6': 0},
        'max_fsum_residual': abs(residue), 'recomputed_objective': growth,
        'n_sparse_rows_above_original_tolerance': int(not gate), 'n_fsum_rows_above_original_tolerance': int(not gate),
        'worst_rows': [{'metabolite': 'a_c', 'sparse_residual': residue, 'fsum_residual': residue,
                        'sum_abs_terms': 5. + growth, 'n_terms': 2,
                        'largest_terms': sorted(terms, key=lambda row: abs(row['term']), reverse=True)}],
        'worst_bounds': bounds, 'largest_absolute_fluxes': sorted(bounds, key=lambda row: abs(row['flux']), reverse=True),
        'passes_original_primal_gate': gate, 'passes_original_case_gate': gate}
    return structure, record


def test_finite_feasible_vector_is_independently_reconstructed():
    structure, record = fixture()
    audit = C.audit_vector(structure, record, RECIPE)
    assert audit['n_fluxes'] == 2 and audit['n_balance_rows'] == 1 and audit['objective'] == 5
    assert C.method_assessment(record, audit, RECIPE)['accepted_for_comparison']


def test_above_tolerance_vector_is_audited_and_explicitly_rejected():
    structure, record = fixture(5.0001)
    audit = C.audit_vector(structure, record, RECIPE)
    decision = C.method_assessment(record, audit, RECIPE)
    assert audit['max_fsum_residual'] > 1e-8
    assert not decision['accepted_for_comparison']
    assert decision['failed_requirements'] == ['original_primal_gate', 'compensated_residual_gate']


@pytest.mark.parametrize('kind', ['residual', 'certificate', 'term', 'bound', 'objective', 'count', 'flag'])
def test_tampered_vector_diagnostics_fail_verification(kind):
    structure, record = fixture()
    if kind == 'residual':
        record['complete_metabolite_residuals']['a_c']['fsum'] = 1
    elif kind == 'certificate':
        record['primal_certificate']['max_bound_violation'] = 1
    elif kind == 'term':
        record['worst_rows'][0]['largest_terms'][0]['coefficient'] = 2
    elif kind == 'bound':
        record['worst_bounds'][0]['bounds'][1] = 100
    elif kind == 'objective':
        record['recomputed_objective'] = 0
    elif kind == 'count':
        record['n_fsum_rows_above_original_tolerance'] = 1
    else:
        record['passes_original_primal_gate'] = False
    with pytest.raises(ValueError):
        C.audit_vector(structure, record, RECIPE)


@pytest.mark.parametrize('fluxes', [None, {'UPT': 5.}, {'UPT': 'nan', 'Growth': 5.}])
def test_missing_or_nonfinite_vector_cannot_be_accepted(fluxes):
    structure, record = fixture()
    record = {'method': record['method'], 'status': 'time limit', 'objective': None,
              'full_primal_fluxes': fluxes, 'finite_primal_vector': False}
    assert C.audit_vector(structure, record, RECIPE) is None
    assert not C.method_assessment(record, None, RECIPE)['accepted_for_comparison']


def test_ranked_diagnostics_allow_ties_but_not_omitted_larger_terms():
    candidates = {'a': {'id': 'a', 'weight': 1}, 'b': {'id': 'b', 'weight': 1}, 'c': {'id': 'c', 'weight': 0}}
    C.ranked_subset([candidates['b']], candidates, 'id', lambda row: row['weight'], 1, 'toy')
    with pytest.raises(ValueError):
        C.ranked_subset([candidates['c']], candidates, 'id', lambda row: row['weight'], 1, 'toy')


def test_case_requires_native_and_highs_and_all_optimal_objectives_agree():
    structure, record = fixture()
    audit = C.audit_vector(structure, record, RECIPE)
    records, audits = [], []
    for key in C.METHODS:
        item = deepcopy(record); item['method']['id'] = key
        records.append(item); audits.append(deepcopy(audit))
    audits[2]['certificate']['max_abs_S_residual'] = 1e-6
    result = C.case_assessment(records, audits, RECIPE)
    assert result['accepted_cross_solver_evidence'] and not result['all_three_methods_accepted']
    records[2]['objective'] = 6
    result = C.case_assessment(records, audits, RECIPE)
    assert not result['accepted_cross_solver_evidence'] and not result['all_reported_optimal_objectives_agree']


@pytest.mark.parametrize('status', ['infeasible', 'unbounded', 'Infeasible or Unbounded'])
def test_contradictory_status_blocks_otherwise_agreeing_methods(status):
    structure, record = fixture()
    audit = C.audit_vector(structure, record, RECIPE)
    records = [deepcopy(record) for _ in C.METHODS]
    for key, row in zip(C.METHODS, records):
        row['method']['id'] = key
    records[2] = {'method': records[2]['method'], 'status': status, 'objective': None, 'full_primal_fluxes': None}
    result = C.case_assessment(records, [audit, audit, None], RECIPE)
    assert result['reported_status_conflict'] and not result['accepted_cross_solver_evidence']


def test_capture_exception_never_accepted_despite_finite_optimum():
    structure, record = fixture()
    audit = C.audit_vector(structure, record, RECIPE)
    record['exception'] = {'type': 'ValueError', 'message': 'capture failed'}
    result = C.method_assessment(record, audit, RECIPE)
    assert not result['accepted_for_comparison']
    assert result['failed_requirements'] == ['no_capture_exception']


def test_infeasible_status_and_tolerance_feasible_witness_are_flagged_separately():
    structure, record = fixture()
    audit = C.audit_vector(structure, record, RECIPE)
    record['status'] = 'infeasible'
    unanimous, conflicts = C.infeasibility_review([record], [audit], RECIPE)
    assert unanimous and len(conflicts) == 1
    assert not C.method_assessment(record, audit, RECIPE)['accepted_for_comparison']
    audit['certificate']['max_bound_violation'] = 8.39
    unanimous, conflicts = C.infeasibility_review([record], [audit], RECIPE)
    assert unanimous and conflicts == []
    record['status'] = 'Infeasible or Unbounded'
    audit['certificate']['max_bound_violation'] = 0.
    unanimous, conflicts = C.infeasibility_review([record], [audit], RECIPE)
    assert not unanimous and conflicts == []


def test_checker_snapshots_and_refuses_overwriting_an_attempt(tmp_path, monkeypatch):
    out = tmp_path / 'attempt'
    run = REPO / 'results/medium_preparation_2026_09_06/curated_numerical_comparison'
    monkeypatch.setattr(C.sys, 'argv', ['verify', '--run', str(run), '--out', str(out)])
    def fail(path):
        assert (out / 'checker_source.py').read_bytes() == Path(C.__file__).read_bytes()
        assert (out / 'oracle_source.py').read_bytes() == Path(C.V.__file__).read_bytes()
        raise ValueError('Toy verification failure')
    monkeypatch.setattr(C, 'verify', fail)
    # Restore the explicit no-optimization guards after testing the CLI.
    monkeypatch.setattr(C.V.cobra.Model, 'optimize', C.V.cobra.Model.optimize)
    monkeypatch.setattr(C.V.cobra.Model, 'slim_optimize', C.V.cobra.Model.slim_optimize)
    with pytest.raises(ValueError, match='Toy verification failure'):
        C.main()
    before = {p.name: p.read_bytes() for p in out.iterdir()}
    with pytest.raises(FileExistsError):
        C.main()
    assert before == {p.name: p.read_bytes() for p in out.iterdir()}
