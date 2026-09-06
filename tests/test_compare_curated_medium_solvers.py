"""Acceptance and evidence-preservation tests; no real-model optimization."""
from copy import deepcopy
import gzip
import json

import cobra
import pytest

from scripts import compare_curated_medium_solvers as C


RECIPE = {'residual_tolerance': 1e-8, 'objective_tolerance': 1e-8}
REACTIONS = {'EX_a_e', 'BIOMASS_TEST'}


def record(method, objective=1.):
    return {'method': {'id': method}, 'status': 'Optimal', 'objective': objective,
            'recomputed_objective': objective, 'full_primal_fluxes': {'EX_a_e': -objective, 'BIOMASS_TEST': objective},
            'primal_certificate': {'max_abs_S_residual': 1e-10, 'max_bound_violation': 0.,
                                  'n_S_rows_over_1e-6': 0, 'n_bounds_over_1e-6': 0},
            'max_fsum_residual': 1e-10, 'passes_original_primal_gate': True}


def records():
    return [record(method) for method in C.METHOD_IDS]


def test_all_three_consistent_methods_accepted():
    assessment = C.assess_case(records(), RECIPE, REACTIONS)
    assert assessment['accepted_cross_solver_evidence']
    assert assessment['all_three_methods_accepted']
    assert len(assessment['pairwise_reported_optimal_objectives']) == 3


def test_rejected_primal_is_retained_and_agreeing_second_highs_can_support_evidence():
    rows = records()
    rows[2]['primal_certificate']['max_abs_S_residual'] = 1e-5
    before = deepcopy(rows)
    assessment = C.assess_case(rows, RECIPE, REACTIONS)
    assert rows == before
    assert not assessment['methods'][C.METHOD_IDS[2]]['accepted_for_comparison']
    assert assessment['accepted_cross_solver_evidence']
    assert not assessment['all_three_methods_accepted']


def test_rejected_vector_objective_still_blocks_cherry_picked_agreement():
    rows = records()
    rows[2] = record(C.METHOD_IDS[2], 2.)
    rows[2]['primal_certificate']['max_abs_S_residual'] = 1e-5
    assessment = C.assess_case(rows, RECIPE, REACTIONS)
    assert assessment['glpk_and_at_least_one_highs_accepted']
    assert not assessment['all_reported_optimal_objectives_agree']
    assert not assessment['accepted_cross_solver_evidence']
    assert len(assessment['pairwise_reported_optimal_objectives']) == 3


def test_individually_feasible_records_are_not_aggregate_success_if_objectives_disagree():
    rows = records()
    rows[2] = record(C.METHOD_IDS[2], 2.)
    assessment = C.assess_case(rows, RECIPE, REACTIONS)
    assert all(value['accepted_for_comparison'] for value in assessment['methods'].values())
    assert not assessment['all_three_methods_accepted']
    assert not assessment['accepted_cross_solver_evidence']


def test_compensated_gate_is_additive_and_does_not_rewrite_original_gate():
    row = record('glpk_native')
    row['max_fsum_residual'] = 2e-8
    assessment = C.assess_record(row, RECIPE, REACTIONS)
    assert assessment['original_primal_gate_recomputed']
    assert row['passes_original_primal_gate'] is True
    assert not assessment['compensated_residual_gate']
    assert not assessment['accepted_for_comparison']


@pytest.mark.parametrize('corruption', ['nan_sparse', 'negative_bound', 'nan_fsum', 'missing_vector',
    'extra_flux', 'nan_flux', 'bool_flux', 'nonfinite_objective', 'objective_vector_disagreement', 'nonoptimal_status'])
def test_invalid_method_record_is_not_accepted(corruption):
    row = record('glpk_native')
    if corruption == 'nan_sparse':
        row['primal_certificate']['max_abs_S_residual'] = float('nan')
    elif corruption == 'negative_bound':
        row['primal_certificate']['max_bound_violation'] = -.1
    elif corruption == 'nan_fsum':
        row['max_fsum_residual'] = float('nan')
    elif corruption == 'missing_vector':
        del row['full_primal_fluxes']['EX_a_e']
    elif corruption == 'extra_flux':
        row['full_primal_fluxes']['unapproved'] = 0.
    elif corruption == 'nan_flux':
        row['full_primal_fluxes']['EX_a_e'] = float('nan')
    elif corruption == 'bool_flux':
        row['full_primal_fluxes']['BIOMASS_TEST'] = True
    elif corruption == 'nonfinite_objective':
        row['objective'] = float('nan')
    elif corruption == 'objective_vector_disagreement':
        row['objective'] = 2.
    else:
        row['status'] = 'Time limit reached'
    assert not C.assess_record(row, RECIPE, REACTIONS)['accepted_for_comparison']


def test_infeasibility_conflict_and_missing_optimal_objective_are_explicit():
    rows = records()
    rows[2] = {'method': {'id': C.METHOD_IDS[2]}, 'status': 'Infeasible', 'objective': None}
    assessment = C.assess_case(rows, RECIPE, REACTIONS)
    assert assessment['reported_feasibility_conflict']
    assert not assessment['accepted_cross_solver_evidence']
    rows[2] = record(C.METHOD_IDS[2])
    rows[2]['objective'] = None
    assessment = C.assess_case(rows, RECIPE, REACTIONS)
    assert not assessment['all_optimal_objectives_finite']
    assert not assessment['accepted_cross_solver_evidence']


def test_timeout_is_a_reported_method_failure_not_fabricated_zero_or_panel_abort():
    rows = records()
    rows[2] = {'method': {'id': C.METHOD_IDS[2]}, 'status': 'Time limit reached', 'objective': None}
    assessment = C.assess_case(rows, RECIPE, REACTIONS)
    assert assessment['accepted_cross_solver_evidence']
    assert not assessment['all_three_methods_accepted']
    assert rows[2]['objective'] is None


@pytest.mark.parametrize('certificate', [None, [], 'unavailable'])
def test_malformed_certificate_is_a_rejection_without_aborting(certificate):
    row = record('glpk_native')
    row['primal_certificate'] = certificate
    assert not C.assess_record(row, RECIPE, REACTIONS)['accepted_for_comparison']


def test_exception_record_cannot_be_accepted_even_if_it_contains_a_good_vector():
    row = record('glpk_native')
    row['exception'] = {'type': 'RuntimeError', 'message': 'later diagnostic failure'}
    assessment = C.assess_record(row, RECIPE, REACTIONS)
    assert assessment['original_primal_gate_recomputed']
    assert not assessment['accepted_for_comparison']


@pytest.mark.parametrize('status', ['Unbounded', 'Infeasible or unbounded', 'Unbounded or infeasible', 'infeasible_or_unbounded'])
def test_infinite_or_infeasible_claim_conflicts_with_reported_finite_optimum(status):
    rows = records()
    rows[2] = {'method': {'id': C.METHOD_IDS[2]}, 'status': status, 'objective': None}
    assessment = C.assess_case(rows, RECIPE, REACTIONS)
    assert assessment['glpk_and_at_least_one_highs_accepted']
    assert assessment['reported_status_conflict']
    assert not assessment['accepted_cross_solver_evidence']


def test_all_infeasible_records_do_not_claim_primal_cross_solver_evidence():
    rows = [{'method': {'id': method}, 'status': 'Infeasible', 'objective': None} for method in C.METHOD_IDS]
    assessment = C.assess_case(rows, RECIPE, REACTIONS)
    assert assessment['n_reported_optimal'] == 0
    assert not assessment['accepted_cross_solver_evidence']
    assert not assessment['all_reported_optimal_objectives_agree']


@pytest.mark.parametrize('kind', ['missing', 'duplicate', 'reordered'])
def test_complete_fixed_method_inventory_required(kind):
    rows = records()
    if kind == 'missing':
        rows.pop()
    elif kind == 'duplicate':
        rows[2] = deepcopy(rows[1])
    else:
        rows.reverse()
    with pytest.raises(ValueError, match='inventory'):
        C.assess_case(rows, RECIPE, REACTIONS)


@pytest.fixture
def toy():
    model = cobra.Model('capture_toy')
    met = cobra.Metabolite('a_e', compartment='e', formula='C', charge=0)
    for rid, coefficient, bounds in [('EX_a_e', -1., (-1., 1000.)), ('BIOMASS_TEST', -1., (0., 1000.))]:
        reaction = cobra.Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
        reaction.add_metabolites({met: coefficient})
        model.add_reactions([reaction])
    model.objective = 'BIOMASS_TEST'
    return model


def read_gzip(path):
    with gzip.open(path, 'rt') as stream:
        return json.load(stream)


def test_invalid_returned_vector_is_persisted_before_acceptance(toy, tmp_path, monkeypatch):
    bad = record('glpk_native')
    bad['primal_certificate']['max_abs_S_residual'] = .1
    monkeypatch.setattr(C.D, 'run_method', lambda *args: bad)
    output = tmp_path / 'raw.json.gz'
    method = {'id': 'glpk_native', 'solver_tolerance': 1e-9}
    captured = C.capture_and_store(toy, method, RECIPE, C.R.signature(toy), output)
    assert read_gzip(output)['full_primal_fluxes'] == bad['full_primal_fluxes']
    assert not C.assess_record(captured, RECIPE, REACTIONS)['accepted_for_comparison']


def test_vector_is_preserved_even_if_post_solve_algebra_guard_fails(toy, tmp_path, monkeypatch):
    signature = C.R.signature(toy)
    def mutate(*args):
        toy.reactions.EX_a_e.lower_bound = -2.
        return record('glpk_native')
    monkeypatch.setattr(C.D, 'run_method', mutate)
    output = tmp_path / 'raw.json.gz'
    with pytest.raises(ValueError, match='raw result retained'):
        C.capture_and_store(toy, {'id': 'glpk_native', 'solver_tolerance': 1e-9}, RECIPE, signature, output)
    assert read_gzip(output)['full_primal_fluxes'] == record('glpk_native')['full_primal_fluxes']


def test_capture_exception_and_nonfinite_values_are_preserved_as_failures(toy, tmp_path, monkeypatch):
    def fail(*args):
        raise RuntimeError('synthetic solver error')
    monkeypatch.setattr(C.D, 'run_method', fail)
    output = tmp_path / 'error.json.gz'
    result = C.capture_and_store(toy, {'id': 'glpk_native', 'solver_tolerance': 1e-9}, RECIPE, C.R.signature(toy), output)
    assert read_gzip(output)['exception']['message'] == 'synthetic solver error'
    assert not C.assess_record(result, RECIPE, REACTIONS)['accepted_for_comparison']
    bad = record('glpk_native')
    bad['full_primal_fluxes']['EX_a_e'] = float('nan')
    monkeypatch.setattr(C.D, 'run_method', lambda *args: bad)
    output = tmp_path / 'nonfinite.json.gz'
    C.capture_and_store(toy, {'id': 'glpk_native', 'solver_tolerance': 1e-9}, RECIPE, C.R.signature(toy), output)
    assert read_gzip(output)['full_primal_fluxes']['EX_a_e'] == 'nan'


def test_recipe_fixes_full_panel_and_distinct_solver_tolerances():
    recipe = C.R.read(C.PLAN)
    assert (recipe['condition_count'], recipe['solve_count']) == (43, 129)
    assert [row['id'] for row in recipe['methods']] == list(C.METHOD_IDS)
    assert [row['solver_tolerance'] for row in recipe['methods']] == [1e-9, 1e-9, 1e-10]
    effective = C.supported_options(recipe['methods'])
    assert effective['highs_primal_simplex']['simplex_strategy'] == 4
    assert effective['highs_simplex_tighter']['simplex_strategy'] == 1
    assert effective['highs_simplex_tighter']['primal_feasibility_tolerance'] == 1e-10
    assert len(recipe['diagnostic_inputs']) == len(set(recipe['diagnostic_inputs']))
    assert all((C.ROOT / path).is_file() for path in recipe['diagnostic_inputs'])
