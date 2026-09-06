"""Differential-study guards and dual solvers on synthetic LPs only."""
from types import SimpleNamespace

import cobra
import numpy as np
import pytest

from gembench.media import Medium
from scripts import run_medium_preparation as R


RECIPE = {'solver_tolerance': 1e-9, 'residual_tolerance': 1e-8, 'objective_tolerance': 1e-8}


@pytest.fixture
def toy():
    model = cobra.Model('generic_objective_toy')
    external = cobra.Metabolite('a_e', compartment='e', formula='C', charge=0)
    internal = cobra.Metabolite('a_c', compartment='c', formula='C', charge=0)
    rows = [('EX_a_e', {external: -1.}, (-5., 1000.), ''),
            ('TRANSPORT', {external: -1., internal: 1.}, (0., 1000.), 'g_transport'),
            ('BIOMASS_TEST', {internal: -1.}, (0., 1000.), '')]
    for rid, coefficients, bounds, rule in rows:
        reaction = cobra.Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
        reaction.add_metabolites(coefficients)
        reaction.gene_reaction_rule = rule
        model.add_reactions([reaction])
    model.objective = 'BIOMASS_TEST'
    return model


@pytest.mark.parametrize('change,section', [('bounds', 'reactions'), ('stoichiometry', 'reactions'),
                                          ('gpr', 'reactions'), ('formula', 'metabolites'),
                                          ('charge', 'metabolites'), ('objective', 'objective')])
def test_signature_detects_changes_material_to_declared_equivalence(toy, change, section):
    baseline = R.signature(toy)
    if change == 'bounds':
        toy.reactions.EX_a_e.lower_bound = -4.
    elif change == 'stoichiometry':
        toy.reactions.TRANSPORT.add_metabolites({toy.metabolites.a_c: 1.})
    elif change == 'gpr':
        toy.reactions.TRANSPORT.gene_reaction_rule = 'other_gene'
    elif change == 'formula':
        toy.metabolites.a_c.formula = 'CH2'
    elif change == 'charge':
        toy.metabolites.a_c.charge = -1
    else:
        toy.objective = {toy.reactions.BIOMASS_TEST: 2.}
    changed = R.signature(toy)
    assert section in R.differences(baseline, changed)
    assert R.digest(changed) != R.digest(baseline)


def test_signature_ignores_display_names_and_compartment_labels_only(toy):
    baseline = R.signature(toy)
    toy.name = 'another display name'
    toy.reactions.TRANSPORT.name = 'another reaction name'
    toy.metabolites.a_e.name = 'another metabolite name'
    toy.metabolites.a_e.compartment = 'C_e'
    toy.metabolites.a_c.compartment = 'C_c'
    assert R.signature(toy) == baseline
    assert R.differences(baseline, R.signature(toy)) == {}


@pytest.mark.parametrize('kind', ['coupling', 'balance', 'minimize'])
def test_guard_rejects_features_lost_by_matrix_translation(toy, kind):
    if kind == 'coupling':
        toy.add_cons_vars([toy.problem.Constraint(toy.reactions.BIOMASS_TEST.flux_expression,
                                                ub=100., name='extra_constraint')])
    elif kind == 'balance':
        toy.metabolites.a_c.constraint.ub = 1.
    else:
        toy.objective_direction = 'min'
    with pytest.raises(ValueError, match='Unsupported|supports maximization'):
        R.signature(toy)
    with pytest.raises(ValueError, match='Unsupported|supports maximization'):
        R.checked_solves(toy, RECIPE)


def test_both_solvers_support_generic_objective_and_retain_full_primal_vectors(toy):
    assert 'Growth' not in toy.reactions
    records = R.checked_solves(toy, RECIPE)
    assert [record['solver'] for record in records] == ['glpk', 'highs']
    for record in records:
        assert record['status'] == 'optimal'
        assert record['objective'] == pytest.approx(5.)
        assert set(record['full_primal_fluxes']) == {reaction.id for reaction in toy.reactions}
        assert record['full_primal_fluxes']['BIOMASS_TEST'] == pytest.approx(5.)
        assert record['primal_certificate']['max_abs_S_residual'] <= RECIPE['residual_tolerance']
        assert record['primal_certificate']['max_bound_violation'] <= RECIPE['residual_tolerance']


def test_infeasibility_is_retained_as_status_with_null_objective(toy):
    toy.reactions.BIOMASS_TEST.lower_bound = 6.
    with pytest.warns(UserWarning, match="Solver status is 'infeasible'"):
        records = R.checked_solves(toy, RECIPE)
    assert [record['status'] for record in records] == ['infeasible', 'infeasible']
    for record in records:
        assert record['objective'] is None
        assert 'full_primal_fluxes' not in record
        assert 'primal_certificate' not in record


@pytest.mark.parametrize('status', ['Time limit reached', 'Numerical error', 'Unbounded'])
def test_solver_failures_are_not_converted_into_no_growth(toy, monkeypatch, status):
    monkeypatch.setattr(R.W, 'solve_highs', lambda *a, **k: SimpleNamespace(status=status, x=None, objective=None))
    with pytest.raises(RuntimeError, match='Unresolved solver status'):
        R.checked_solves(toy, RECIPE)


def test_solver_feasibility_disagreement_is_an_error(toy, monkeypatch):
    monkeypatch.setattr(R.W, 'solve_highs', lambda *a, **k: SimpleNamespace(status='Infeasible', x=None, objective=None))
    with pytest.raises(RuntimeError, match='classifications disagree'):
        R.checked_solves(toy, RECIPE)


@pytest.mark.parametrize('invalid', [float('nan'), float('inf'), -.1, .1])
def test_invalid_residual_certificate_is_rejected(toy, monkeypatch, invalid):
    monkeypatch.setattr(R.W, 'certify', lambda *a, **k: {'max_abs_S_residual': invalid,
        'max_bound_violation': 0., 'n_S_rows_over_1e-6': 0, 'n_bounds_over_1e-6': 0})
    with pytest.raises(RuntimeError, match='certificate'):
        R.checked_solves(toy, RECIPE)


def test_nonfinite_highs_vector_is_rejected(toy, monkeypatch):
    monkeypatch.setattr(R.W, 'solve_highs', lambda *a, **k:
        SimpleNamespace(status='Optimal', x=np.array([-5., float('nan'), 5.]), objective=5.))
    with pytest.raises(RuntimeError, match='Nonfinite'):
        R.checked_solves(toy, RECIPE)


def test_condition_carbon_overrides_trace_capacity_without_mutating_media(monkeypatch):
    medium = Medium('synthetic', '', {'EX_a_e': -.001, 'EX_pi_e': -1000.})
    monkeypatch.setattr(R, 'base_medium', lambda name: medium)
    condition = SimpleNamespace(media='synthetic', exchanges=['EX_a_e', 'EX_b_e'])
    assembled = R.medium_for(condition)
    assert assembled.uptakes == {'EX_a_e': -10., 'EX_b_e': -10., 'EX_pi_e': -1000.}
    assert medium.uptakes == {'EX_a_e': -.001, 'EX_pi_e': -1000.}
