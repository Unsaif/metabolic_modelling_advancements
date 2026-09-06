from types import SimpleNamespace

import cobra
import pandas as pd
import pytest

from gembench.cofactor import pool_balance_certificate, probe_metabolite_production


def cyclic_model():
    model = cobra.Model('catalytic_pool')
    mets = {mid: cobra.Metabolite(mid, compartment='c') for mid in ('food', 'product', 'q', 'qh')}
    for rid, stoich in [('FOOD', {'food': 1}), ('REDUCE', {'food': -1, 'q': -1, 'qh': 1, 'product': 1}),
                        ('OXIDIZE', {'qh': -1, 'q': 1}), ('Growth', {'product': -1})]:
        reaction = cobra.Reaction(rid, lower_bound=0, upper_bound=10)
        reaction.add_metabolites({mets[mid]: value for mid, value in stoich.items()})
        model.add_reactions([reaction])
    model.objective = 'Growth'
    return model


def test_growth_does_not_prove_cofactor_production_and_demand_blocks_closed_pool():
    model = cyclic_model()
    assert model.slim_optimize() == pytest.approx(10)
    report = pool_balance_certificate(model, {'q': 1, 'qh': 1})
    assert report['exactly_conserved_in_stoichiometry']
    assert report['direction_compatible_supply_reactions'] == []
    assert all(r['maximum_demand_flux'] == 0 for r in probe_metabolite_production(model, ['q', 'qh']))
    model.reactions.Growth.add_metabolites({model.metabolites.qh: -0.0001})
    assert model.slim_optimize() == pytest.approx(0)
    assert not pool_balance_certificate(model, {'q': 1, 'qh': 1})['exactly_conserved_in_stoichiometry']


def test_real_precursor_route_opens_pool_but_closed_medium_still_blocks_production():
    model = cyclic_model()
    synthesis = cobra.Reaction('SYNTHESIS', lower_bound=0, upper_bound=10)
    synthesis.add_metabolites({model.metabolites.food: -1, model.metabolites.q: 1})
    model.add_reactions([synthesis])
    report = pool_balance_certificate(model, {'q': 1, 'qh': 1})
    assert report['direction_compatible_supply_reactions'] == ['SYNTHESIS']
    assert probe_metabolite_production(model, ['q'])[0]['maximum_demand_flux'] == pytest.approx(10)
    model.reactions.FOOD.upper_bound = 0
    assert not probe_metabolite_production(model, ['q'])[0]['producible_at_threshold']
    assert not pool_balance_certificate(model, {'q': 1, 'qh': 1})['exactly_conserved_in_stoichiometry']


def test_small_stoichiometric_leak_is_not_rounded_to_conservation():
    model = cyclic_model()
    model.reactions.OXIDIZE.add_metabolites({model.metabolites.q: 1e-12})
    report = pool_balance_certificate(model, {'q': 1, 'qh': 1})
    assert not report['exactly_conserved_in_stoichiometry']
    assert report['direction_compatible_supply_reactions'] == ['OXIDIZE']


def test_reverse_source_and_disabled_source_are_distinguished():
    model = cyclic_model()
    source = cobra.Reaction('REVERSE_SOURCE', lower_bound=-3, upper_bound=0)
    source.add_metabolites({model.metabolites.q: -1})
    model.add_reactions([source])
    assert pool_balance_certificate(model, {'q': 1, 'qh': 1})['direction_compatible_supply_reactions'] == ['REVERSE_SOURCE']
    assert probe_metabolite_production(model, ['q'])[0]['maximum_demand_flux'] == pytest.approx(3)
    source.lower_bound = 0
    assert pool_balance_certificate(model, {'q': 1, 'qh': 1})['direction_compatible_supply_reactions'] == []


def test_probe_restores_objective_bounds_and_reactions_and_respects_maintenance():
    model = cyclic_model()
    model.objective_direction = 'min'
    model.reactions.Growth.lower_bound = 1
    before = (str(model.objective.expression), model.objective_direction, {r.id: r.bounds for r in model.reactions})
    assert probe_metabolite_production(model, ['product'])[0]['maximum_demand_flux'] == pytest.approx(9)
    assert (str(model.objective.expression), model.objective_direction, {r.id: r.bounds for r in model.reactions}) == before


@pytest.mark.parametrize('status,objective', [('time_limit', 0), ('infeasible', None), ('optimal', float('nan'))])
def test_failure_is_unknown_and_temporary_probe_is_removed(monkeypatch, status, objective):
    model = cyclic_model()
    before = {r.id for r in model.reactions}
    monkeypatch.setattr(model, 'optimize', lambda: SimpleNamespace(status=status, objective_value=objective))
    with pytest.raises(RuntimeError, match='Production probe.*failed'):
        probe_metabolite_production(model, ['q'])
    assert {r.id for r in model.reactions} == before


@pytest.mark.parametrize('weights', [{}, {'missing': 1}, {'q': 0}, {'q': -1}, {'q': float('inf')}, {'q': float('nan')}])
def test_invalid_or_unrepresented_pool_is_not_a_clean_result(weights):
    with pytest.raises(ValueError):
        pool_balance_certificate(cyclic_model(), weights)


def test_nonzero_balance_constraint_does_not_support_closed_pool_certificate():
    model = cyclic_model()
    model.metabolites.q.constraint.ub = 1
    with pytest.raises(ValueError, match='zero steady-state balance'):
        pool_balance_certificate(model, {'q': 1, 'qh': 1})


def test_custom_solver_supply_does_not_support_stoichiometric_pool_certificate():
    model = cyclic_model()
    supply = model.problem.Variable('extra_supply', lb=0, ub=1)
    model.add_cons_vars(supply)
    model.metabolites.q.constraint.set_linear_coefficients({supply: 1})
    # The current solver can supply the pool, although the reaction metadata
    # still contains only one-for-one quinone interconversion.
    assert probe_metabolite_production(model, ['q'])[0]['maximum_demand_flux'] == pytest.approx(1)
    with pytest.raises(ValueError, match='solver balance differs from stored stoichiometry'):
        pool_balance_certificate(model, {'q': 1, 'qh': 1})


def test_altered_reaction_coefficient_in_solver_row_is_rejected():
    model = cyclic_model()
    model.metabolites.q.constraint.set_linear_coefficients({model.reactions.REDUCE.forward_variable: -2})
    with pytest.raises(ValueError, match='solver balance differs from stored stoichiometry'):
        pool_balance_certificate(model, {'q': 1, 'qh': 1})


def test_probe_respects_additional_constraint_and_preserves_custom_objective():
    model = cyclic_model()
    limit = model.problem.Constraint(model.reactions.REDUCE.flux_expression, ub=2, name='production_limit')
    model.add_cons_vars(limit)
    model.objective = model.problem.Objective(model.reactions.Growth.flux_expression + 2 * model.reactions.OXIDIZE.flux_expression,
                                              direction='min')
    before = (str(model.objective.expression), model.objective_direction,
              str(limit.expression), limit.lb, limit.ub, len(model.constraints))
    assert probe_metabolite_production(model, ['product'])[0]['maximum_demand_flux'] == pytest.approx(2)
    assert (str(model.objective.expression), model.objective_direction,
            str(limit.expression), limit.lb, limit.ub, len(model.constraints)) == before
    assert pool_balance_certificate(model, {'q': 1, 'qh': 1})['exactly_conserved_in_stoichiometry']


def test_nonfinite_flux_is_unknown_and_probe_state_is_restored(monkeypatch):
    model = cyclic_model()
    before = (str(model.objective.expression), model.objective_direction, {r.id: r.bounds for r in model.reactions})

    def nonfinite_solution():
        return SimpleNamespace(status='optimal', objective_value=0.0,
                               fluxes=pd.Series(float('nan'), index=[r.id for r in model.reactions]))

    monkeypatch.setattr(model, 'optimize', nonfinite_solution)
    with pytest.raises(RuntimeError, match='nonfinite fluxes'):
        probe_metabolite_production(model, ['q'])
    assert (str(model.objective.expression), model.objective_direction, {r.id: r.bounds for r in model.reactions}) == before


def test_zero_production_does_not_reach_capacity_when_threshold_equals_capacity():
    model = cyclic_model()
    result = probe_metabolite_production(model, ['q'], capacity=1e-6, threshold=1e-6)[0]
    assert not result['capacity_reached']
    assert not result['producible_at_threshold']
    saturated = probe_metabolite_production(model, ['product'], capacity=1e-6, threshold=1e-6)[0]
    assert saturated['capacity_reached']
