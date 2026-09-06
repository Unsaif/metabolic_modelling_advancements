import cobra
import numpy as np
import pytest
from scipy.optimize import LinearConstraint
from scipy.sparse import csr_matrix

from gembench.gapfill import _solve_highs, apply_gapfill, gapfill
from gembench.media import Medium, apply_medium


def growth_models(capacity=1000, yield_coefficient=1):
    model = cobra.Model("draft")
    ae = cobra.Metabolite("a_e", compartment="e")
    ac = cobra.Metabolite("a_c", compartment="c")
    bc = cobra.Metabolite("b_c", compartment="c")
    exchange = cobra.Reaction("EX_a_e", lower_bound=-1000)
    exchange.add_metabolites({ae: -1})
    transport = cobra.Reaction("TRANSPORT")
    transport.add_metabolites({ae: -1, ac: 1})
    growth = cobra.Reaction("Growth")
    growth.add_metabolites({bc: -1})
    model.add_reactions([exchange, transport, growth])
    model.objective = growth
    universe = cobra.Model("universe")
    fill = cobra.Reaction("FILL", upper_bound=capacity)
    fill.add_metabolites({ac: -1, bc: yield_coefficient})
    universe.add_reactions([fill])
    medium = Medium("test", "test", {"EX_a_e": -1000})
    return model, universe, medium


def test_candidate_capacity_cannot_be_inflated_to_claim_success():
    model, universe, medium = growth_models(capacity=0.01)
    result = gapfill(model, universe, medium, min_growth=0.05, verbose=False)
    assert result.status.startswith("failed:")
    assert result.added_reactions == []


def test_candidate_flux_above_arbitrary_big_m_remains_available():
    model, universe, medium = growth_models(yield_coefficient=0.0001)
    result = gapfill(model, universe, medium, min_growth=0.05, verbose=False)
    assert result.status == "Optimal"
    assert result.added_reactions == ["FILL"]
    assert result.growth_after == pytest.approx(0.1)


def test_gapfill_uses_requested_biomass_not_original_objective():
    model, universe, medium = growth_models()
    waste = cobra.Reaction("WASTE")
    waste.add_metabolites({model.metabolites.a_c: -1})
    model.add_reactions([waste])
    model.objective = waste
    result = gapfill(model, universe, medium, verbose=False)
    assert result.added_reactions == ["FILL"]
    assert result.growth_before == 0
    assert result.growth_after == pytest.approx(1000)
    assert "WASTE" in str(model.objective.expression)


def test_unconventionally_named_boundary_source_is_not_a_candidate():
    model, universe, medium = growth_models(capacity=0.01)
    free = cobra.Reaction("MAGIC_SOURCE")
    free.add_metabolites({model.metabolites.b_c: 1})
    universe.add_reactions([free])
    result = gapfill(model, universe, medium, verbose=False)
    assert result.status.startswith("failed:")
    assert result.n_candidates == 1


def test_oxygen_transport_is_allowed_and_oxygen_release_stays_blocked():
    model, universe, medium = growth_models()
    ae = model.metabolites.a_e
    ac = model.metabolites.a_c
    ae.id = "o2_e"
    ac.id = "o2_c"
    model.reactions.EX_a_e.id = "EX_o2_e"
    # The only missing reaction is a reversible O2 transporter. Movement between
    # compartments has zero net oxygen production and must not be disabled.
    tr = model.reactions.TRANSPORT.copy()
    tr.bounds = (-1000, 1000)
    model.remove_reactions([model.reactions.TRANSPORT])
    fill = universe.reactions.FILL.copy()
    fill.add_metabolites({next(m for m in fill.metabolites if m.id == "a_c"): 1})
    fill.add_metabolites({ac: -1})
    model.add_reactions([fill])
    universe = cobra.Model("oxygen_universe")
    universe.add_reactions([tr])
    medium = Medium("oxygen", "", {"EX_o2_e": -10})
    result = gapfill(model, universe, medium, verbose=False)
    assert result.added_reactions == ["TRANSPORT"]
    filled = apply_gapfill(model, universe, result.added_reactions)
    assert filled.reactions.TRANSPORT.bounds == (-1000, 1000)

    consuming = cobra.Reaction("OXYGENASE", lower_bound=-500)
    consuming.add_metabolites({ac: -1, model.metabolites.b_c: 1})
    universe.add_reactions([consuming])
    filled = apply_gapfill(model, universe, ["OXYGENASE"])
    assert filled.reactions.OXYGENASE.bounds == (0, 1000)
    unrestricted = apply_gapfill(model, universe, ["OXYGENASE"], forbid_o2_production=False)
    assert unrestricted.reactions.OXYGENASE.bounds == (-500, 1000)


def test_candidate_mandatory_flux_is_conditional_on_selection():
    model, universe, medium = growth_models()
    universe.reactions.FILL.lower_bound = 10
    medium.uptakes["EX_a_e"] = -1
    result = gapfill(model, universe, medium, verbose=False)
    assert result.status.startswith("failed:")


def test_highs_time_limit_without_incumbent_is_not_a_solution():
    # A knapsack with no preprocessing simplification and a zero solve budget.
    rng = np.random.default_rng(12)
    cost = -rng.integers(1, 100, size=100).astype(float)
    constraint = LinearConstraint(csr_matrix(rng.integers(1, 100, size=(30, 100))),
                                  np.full(30, -np.inf), np.full(30, 1200))
    solution, status = _solve_highs(cost, [constraint], np.zeros(100), np.ones(100), 0, 0, 0)
    assert solution is None
    assert "time limit" in status.lower()


def test_successful_gapfill_reproduces_growth_after_application():
    model, universe, medium = growth_models(capacity=2)
    result = gapfill(model, universe, medium, verbose=False)
    filled = apply_gapfill(model, universe, result.added_reactions)
    apply_medium(filled, medium)
    assert filled.slim_optimize() == pytest.approx(result.growth_after)
    assert filled.reactions.FILL.upper_bound == 2
    assert "FILL" not in model.reactions
