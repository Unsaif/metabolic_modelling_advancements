"""Only synthetic networks; no real study outcomes are used by these tests."""
from types import SimpleNamespace

import cobra
import pytest

from scripts import verify_quinone_repair_solvers as V


@pytest.fixture
def model():
    m = cobra.Model('synthetic_quinone')
    a = cobra.Metabolite('a_c', compartment='c')
    q = cobra.Metabolite('q8_c', compartment='c')
    qh = cobra.Metabolite('q8h2_c', compartment='c')
    supply = cobra.Reaction('SUPPLY', lower_bound=0, upper_bound=2)
    supply.add_metabolites({a: 1})
    growth = cobra.Reaction('Growth', lower_bound=0, upper_bound=1000)
    growth.add_metabolites({a: -1})
    growth.gene_reaction_rule = 'gene_a'
    redox = cobra.Reaction('REDOX', lower_bound=-1000, upper_bound=1000)
    redox.add_metabolites({q: -1, qh: 1})
    m.add_reactions([supply, growth, redox])
    m.objective = growth
    return m


def test_two_solvers_growth_and_closed_pool(model):
    result = V.check_case(model, 2.0, 0.001)
    assert result['valid'], result
    with model:
        demand = cobra.Reaction('DIAG_DEMAND_q8_c', lower_bound=0, upper_bound=1000)
        demand.add_metabolites({model.metabolites.q8_c: -1})
        model.add_reactions([demand])
        model.objective = demand
        result = V.check_case(model, 0.0, 1e-6)
        assert result['valid'], result
    with model:
        model.genes.gene_a.knock_out()
        result = V.check_case(model, 0.0, 0.001)
        assert result['valid'], result
    assert V.check_case(model, 2.0, 0.001)['valid']


def test_rejects_wrong_saved_optimum(model):
    result = V.check_case(model, 0.0, 0.001)
    assert not result['valid']
    assert any('saved result' in p for p in result['problems'])


@pytest.mark.parametrize('customization', ['row', 'variable', 'rhs', 'stoichiometry', 'objective', 'minimize', 'split_bound'])
def test_rejects_features_dropped_by_matrix_translation(model, customization):
    if customization == 'row':
        model.add_cons_vars(model.problem.Constraint(model.reactions.Growth.flux_expression, ub=1, name='custom'))
    elif customization == 'variable':
        model.add_cons_vars(model.problem.Variable('custom', lb=0))
    elif customization == 'rhs':
        model.metabolites.a_c.constraint.ub = 1
    elif customization == 'stoichiometry':
        model.metabolites.a_c.constraint.set_linear_coefficients({model.reactions.Growth.forward_variable: -2})
    elif customization == 'objective':
        model.objective = model.problem.Objective(model.reactions.Growth.forward_variable)
    elif customization == 'minimize':
        model.objective_direction = 'min'
    else:
        model.reactions.Growth.forward_variable.ub = 1
    with pytest.raises(ValueError, match='Unsupported|supports maximization'):
        V.assert_ordinary_lp(model)


def test_rejects_bad_primal_certificate_even_when_objectives_agree(model, monkeypatch):
    def bad(*_):
        return {'status': 'optimal', 'objective': 2.0,
                'primal_certificate': {'max_abs_S_residual': 1e-4, 'max_bound_violation': 0.0}}
    monkeypatch.setattr(V, 'solve_case', bad)
    result = V.check_case(model, 2.0, 0.001)
    assert not result['valid']
    assert all('residual' in p for p in result['problems'])


def test_target_skips_require_actual_absence_or_orphan(model):
    mapping = {'model_to_browser': {'gene_a': 'locus_a'}}
    recipe = {'targeted_gene_loci': ['locus_a', 'locus_missing']}
    saved = {'targeted_gene_predictions': [
        {'locus': 'locus_a', 'model_gene_id': 'gene_a', 'associated_reactions': ['Growth'], 'status': 'optimal', 'growth': 0},
        {'locus': 'locus_missing', 'model_gene_id': None, 'associated_reactions': [], 'status': 'unrepresented_gene', 'growth': None}]}
    assert [r[2] for r in V.resolve_targets(model, mapping, saved, recipe)] == [None, 'unrepresented_gene']
    saved['targeted_gene_predictions'][0]['status'] = 'unrepresented_function'
    with pytest.raises(ValueError, match='Invalid saved targeted gene'):
        V.resolve_targets(model, mapping, saved, recipe)


def test_physical_context_requires_matching_bounds(model):
    recipe = {'physical_settings': {}, 'protocol': {}}
    condition = SimpleNamespace(key='toy', bigg_ids=[], exchanges=[])
    saved = {'condition': 'toy', 'settings': {},
             'objective': {'expression': str(model.objective.expression), 'direction': 'max'},
             'growth_bounds': [0, 1000], 'atp_maintenance_bounds': None, 'mandatory_flux_bounds': {},
             'final_exchange_bounds': {}, 'declared_carbon_source_ids': [], 'declared_carbon_exchanges': [],
             'declared_carbon_uptake': -10, 'medium_completion_added': [], 'missing_medium_components': [],
             'effective_feasibility_tolerance': 1e-9}
    # The toy's one-metabolite boundaries may be auto-classified as exchanges.
    saved['final_exchange_bounds'] = {r.id: list(r.bounds) for r in model.exchanges}
    V.assert_physical_context(model, saved, recipe, condition, [], [])
    saved['growth_bounds'] = [1, 1000]
    with pytest.raises(ValueError, match='growth_bounds'):
        V.assert_physical_context(model, saved, recipe, condition, [], [])


def test_changed_inputs_and_duplicate_json_are_rejected(tmp_path):
    path = tmp_path / 'input.json'
    path.write_text('{"a":1}')
    records = [V.file_record(path)]
    V.unchanged(records)
    path.write_text('{"a":2}')
    with pytest.raises(ValueError, match='changed'):
        V.unchanged(records)
    path.write_text('{"a":1,"a":2}')
    with pytest.raises(ValueError, match='Duplicate JSON'):
        V.read_json(path)
