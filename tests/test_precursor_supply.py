import math

import cobra
import pytest

from scripts import run_precursor_supply as R


@pytest.fixture
def toy():
    model = cobra.Model('precursor_toy')
    ids = ['carbon_c', '4hbz_c', '4hbz_e', 'q8_c', 'q8h2_c', 'T4hcinnm_e']
    mets = {mid: cobra.Metabolite(mid, compartment='c') for mid in ids}
    model.add_metabolites(list(mets.values()))
    rows = [
        ('EX_carbon', {'carbon_c': -1}, (-10, 1000), ''),
        ('CHRPL', {'carbon_c': -1, '4hbz_c': 1}, (0, 1000), 'g_ubic'),
        ('HBZOPT', {'4hbz_c': -1, 'q8h2_c': 1}, (0, 1000), 'g_ubia'),
        ('UHBZ1t_pp', {'4hbz_e': -1, '4hbz_c': 1}, (-1000, 1000), 'g_transport'),
        ('EX_4hbz_e', {'4hbz_e': -1}, (0, 1000), ''),
        ('EX_T4hcinnm_e', {'T4hcinnm_e': -1}, (0, 1000), ''),
        ('4HBHYOX', {'4hbz_c': -1, 'carbon_c': 1}, (0, 1000), ''),
        ('sink_2ohph_c', {'4hbz_c': -1}, (0, 1000), ''),
        ('Growth', {'carbon_c': -1, 'q8h2_c': -.1}, (0, 1000), '')]
    for rid, stoich, bounds, rule in rows:
        r = cobra.Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
        r.add_metabolites({mets[mid]: c for mid, c in stoich.items()})
        r.gene_reaction_rule = rule
        model.add_reactions([r])
    model.objective = 'Growth'
    recipe = {'solver_tolerance': 1e-9, 'residual_tolerance': 1e-8, 'objective_tolerance': 1e-8,
              'witness_reactions': ['Growth', 'CHRPL', 'EX_4hbz_e', 'EX_T4hcinnm_e', '4HBHYOX', 'sink_2ohph_c']}
    aliases = {'PP_5317': 'g_ubic', 'PP_5318': 'g_ubia', 'PP_1376': 'g_transport'}
    return model, recipe, aliases


def test_minimum_import_matches_known_material_requirement_and_restores_context(toy):
    model, recipe, aliases = toy
    wt = model.slim_optimize()
    bounds = {r.id: r.bounds for r in model.reactions}
    for fraction in (.1, .5, .95):
        spec = R.minimum_spec('toy', fraction, wt, 2)
        with model:
            R.apply_case(model, spec, aliases)
            answer = R.checked_solve(model, recipe)
            assert answer['growth'] == pytest.approx(fraction * wt)
            assert -answer['objective'] == pytest.approx(.1 * fraction * wt)
            assert model.reactions.EX_4hbz_e.upper_bound == 0
    assert {r.id: r.bounds for r in model.reactions} == bounds
    assert model.slim_optimize() == pytest.approx(wt)


@pytest.mark.parametrize('blocked', ['UHBZ1t_pp', 'HBZOPT'])
def test_supply_cannot_bypass_transport_or_downstream_synthesis(toy, blocked):
    model, recipe, aliases = toy
    R.apply_case(model, {'loci': ['PP_5317'], 'exchange_bounds': [-2, 1000], 'close_reactions': [blocked]}, aliases)
    assert R.checked_solve(model, recipe)['growth'] == pytest.approx(0)


def test_precursor_rescue_does_not_require_catabolism(toy):
    model, recipe, aliases = toy
    R.apply_case(model, {'loci': ['PP_5317'], 'exchange_bounds': [-2, 1000], 'close_reactions': ['4HBHYOX']}, aliases)
    assert R.checked_solve(model, recipe)['growth'] == pytest.approx(10)


def test_possible_donor_secretion_respects_growth_floor_and_transport(toy):
    model, recipe, aliases = toy
    wt = model.slim_optimize()
    for closure, expected in [([], .5), (['UHBZ1t_pp'], 0)]:
        with model:
            R.apply_case(model, {'objective': 'EX_4hbz_e', 'exchange_bounds': [0, 1000],
                                'growth_bounds': [.95 * wt, 1000], 'close_reactions': closure}, aliases)
            assert R.checked_solve(model, recipe)['objective'] == pytest.approx(expected)


@pytest.mark.parametrize('fraction,wild_type,capacity', [(0, 1, 1), (1, 1, 1), (.5, 0, 1), (.5, math.inf, 1), (.5, 1, math.nan), (.5, 1, -1)])
def test_invalid_targets_are_rejected(fraction, wild_type, capacity):
    with pytest.raises(ValueError):
        R.minimum_spec('toy', fraction, wild_type, capacity)


def test_unknown_case_fields_do_not_silently_change_experiment(toy):
    model, _, aliases = toy
    with pytest.raises(ValueError, match='Unknown'):
        R.apply_case(model, {'add_reaction': 'unapproved'}, aliases)


def test_failed_solver_status_is_never_reported_as_zero(toy, monkeypatch):
    model, recipe, _ = toy
    monkeypatch.setattr(R, 'solve_case', lambda *args: {'status': 'infeasible', 'objective': None})
    with pytest.raises(RuntimeError, match='did not solve optimally'):
        R.checked_solve(model, recipe)


def test_invalid_primal_certificate_is_rejected(toy, monkeypatch):
    model, recipe, _ = toy
    monkeypatch.setattr(R, 'solve_case', lambda *args: {'status': 'optimal', 'objective': 0.,
        'primal_certificate': {'max_abs_S_residual': 1., 'max_bound_violation': 0.}})
    with pytest.raises(RuntimeError, match='Failed numerical certificate'):
        R.checked_solve(model, recipe)


def test_both_audit_provenance_schemas_are_covered():
    record = {'input_sha256': {'identity.py': 'a'},
              'provenance': {'inputs': [{'path': 'fitness.py', 'sha256': 'b'}]}}
    assert R.evidence_hashes(record) == {'identity.py': 'a', 'fitness.py': 'b'}
    with pytest.raises(ValueError, match='no source hashes'):
        R.evidence_hashes({})
    record['provenance']['inputs'].append({'path': 'identity.py', 'sha256': 'different'})
    with pytest.raises(ValueError, match='Conflicting'):
        R.evidence_hashes(record)
