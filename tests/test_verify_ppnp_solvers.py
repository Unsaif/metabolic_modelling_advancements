"""Synthetic cases for the new solver checker; no real outcome values used."""
import copy
import json

import cobra
import pytest

from scripts import verify_ppnp_solvers as V


@pytest.fixture
def case():
    model = cobra.Model('synthetic_ppnp')
    nutrient = cobra.Metabolite('a_c', compartment='c')
    q = cobra.Metabolite('q8_c', compartment='c')
    qh = cobra.Metabolite('q8h2_c', compartment='c')
    supply = cobra.Reaction('SUPPLY', lower_bound=0, upper_bound=2)
    supply.add_metabolites({nutrient: 1})
    synthesis = cobra.Reaction('DMQMT', lower_bound=0, upper_bound=1000)
    synthesis.add_metabolites({q: 1})
    redox = cobra.Reaction('REDOX', lower_bound=-1000, upper_bound=1000)
    redox.add_metabolites({q: -1, qh: 1})
    growth = cobra.Reaction('Growth', lower_bound=0, upper_bound=1000)
    growth.add_metabolites({nutrient: -1, qh: -0.5})
    growth.gene_reaction_rule = 'rbk'
    model.add_reactions([supply, synthesis, redox, growth])
    model.objective = growth
    recipe = {'targeted_gene_loci': ['PP_2458', 'PP_4248'],
              'physical_settings': {'residual_tolerance': 1e-8, 'production_capacity': 1000,
                                    'production_threshold': 1e-6},
              'protocol': {'growth_threshold': 0.001},
              'witness_reactions': ['Growth', 'DMQMT'],
              'pool_weights': {'q8_c': 1, 'q8h2_c': 1},
              'mechanistic_cases': [{'id': 'wild_type', 'loci': []}]}
    row = {'id': 'wild_type', 'loci': [], 'close_added_reactions': False,
           'status': 'optimal', 'growth': 2.0, 'mass_balance_max_abs': 0.,
           'bound_violation_max': 0., 'selected_fluxes': {'Growth': 2.0, 'DMQMT': 1.0}}
    mapping = {'org_id': 'Putida', 'model_to_browser': {'rbk': 'PP_2458'}}
    return model, recipe, row, mapping


def test_aggregate_pool_is_checked_when_redox_is_not_in_partial_witness(case):
    model, recipe, row, _ = case
    report = V.selected_flux_checks(model, row, recipe)
    assert report['aggregate_quinone_pool_residual'] == 0
    assert report['complete_metabolite_balance_residuals'] == {}
    row['selected_fluxes']['DMQMT'] = 1.01
    with pytest.raises(ValueError, match='aggregate quinone'):
        V.selected_flux_checks(model, row, recipe)


def test_complete_selected_balance_catches_nonquinone_failure(case):
    model, recipe, row, _ = case
    recipe['witness_reactions'].append('SUPPLY')
    row['selected_fluxes']['SUPPLY'] = 1.0
    with pytest.raises(ValueError, match='completely covered'):
        V.selected_flux_checks(model, row, recipe)


@pytest.mark.parametrize('corruption', ['missing', 'extra', 'nan', 'growth', 'claimed_residual', 'bounds'])
def test_corrupt_partial_witness_rejected(case, corruption):
    model, recipe, row, _ = case
    if corruption == 'missing':
        del row['selected_fluxes']['DMQMT']
    elif corruption == 'extra':
        row['selected_fluxes']['REDOX'] = 1.
    elif corruption == 'nan':
        row['selected_fluxes']['DMQMT'] = float('nan')
    elif corruption == 'growth':
        row['growth'] = 1.
    elif corruption == 'claimed_residual':
        row['mass_balance_max_abs'] = 1e-3
    else:
        model.reactions.Growth.upper_bound = 1.
    with pytest.raises(ValueError):
        V.selected_flux_checks(model, row, recipe)


def test_aliases_allow_true_absence_but_reject_collisions_and_ambiguity(case):
    model, recipe, _, mapping = case
    assert V.unique_aliases(model, mapping, recipe['targeted_gene_loci']) == {'PP_2458': 'rbk', 'PP_4248': 'PP_4248'}
    mapping['model_to_browser']['another'] = 'PP_2458'
    with pytest.raises(ValueError, match='Ambiguous'):
        V.unique_aliases(model, mapping, recipe['targeted_gene_loci'])
    del mapping['model_to_browser']['another']
    model.reactions.DMQMT.gene_reaction_rule = 'PP_4248'
    with pytest.raises(ValueError, match='collision'):
        V.unique_aliases(model, mapping, recipe['targeted_gene_loci'])


def test_missing_case_cannot_carry_fake_no_growth_prediction(case):
    model, recipe, row, mapping = case
    specification = {'id': 'rbk_ppnp', 'loci': ['PP_2458', 'PP_4248']}
    recipe['mechanistic_cases'] = [specification]
    absent = {'id': 'rbk_ppnp', 'loci': specification['loci'], 'close_added_reactions': False,
              'status': 'unrepresented_gene_or_function', 'growth': None, 'unrepresented': ['PP_4248']}
    _, checked = V.mechanism_inventory(model, mapping, {'cases': [absent]}, recipe)
    assert checked[0][2] == ['PP_4248']
    absent['growth'] = 0.
    with pytest.raises(ValueError, match='absent/orphan'):
        V.mechanism_inventory(model, mapping, {'cases': [absent]}, recipe)


def test_mechanism_closures_restore_parent_context(case):
    model, _, _, _ = case
    addition = cobra.Reaction('PUNP1', lower_bound=-1000, upper_bound=1000)
    model.add_reactions([addition])
    before = {r.id: r.bounds for r in model.reactions}
    with model:
        V.apply_mechanism(model, {'PP_2458': 'rbk'}, {'loci': ['PP_2458'], 'close_added_reactions': True})
        assert model.reactions.Growth.bounds == (0, 0)
        assert model.reactions.PUNP1.bounds == (0, 0)
    assert {r.id: r.bounds for r in model.reactions} == before


def test_numerical_mechanism_checks_do_not_require_identical_degenerate_fluxes(case):
    model, recipe, row, mapping = case
    # Disconnected cycle has infinitely many optimal feasible flux choices.
    b, c = (cobra.Metabolite(name, compartment='c') for name in ('b_c', 'c_c'))
    for rid, stoichiometry in [('CYCLE1', {b: -1, c: 1}), ('CYCLE2', {c: -1, b: 1})]:
        reaction = cobra.Reaction(rid, lower_bound=0, upper_bound=1000)
        reaction.add_metabolites(stoichiometry)
        model.add_reactions([reaction])
        recipe['witness_reactions'].append(rid)
        row['selected_fluxes'][rid] = 7.
    results = []
    V.check_mechanisms(model, mapping, {'cases': [row]}, recipe, lambda name, result: results.append(result))
    assert len(results) == 1
    assert results[0]['valid'], results
    assert results[0]['independent_glpk_selected_fluxes']['CYCLE1'] != 7.


def test_failed_solver_summary_preserves_failure_as_json_without_infinity():
    case = {'valid': False, 'solves': [{'status': 'infeasible', 'objective': None},
                                     {'status': 'infeasible', 'objective': None}]}
    report = V.summarize([case], 1, 'synthetic')
    assert not report['valid']
    assert report['failed_cases'] == [case]
    assert report['n_solver_runs_without_full_finite_certificate'] == 2
    assert report['max_independent_mass_balance_residual'] is None
    json.dumps(report, allow_nan=False)


def test_failure_report_does_not_overwrite_existing_outputs(tmp_path):
    study = tmp_path / 'study'
    study.mkdir()
    (study / 'manifest.json').write_text('{}')
    for arm in V.ARM_IDS:
        directory = study / arm
        directory.mkdir()
        for name in ('model.xml.gz', 'physical.json', 'gene_map.json', 'intervention.json', 'mechanism.json'):
            (directory / name).write_text('{}')
    out = tmp_path / 'verification'
    with pytest.raises(KeyError):
        V.main(['--study', str(study), '--out', str(out)])
    failure = (out / 'failure.json').read_bytes()
    assert json.loads(failure)['completed_cases'] == 0
    with pytest.raises(FileExistsError):
        V.main(['--study', str(study), '--out', str(out)])
    assert (out / 'failure.json').read_bytes() == failure
