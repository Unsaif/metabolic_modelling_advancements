import json
from types import SimpleNamespace

import cobra
import pytest

from gembench.cofactor import probe_metabolite_production
from gembench.gene_mapping import GeneMap
from scripts import run_quinone_repair as runner


@pytest.fixture
def toy():
    source = cobra.Model('source')
    mids = ['a_c', 'b_c', 'c_c', 'd_c', 'e_c', 'q8h2_c']
    mets = {mid: cobra.Metabolite(mid, formula='C', charge=0, compartment='c') for mid in mids}
    rules = ['PP_1765', 'PP_5199', 'PP_5011', 'PP_0427 or PP_5197', 'PP_1765']
    for i, (rid, rule) in enumerate(zip(runner.CANDIDATES, rules)):
        r = cobra.Reaction(rid, lower_bound=0, upper_bound=999999)
        r.add_metabolites({mets[mids[i]]: -1, mets[mids[i + 1]]: 1})
        r.gene_reaction_rule = rule
        source.add_reactions([r])
    base = cobra.Model('base')
    base.add_metabolites([mets['a_c'].copy(), mets['q8h2_c'].copy(),
                         cobra.Metabolite('mql8_c', formula='C', charge=0, compartment='c')])
    for rid, stoich, rule in [('FOOD', {'a_c': 1}, ''), ('Growth', {'a_c': -1}, ''),
                             ('OTHER_UBIE_ROLE', {'a_c': -1, 'mql8_c': 1}, 'NP_747113_1')]:
        r = cobra.Reaction(rid, lower_bound=0, upper_bound=10)
        r.add_metabolites({base.metabolites.get_by_id(mid): value for mid, value in stoich.items()})
        r.gene_reaction_rule = rule
        base.add_reactions([r])
    base.objective = 'Growth'
    mapping = GeneMap('Putida', {'NP_747113_1': 'PP_5011'}, [])
    aliases = {g.id: 'NP_747113_1' if g.id == 'PP_5011' else g.id for g in source.genes}
    candidate = {'reaction_ids': runner.CANDIDATES.copy(),
                 'reaction_records': [runner.reaction_record(r) for r in source.reactions],
                 'all_involved_metabolites': {m.id: runner.met_record(m) for m in source.metabolites},
                 'new_metabolites': {m.id: runner.met_record(m) for m in source.metabolites if m.id not in base.metabolites},
                 'gene_aliases': aliases,
                 'gene_identity_records': [{'source_locus': k, 'target_model_gene_id': v,
                                             'identity_mapping_ambiguous': False} for k, v in aliases.items()],
                 'uncertainties': ['Provisional toy chemistry.']}
    return base, source, candidate, mapping


def mapped(model):
    return GeneMap('Putida', {g.id: 'PP_5011' if g.id == 'NP_747113_1' else g.id for g in model.genes}, [])


def test_exact_transfer_reuses_existing_gene_and_enables_production(toy):
    model, source, candidate, mapping = toy
    assert probe_metabolite_production(model, ['q8h2_c'])[0]['maximum_demand_flux'] == 0
    before_source = [runner.reaction_record(r) for r in source.reactions]
    report = runner.transfer_candidate(model, source, candidate, mapping)
    assert len(report['source_reaction_ids']) == 5
    assert report['new_metabolites'] == ['b_c', 'c_c', 'd_c', 'e_c']
    assert 'PP_5011' not in model.genes
    assert model.reactions.OMBZLM.gene_reaction_rule == 'NP_747113_1'
    assert len(model.genes.NP_747113_1.reactions) == 2
    assert [runner.reaction_record(r) for r in source.reactions] == before_source
    for rid in runner.CANDIDATES:
        assert model.reactions.get_by_id(rid).bounds == source.reactions.get_by_id(rid).bounds
        assert {m.id: c for m, c in model.reactions.get_by_id(rid).metabolites.items()} == candidate['reaction_records'][runner.CANDIDATES.index(rid)]['metabolites']
    assert probe_metabolite_production(model, ['q8h2_c'])[0]['maximum_demand_flux'] == pytest.approx(10)


@pytest.mark.parametrize('error', ['existing_reaction', 'source_record', 'metabolite_conflict', 'aliases', 'ambiguous_identity'])
def test_transfer_rejects_conflicts_before_any_mutation(toy, error):
    model, source, candidate, mapping = toy
    if error == 'existing_reaction':
        model.add_reactions([cobra.Reaction(runner.CANDIDATES[0])])
    elif error == 'source_record':
        candidate['reaction_records'][0]['upper_bound'] = 1000
    elif error == 'metabolite_conflict':
        model.metabolites.a_c.formula = 'H'
    elif error == 'aliases':
        candidate['gene_aliases']['PP_5011'] = 'PP_5011'
    else:
        candidate['gene_identity_records'][0]['identity_mapping_ambiguous'] = True
    before = ({r.id for r in model.reactions}, {m.id for m in model.metabolites}, {g.id for g in model.genes})
    with pytest.raises(ValueError):
        runner.transfer_candidate(model, source, candidate, mapping)
    assert ({r.id for r in model.reactions}, {m.id for m in model.metabolites}, {g.id for g in model.genes}) == before


def test_ambiguous_existing_locus_is_not_silently_duplicated(toy):
    model, _, _, mapping = toy
    mapping.model_to_browser['another_id'] = 'PP_5011'
    with pytest.raises(ValueError, match='Ambiguous model gene identity'):
        runner.resolve_aliases(model, mapping, ['PP_5011'])


def test_explicit_metadata_policy_preserves_target_and_normalizes_only_new_compartments(toy):
    model, source, candidate, mapping = toy
    for met in source.metabolites:
        met.charge = 1  # Every one-for-one source reaction remains charge balanced.
    for met in model.metabolites:
        met.compartment = 'C_c'
    candidate['reaction_records'] = [runner.reaction_record(r) for r in source.reactions]
    candidate['all_involved_metabolites'] = {m.id: runner.met_record(m) for m in source.metabolites}
    candidate['new_metabolites'] = {m.id: runner.met_record(m) for m in source.metabolites if m.id not in model.metabolites}
    policy = {'compartment_aliases': {'c': 'C_c'},
              'expected_shared_charge_differences': {'a_c': {'source': 1, 'target': 0}, 'q8h2_c': {'source': 1, 'target': 0}}}
    with pytest.raises(ValueError, match='formula/compartment conflict'):
        runner.transfer_candidate(model, source, candidate, mapping)
    with pytest.raises(ValueError, match='Shared charge differences'):
        runner.transfer_candidate(model, source, candidate, mapping, {'compartment_aliases': {'c': 'C_c'}})
    report = runner.transfer_candidate(model, source, candidate, mapping, policy)
    assert model.metabolites.a_c.charge == 0 and model.metabolites.q8h2_c.charge == 0
    assert model.metabolites.b_c.charge == 1 and model.metabolites.b_c.compartment == 'C_c'
    assert source.metabolites.b_c.compartment == 'c'
    assert all(not value for value in report['source_balance'].values())
    assert report['applied_reactions'][0]['mass_charge_imbalance'] == {'charge': 1}
    assert report['preserved_charge_differences'] == policy['expected_shared_charge_differences']


def test_sequence_skip_requires_same_lp_and_only_unscored_changed_genes(toy):
    model, source, candidate, mapping = toy
    runner.transfer_candidate(model, source, candidate, mapping)
    variant = model.copy()
    variant.reactions.OMMBLHX.gene_reaction_rule = 'PP_0427'
    report = runner.verify_sequence_equivalence(model, variant, mapped(model), mapped(variant), {'PP_5317'})
    assert report['benchmark_omitted_for_sequence_arm']
    assert report['changed_reactions'] == [{'reaction': 'OMMBLHX', 'loci': ['PP_0427', 'PP_5197']}]
    with pytest.raises(ValueError, match='benchmark-covered genes'):
        runner.verify_sequence_equivalence(model, variant, mapped(model), mapped(variant), {'PP_5197'})
    variant.reactions.FOOD.upper_bound = 5
    with pytest.raises(ValueError, match='stoichiometry, bounds'):
        runner.verify_sequence_equivalence(model, variant, mapped(model), mapped(variant), set())


def test_targeted_knockouts_distinguish_absent_orphan_and_represented_genes(toy):
    model, source, candidate, mapping = toy
    runner.transfer_candidate(model, source, candidate, mapping)
    model.reactions.Growth.add_metabolites({model.metabolites.q8h2_c: -0.001})
    model.reactions.OMMBLHX.gene_reaction_rule = 'PP_0427'  # PP_5197 remains as an orphan.
    assert 'PP_5197' in model.genes and not model.genes.PP_5197.reactions
    before = {r.id: r.bounds for r in model.reactions}
    results = runner.targeted_knockouts(model, mapped(model), ['absent', 'PP_5197', 'PP_1765'], set())
    assert results[0]['status'] == 'unrepresented_gene' and results[0]['growth'] is None
    assert results[1]['status'] == 'unrepresented_function' and results[1]['growth'] is None
    assert results[2]['status'] == 'optimal' and results[2]['growth'] == pytest.approx(0)
    assert {r.id: r.bounds for r in model.reactions} == before
    assert all(r['role'].startswith('model prediction only') for r in results)


def test_targeted_zero_is_uninterpretable_when_wild_type_cannot_grow(toy):
    model, _, _, mapping = toy
    model.reactions.FOOD.upper_bound = 0
    result = runner.targeted_knockouts(model, mapping, ['PP_5011'], set())[0]
    assert result['status'] == 'optimal' and result['growth'] == pytest.approx(0)
    assert result['wild_type_growth'] == 0 and not result['wild_type_grows_at_threshold']
    assert 'Uninterpretable for gene essentiality' in result['interpretation']


def test_configuration_changes_only_declared_biomass_coefficient(toy, monkeypatch):
    base, source, candidate, mapping = toy
    monkeypatch.setattr(runner, 'gene_map_for', lambda org, model, names, variant: mapped(model))
    config = {'id': 'demand_only', 'add_pathway': False, 'biomass_coefficient_source': 'template', 'gpr_overrides': {}}
    model, _, report = runner.apply_configuration(base, source, candidate, mapping, config,
                                                  {'template': -0.0001}, set(candidate['gene_aliases']))
    assert report['applied_biomass_coefficient'] == -0.0001
    assert {m.id: c for m, c in model.reactions.Growth.metabolites.items()} == {'a_c': -1, 'q8h2_c': -0.0001}
    assert {m.id: c for m, c in base.reactions.Growth.metabolites.items()} == {'a_c': -1}
    with pytest.raises(ValueError, match='no biomass quinone demand'):
        runner.apply_configuration(model, source, candidate, mapping, config, {'template': -0.0001}, set())


@pytest.mark.parametrize('status,objective', [('infeasible', None), ('time_limit', 0), ('optimal', float('nan'))])
def test_failed_growth_is_never_a_zero_prediction(toy, monkeypatch, status, objective):
    model = toy[0]
    monkeypatch.setattr(model, 'optimize', lambda: SimpleNamespace(status=status, objective_value=objective))
    with pytest.raises(RuntimeError, match='Growth solve failed'):
        runner.growth_solution(model)


def test_existing_output_is_rejected_before_loading_inputs(tmp_path, monkeypatch):
    out = tmp_path / 'already_completed'
    out.mkdir()
    sentinel = out / 'sentinel'
    sentinel.write_text('preserve')
    monkeypatch.setattr('sys.argv', ['run_quinone_repair', '--out', str(out), '--plan', str(tmp_path / 'missing.json')])
    with pytest.raises(FileExistsError, match='Output exists'):
        runner.main()
    assert list(out.iterdir()) == [sentinel]
    assert sentinel.read_text() == 'preserve'


def test_fresh_failure_is_preserved(tmp_path, monkeypatch):
    out = tmp_path / 'new'
    monkeypatch.setattr('sys.argv', ['run_quinone_repair', '--out', str(out), '--plan', str(tmp_path / 'missing.json')])
    with pytest.raises(FileNotFoundError):
        runner.main()
    assert json.loads((out / 'failure.json').read_text())['type'] == 'FileNotFoundError'


def test_manifest_is_written_before_any_source_model_is_loaded(tmp_path, monkeypatch):
    out = tmp_path / 'fresh'
    source = tmp_path / 'source.xml'
    source.write_text('toy source bytes')
    candidate = {'provenance': {'input_sha256': {}},
                 'source_model': {'path': 'source.xml', 'sha256': runner.sha256_of(source)}}
    (tmp_path / 'candidate.json').write_text(json.dumps(candidate))
    recipe = {'evaluation_role': 'development', 'organism': 'Putida', 'base_model_patch_version': '0.4',
              'configurations': [{'id': name} for name in runner.ARM_IDS], 'candidate_path': 'candidate.json'}
    plan = tmp_path / 'plan.json'
    plan.write_text(json.dumps(recipe))
    monkeypatch.setattr(runner, 'ROOT', tmp_path)
    monkeypatch.setattr(runner, 'freeze_inputs', lambda *args: {'content_fingerprint': 'test frozen inputs'})
    monkeypatch.setattr(runner, 'verify_inputs', lambda *args: {'valid': True})

    def first_model_load(*args):
        assert json.loads((out / 'manifest.json').read_text())['content_fingerprint'] == 'test frozen inputs'
        raise RuntimeError('Stop after checking freeze order')

    monkeypatch.setattr(cobra.io, 'read_sbml_model', first_model_load)
    monkeypatch.setattr('sys.argv', ['run_quinone_repair', '--plan', str(plan), '--out', str(out)])
    with pytest.raises(RuntimeError, match='Stop after checking freeze order'):
        runner.main()
    assert (out / 'failure.json').exists()
