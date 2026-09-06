"""Independent checker regression tests on synthetic models and saved records."""
import copy
import json

import cobra
import numpy as np
import pandas as pd
import pytest

from scripts import verify_ppnp_repair as C


@pytest.fixture
def models(tmp_path):
    base = cobra.Model('toy_parent')
    definitions = {'adn_c': ('C10H13N5O4', 0), 'ade_c': ('C5H5N5', 0),
                   'ins_c': ('C10H12N4O5', 0), 'hxan_c': ('C5H4N4O', 0),
                   'pi_c': ('HO4P', 0), 'r1p_c': ('C5H9O8P', 0),
                   'q8h2_c': ('C49H76O4', 0), 'rib__D_c': ('C5H10O5', 0)}
    for mid, (formula, charge) in definitions.items():
        base.add_metabolites([cobra.Metabolite(mid, name=mid, formula=formula, charge=charge, compartment='C_c')])
    growth = cobra.Reaction('Growth', lower_bound=0, upper_bound=1000)
    growth.add_metabolites({base.metabolites.q8h2_c: -0.001})
    growth.gene_reaction_rule = 'NP_100_1'
    rhcys = cobra.Reaction('RHCYS', lower_bound=0, upper_bound=1000)
    rhcys.add_metabolites({base.metabolites.rib__D_c: 1})
    base.add_reactions([growth, rhcys])
    base.objective = growth
    # The real parent is already an SBML artifact. Reading it once also applies
    # SBML's default gene labels before testing a subsequent export/read cycle.
    cobra.io.write_sbml_model(base, str(tmp_path / 'parent.xml'))
    base = cobra.io.read_sbml_model(str(tmp_path / 'parent.xml'))
    source = cobra.Model('toy_source')
    for rid, coefficients in C.STOICHIOMETRY.items():
        reaction = cobra.Reaction(rid, name=rid, subsystem='salvage', lower_bound=-1000, upper_bound=1000)
        mets = {}
        for mid, value in coefficients.items():
            if mid not in source.metabolites:
                m = base.metabolites.get_by_id(mid).copy()
                if mid in ('pi_c', 'r1p_c'):
                    m.charge = -2
                source.add_metabolites([m])
            mets[source.metabolites.get_by_id(mid)] = value
        reaction.add_metabolites(mets)
        reaction.annotation = {'ec-code': '2.4.2.1'}
        source.add_reactions([reaction])
    cobra.io.write_sbml_model(source, str(tmp_path / 'source.xml'))
    source = cobra.io.read_sbml_model(str(tmp_path / 'source.xml'))
    return base, source


def configuration(arm='both_forward', direction='forward', close=False, substrates=None):
    return {'id': arm, 'substrates': ['adenosine', 'inosine'] if substrates is None else substrates,
            'direction': direction, 'close_rhcys': close}


def exported_arm(base, source, config):
    # A small fixture builder, deliberately not the production intervention function.
    model = base.copy()
    for substrate in config['substrates']:
        r = source.reactions.get_by_id(C.REACTIONS[substrate]).copy()
        r.bounds = (0 if config['direction'] == 'forward' else -1000, 1000)
        r.gene_reaction_rule = 'PP_4248'
        model.add_reactions([r])
    if config['close_rhcys']:
        model.reactions.RHCYS.bounds = (0, 0)
    model.id = base.id + '__ppnp_' + config['id']
    return model


@pytest.mark.parametrize('config', [configuration('parent', substrates=[]),
    configuration('adenosine_forward', substrates=['adenosine']),
    configuration('inosine_forward', substrates=['inosine']), configuration(),
    configuration('parent_rhcys_closed', substrates=[], close=True),
    configuration('both_forward_rhcys_closed', close=True),
    configuration('both_reversible', direction='reversible')])
def test_all_seven_declared_model_variants_and_sbml_roundtrip(models, config, tmp_path, monkeypatch):
    base, source = models
    model = exported_arm(base, source, config)
    path = tmp_path / 'model.xml.gz'
    cobra.io.write_sbml_model(model, str(path))
    loaded = cobra.io.read_sbml_model(str(path))
    def forbidden(*args, **kwargs):
        raise AssertionError('Artifact verification must not optimize')
    monkeypatch.setattr(cobra.Model, 'optimize', forbidden)
    before = {m.id: {r.id for r in m.reactions} for m in base.metabolites}
    assert C.verify_declared_arm(base, loaded, source, config, 'PP_4248') == [C.REACTIONS[s] for s in config['substrates']]
    assert before == {m.id: {r.id for r in m.reactions} for m in base.metabolites}


@pytest.mark.parametrize('change', ['tiny_biomass', 'phosphate', 'reverse_bound', 'wrong_gene',
                                  'extra_reaction', 'metabolite_charge', 'rhcys', 'objective'])
def test_rejects_undeclared_model_changes(models, change):
    base, source = models
    config = configuration('both_forward_rhcys_closed', close=True)
    model = exported_arm(base, source, config)
    if change == 'tiny_biomass':
        model.reactions.Growth.add_metabolites({model.metabolites.q8h2_c: -1e-14})
    elif change == 'phosphate':
        model.reactions.PUNP1.add_metabolites({model.metabolites.pi_c: 1})
    elif change == 'reverse_bound':
        model.reactions.PUNP5.lower_bound = -1000
    elif change == 'wrong_gene':
        model.reactions.PUNP1.gene_reaction_rule = 'NP_100_1'
    elif change == 'extra_reaction':
        model.add_reactions([cobra.Reaction('UNDECLARED')])
    elif change == 'metabolite_charge':
        model.metabolites.pi_c.charge = -2
    elif change == 'rhcys':
        model.reactions.RHCYS.upper_bound = 1000
    else:
        model.objective = model.reactions.PUNP1
    with pytest.raises(C.V.VerificationError, match='differ'):
        C.verify_declared_arm(base, model, source, config, 'PP_4248')


def test_independent_annotation_mapping_adds_ppnp_and_rejects_ambiguous_alias(models):
    base, source = models
    table = pd.DataFrame([{'version': 'NP_100.1', 'locus_tags': 'PP_0001', 'old_locus_tags': ''}])
    model = exported_arm(base, source, configuration())
    assert C.mapping_from_annotations(model, table, {'PP_0001', 'PP_4248'}) == {'NP_100_1': 'PP_0001', 'PP_4248': 'PP_4248'}
    table.loc[0, 'old_locus_tags'] = 'PP_0002'
    with pytest.raises(C.V.VerificationError, match='Ambiguous'):
        C.mapping_from_annotations(model, table, {'PP_0001', 'PP_0002', 'PP_4248'})


def test_raw_replicate_averaging_and_new_gene_coverage_are_separate(tmp_path):
    folder = tmp_path / 'data/fitness_browser/Putida'
    folder.mkdir(parents=True)
    raw = pd.DataFrame({'orgId': ['Putida'] * 2, 'locusId': ['locus1', 'locus2'],
        'sysName': ['PP_0001', 'PP_4248'], 'geneName': ['', 'ppnP'], 'desc': ['', ''],
        'exp1 label': [-4, -1], 'exp2 label': [0, 1]})
    raw.to_csv(folder / 'fit_logratios.tsv', sep='\t', index=False)
    pd.DataFrame({'expName': ['exp1', 'exp2']}).to_csv(folder / 'experiments.tsv', sep='\t', index=False)
    pd.DataFrame({'sysName': ['PP_0001', 'PP_4248', 'unscored']}).to_csv(folder / 'genes.tsv', sep='\t', index=False)
    fitness, _, _ = C.V.read_fitness(tmp_path)
    assert fitness.mean(axis=1).to_dict() == {'PP_0001': -2.0, 'PP_4248': 0.0}
    parent = {'browser_genes': np.array(['PP_0001'])}
    other = {'browser_genes': np.array(['PP_0001', 'PP_4248']), 'sim_growth': np.array([[1, 0], [0, 0]]),
             'fitness': np.array([[-2, 1], [0, np.nan]]), 'wt_growth': np.array([1, 0]), 'st': 0.001}
    coverage = C.coverage_change(parent, other)
    assert coverage == {'new_browser_genes': ['PP_4248'], 'lost_browser_genes': [],
                        'new_gene_coverage': [{'locus': 'PP_4248', 'n_finite_pairs': 1,
                            'n_finite_pairs_where_wt_grows': 1, 'excluded_from_common_gene_comparison': True}]}


def test_reordered_axes_paired_intervals_and_changed_predictions_exclude_new_gene():
    a = {'browser_genes': np.array(['a', 'b']), 'conditions': np.array(['c1', 'c2']),
         'sim_growth': np.array([[1, 1], [1, 1.0]]), 'fitness': np.array([[1, 1], [-3, -3.0]]),
         'wt_growth': np.ones(2), 'st': 0.001, 'ft': -2.0}
    b = copy.deepcopy(a)
    b.update(browser_genes=np.array(['PP_4248', 'b', 'a']),
             sim_growth=np.array([[0, 0], [0, 0], [1, 1.0]]),
             fitness=np.array([[1, 1], [-3, -3], [1, 1.0]]))
    comparison = C.V.paired_comparison(a, b, resamples=20, seed=0)
    assert comparison['n_common_genes'] == 2
    assert comparison['mcc'] == {'A': 0.0, 'B': 1.0}
    assert comparison['paired_mcc_difference_B_minus_A']['point'] == 1
    assert comparison['paired_mcc_difference_B_minus_A']['ci95'] == [0, 1]
    assert {r['gene'] for r in C.V.changed_predictions(a, b, 'variant')} == {'b'}


def test_portable_provenance_accepts_historical_prefix_but_not_wrong_arm(models, tmp_path):
    base, _ = models
    directory = tmp_path / 'parent'
    directory.mkdir()
    cobra.io.write_sbml_model(base, str(directory / 'model.xml.gz'))
    record = {'file': '/historical/checkout/results/parent/model.xml.gz', 'model_id': base.id,
              'n_reactions': len(base.reactions), 'n_metabolites': len(base.metabolites), 'n_genes': len(base.genes),
              'sha256': C.V.sha(directory / 'model.xml.gz')}
    C.V.verify_model_provenance(record, directory, base)
    record['file'] = '/historical/checkout/results/wrong_arm/model.xml.gz'
    with pytest.raises(C.V.VerificationError, match='wrong arm'):
        C.V.verify_model_provenance(record, directory, base)


def test_failure_report_is_saved_and_never_overwritten(tmp_path):
    report = tmp_path / 'failure_report.json'
    assert C.main(['--run-dir', str(tmp_path / 'missing'), '--out', str(report)]) == 1
    saved = json.loads(report.read_text())
    assert not saved['passed'] and saved['failure']['type'] == 'FileNotFoundError'
    assert len(saved['source_sha256']) == 2
    original = report.read_bytes()
    with pytest.raises(C.V.VerificationError, match='already exists'):
        C.main(['--run-dir', str(tmp_path / 'missing'), '--out', str(report)])
    assert report.read_bytes() == original


@pytest.mark.parametrize('corruption', [None, 'replicate_mean', 'ppnp_coverage', 'metric_interval'])
def test_saved_benchmark_detects_averaging_coverage_and_interval_errors(models, tmp_path, corruption):
    base, source = models
    model = exported_arm(base, source, configuration())
    directory = tmp_path / 'both_forward'
    directory.mkdir()
    cobra.io.write_sbml_model(model, str(directory / 'model.xml.gz'))
    model = cobra.io.read_sbml_model(str(directory / 'model.xml.gz'))
    params = {'processes': 2, 'solver': 'glpk', 'growth_threshold': 0.001, 'fitness_threshold': -2.0,
              'carbon_uptake': -10, 'complete_medium_transport': True, 'drop_rich_medium_essentials': False,
              'max_conditions': None, 'genes_subset': None, 'knockout_genes': [], 'rich_medium_uptake': -1000,
              'medium_completion_exclude': ['pnto__R', 'fol', 'hco3']}
    recipe = {'protocol': params, 'caveats': ['synthetic test']}
    arrays = {'sim_growth': np.array([[1.0], [0.0]]), 'fitness': np.array([[-2.0], [-3.0]]),
              'wt_growth': np.array([1.0]), 'model_genes': np.array(['NP_100_1', 'PP_4248']),
              'browser_genes': np.array(['PP_0001', 'PP_4248']), 'conditions': np.array(['glucose | M'])}
    run = {**arrays, 'st': 0.001, 'ft': -2.0}
    native = C.V.native_score(run)
    assert native['mcc']['point'] == 1.0
    counts = {'model_genes': 2, 'model_genes_mapped': 2, 'genes_with_fitness': 2, 'genes_after_adjustment': 2,
              'conditions_total': 1, 'conditions_mapped': 1, 'conditions_wt_grows': 1, 'medium_completion_exchanges_added': 0}
    card = {'protocol': {'params': params, 'intervention': {}, 'study_fingerprint': 'toy', 'evaluation_role': 'development'},
            'warnings': recipe['caveats'], 'dataset_provenance': {'n_genes_with_fitness': '2', 'n_experiments': '2'},
            'model': {'file': '/old/checkout/both_forward/model.xml.gz', 'model_id': model.id,
                'n_genes': 2, 'n_reactions': len(model.reactions), 'n_metabolites': len(model.metabolites),
                'sha256': C.V.sha(directory / 'model.xml.gz')},
            'results': {'condition_level': {'n_conditions_mapped': 1, 'n_conditions_wt_grows': 1},
                        'gene_level_conditions_where_wt_grows': native, 'counts': counts}}
    if corruption == 'replicate_mean':
        arrays['fitness'][0, 0] = -1.9  # Same binary label; only independent raw averaging catches this.
    elif corruption == 'ppnp_coverage':
        for key in ('fitness', 'sim_growth', 'model_genes', 'browser_genes'):
            arrays[key] = arrays[key][:1]
    elif corruption == 'metric_interval':
        card['results']['gene_level_conditions_where_wt_grows']['mcc']['ci95'] = [0.1, 0.2]
    (directory / 'card.json').write_text(json.dumps(C.V.safe(card)))
    np.savez_compressed(directory / 'matrices.npz', **arrays)
    pd.DataFrame([{'condition': 'glucose', 'media': 'M', 'bigg_ids': 'glc__D', 'mapping_confidence': 'high',
                   'n_experiments': 2, 'missing_exchanges': 'EX_glc__D_e', 'wt_growth': 1,
                   'wt_grows': True, 'n_ko_no_growth': 1, 'n_genes_fitness_below_threshold': 1}]).to_csv(
                       directory / 'conditions.tsv', sep='\t', index=False)
    pd.DataFrame([{'condition': 'glucose', 'media': 'M', 'wt_growth': 1, 'n_genes': 2,
                   'tp': 1, 'tn': 1, 'fp': 0, 'fn': 0, 'mcc': 1, 'aucpr_bernstein': 1, 'auroc_standard': 1}]).to_csv(
                       directory / 'per_condition_metrics.tsv', sep='\t', index=False)
    fitness = pd.DataFrame({'exp1': [-4, -4], 'exp2': [0, -2]}, index=['PP_0001', 'PP_4248'])
    experiments = pd.DataFrame({'expGroup': ['carbon source'] * 2, 'condition_1': ['glucose'] * 2,
                                'media': ['M'] * 2, 'condition_2': ['', ''], 'expName': ['exp1', 'exp2']})
    conditions = {'glucose | M': {'bigg_ids': 'glc__D', 'confidence': 'high', 'experiments': ['exp1', 'exp2']}}
    report = {'model': model, 'intervention': {}, 'mapping': {'model_to_browser': {'NP_100_1': 'PP_0001', 'PP_4248': 'PP_4248'}},
              'physical': {'medium_completion_added': [], 'condition': 'glucose | M',
                           'wild_type': {'growth': 1}, 'targeted_gene_predictions': []}}
    arguments = (directory, report, recipe, {'content_fingerprint': 'toy'}, fitness, experiments, conditions)
    if corruption:
        with pytest.raises(C.V.VerificationError, match={'replicate_mean': 'replicate-average',
            'ppnp_coverage': 'model-gene coverage', 'metric_interval': 'ci95'}[corruption]):
            C.V.verify_benchmark(*arguments)
    else:
        _, checked = C.V.verify_benchmark(*arguments)
        assert checked['mcc']['point'] == 1.0 and checked['n_gene_condition_pairs'] == 2
