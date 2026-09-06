import json

import cobra
import numpy as np
import pytest
from sklearn.metrics import auc, average_precision_score, matthews_corrcoef, precision_recall_curve, roc_auc_score

from scripts import verify_quinone_repair as verifier


def run(sim, fit, wt, genes=None, conditions=None):
    sim = np.asarray(sim, dtype=float)
    return {'sim_growth': sim, 'fitness': np.asarray(fit, dtype=float), 'wt_growth': np.asarray(wt, dtype=float),
            'browser_genes': np.asarray(genes or [f'g{i}' for i in range(sim.shape[0])]),
            'conditions': np.asarray(conditions or [f'c{i}' for i in range(sim.shape[1])]), 'st': 0.001, 'ft': -2.0}


def test_elementary_confusion_masks_before_thresholding_and_includes_equality():
    sim = np.array([0.001, 0, 1, 0, np.nan, np.inf, 0])
    fit = np.array([-2, -3, -3, -1, -3, -3, -np.inf])
    counts = verifier.confusion(sim, fit)
    assert counts == {'tp': 1, 'tn': 1, 'fp': 1, 'fn': 1}
    assert verifier.mcc_from_counts(counts) == 0


def test_no_observations_undefined_but_nonempty_degenerate_mcc_zero():
    assert np.isnan(verifier.mcc_from_counts({'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0}))
    assert verifier.mcc_from_counts({'tp': 4, 'tn': 0, 'fp': 0, 'fn': 0}) == 0
    assert verifier.mcc_from_counts({'tp': 3, 'tn': 4, 'fp': 0, 'fn': 0}) == 1
    assert verifier.mcc_from_counts({'tp': 0, 'tn': 0, 'fp': 3, 'fn': 4}) == -1


def test_tie_grouped_rank_metrics_and_elementary_mcc_against_external_reference():
    rng = np.random.default_rng(9)
    # Ties and finite masks exercise details that easy perfect-ranking examples miss.
    sim = rng.choice([0, 0.001, 0.2, 1, np.nan], size=(17, 4))
    fit = rng.choice([-4, -2, 0, 1, np.nan], size=sim.shape)
    points = verifier.metric_points(sim, fit, 0.001, -2)
    finite = np.isfinite(sim) & np.isfinite(fit)
    s, f = sim[finite], fit[finite]
    precision, recall, _ = precision_recall_curve(s >= 0.001, -f, pos_label=0)
    assert points['aucpr_bernstein'] == pytest.approx(auc(recall, precision), abs=1e-14)
    assert points['aucpr_standard'] == pytest.approx(average_precision_score(f < -2, -s), abs=1e-14)
    assert points['auroc_standard'] == pytest.approx(roc_auc_score(f < -2, -s), abs=1e-14)
    assert points['mcc'] == pytest.approx(matthews_corrcoef(f >= -2, s >= 0.001), abs=1e-14)


def test_native_no_growing_arm_retains_zero_coverage_without_fabricating_score():
    result = verifier.native_score(run([[0, 0], [0, 0]], [[-3, 1], [1, -3]], [0, 0]), resamples=8)
    assert result == {'n_genes': 0, 'n_conditions': 0, 'n_gene_condition_pairs': 0,
                      'n_missing_fitness': 0, 'n_nonfinite_simulation': 0}
    assert 'mcc' not in result


def test_native_coverage_masks_missing_rows_and_retains_native_gene_axis():
    result = verifier.native_score(run([[1, 0], [0, 0], [np.nan, 0]], [[1, -3], [np.nan, 1], [-3, 1]], [1, 0]), resamples=8)
    assert result['n_genes'] == 2  # Has observed fitness, even if one prediction is unknown.
    assert result['n_conditions'] == 1 and result['n_gene_condition_pairs'] == 1
    assert result['n_nonfinite_simulation'] == 1
    assert result['confusion'] == {'tp': 1, 'tn': 0, 'fp': 0, 'fn': 0}


def test_alignment_uses_common_genes_finite_both_wt_and_identical_observations():
    a = run([[1, 0, 0], [0, 1, 0]], [[1, -3, 1], [-3, 1, 1]], [1, 1, np.inf], ['a', 'b'])
    b = run([[1, 1, 0], [0, 0, 0]], [[-3, 1, 1], [1, -3, 1]], [1, 0, 1], ['b', 'extra'])
    result = verifier.paired_comparison(a, b, resamples=8)
    assert result['n_common_genes'] == 1 and result['conditions_both_grow'] == ['c0']
    assert result['n_common_conditions'] == 3 and result['n_shared_finite_pairs'] == 1
    assert result['confusion']['A'] == {'tp': 0, 'tn': 1, 'fp': 0, 'fn': 0}
    assert result['confusion']['B'] == {'tp': 0, 'tn': 0, 'fp': 1, 'fn': 0}
    assert verifier.changed_predictions(a, b, 'other') == [
        {'arm': 'other', 'gene': 'b', 'condition': 'c0', 'before': 0.0, 'after': 1.0, 'fitness': -3.0}]
    b['fitness'][0, 0] = -4
    with pytest.raises(verifier.VerificationError, match='observations disagree'):
        verifier.align(a, b)


def test_empty_alignment_and_nonfinite_simulations_cannot_be_scored_as_no_growth():
    a = run([[0, 1]], [[-3, 1]], [0, 1])
    b = run([[1, 0]], [[-3, 1]], [1, 0])
    compared = verifier.paired_comparison(a, b, resamples=8)
    assert compared['n_conditions_both_grow'] == 0 and compared['n_shared_finite_pairs'] == 0
    assert np.isnan(compared['mcc']['A']) and np.isnan(compared['mcc']['B'])
    assert compared['paired_mcc_difference_B_minus_A']['finite_bootstrap_samples'] == 0
    assert verifier.safe(compared)['mcc'] == {'A': None, 'B': None}
    a['wt_growth'][:] = b['wt_growth'][:] = 1
    b['sim_growth'][:] = np.nan
    assert verifier.paired_comparison(a, b, resamples=8)['n_shared_finite_pairs'] == 0
    assert verifier.changed_predictions(a, b, 'other') == []


def test_pool_certificate_detects_every_stored_nonzero_without_rounding():
    model = cobra.Model('pool')
    a, b = (cobra.Metabolite(name, formula='C', charge=0, compartment='c') for name in ['a_c', 'b_c'])
    cycle = cobra.Reaction('CYCLE', lower_bound=-100, upper_bound=100)
    cycle.add_metabolites({a: -1, b: 1})
    demand = cobra.Reaction('Growth', lower_bound=0, upper_bound=100)
    demand.add_metabolites({b: -1e-20})
    model.add_reactions([cycle, demand])
    result = verifier.pool_certificate(model, {'a_c': 1, 'b_c': 1})
    assert not result['exactly_conserved_in_stoichiometry']
    assert result['direction_compatible_supply_reactions'] == []
    assert result['reactions'][1]['exact_numerator'] != '0'


@pytest.mark.parametrize('field,value', [('lower_bound', -100), ('upper_bound', 1),
                                        ('metabolites', {'a_c': -1, 'b_c': 2}),
                                        ('gene_reaction_rule', 'PP_1 and PP_2')])
def test_reaction_comparison_rejects_undeclared_chemistry_direction_and_gpr(field, value):
    expected = {'id': 'reaction', 'lower_bound': 0, 'upper_bound': 1000,
                'metabolites': {'a_c': -1, 'b_c': 1}, 'gene_reaction_rule': 'PP_1 or PP_2'}
    actual = {**expected, field: value}
    with pytest.raises(verifier.VerificationError):
        verifier.compare_reaction(actual, expected, 'candidate')


def test_gpr_serialization_order_is_not_mistaken_for_a_different_rule():
    verifier.compare_reaction({'gene_reaction_rule': 'PP_2 or PP_1', 'equation': 'serialized'},
                              {'gene_reaction_rule': 'PP_1 or PP_2', 'equation': 'source'}, 'rule')
    assert verifier.renamed_rule('PP_1 and PP_11', {'PP_1': 'NP_1_1'}) == 'NP_1_1 and PP_11'


@pytest.fixture
def relocated_model(tmp_path):
    directory = tmp_path / 'another_checkout' / 'results' / 'curated_path_template'
    directory.mkdir(parents=True)
    model = cobra.Model('relocated_example')
    met = cobra.Metabolite('a_c', formula='C', charge=0, compartment='c')
    reaction = cobra.Reaction('Growth')
    reaction.add_metabolites({met: -1})
    model.add_reactions([reaction])
    model.objective = reaction
    cobra.io.write_sbml_model(model, str(directory / 'model.xml.gz'))
    card = {'file': '/original/machine/checkout/results/curated_path_template/model.xml.gz',
            'model_id': model.id, 'sha256': verifier.sha(directory / 'model.xml.gz'),
            'n_reactions': len(model.reactions), 'n_metabolites': len(model.metabolites), 'n_genes': len(model.genes)}
    return directory, model, card


def test_model_provenance_accepts_relocated_checkout_with_identical_arm_and_bytes(relocated_model):
    directory, model, card = relocated_model
    verifier.verify_model_provenance(card, directory, model)
    card['file'] = 'original_relative_checkout/curated_path_template/model.xml.gz'
    verifier.verify_model_provenance(card, directory, model)


@pytest.mark.parametrize('field,value,match', [
    ('file', '/old/checkout/baseline/model.xml.gz', 'wrong arm'),
    ('file', '/old/checkout/curated_path_template/another.xml.gz', 'wrong arm'),
    ('sha256', '0' * 64, 'sha256'),
    ('model_id', 'another_model', 'model_id'),
])
def test_relocation_never_weakens_arm_hash_or_model_identity(relocated_model, field, value, match):
    directory, model, card = relocated_model
    card[field] = value
    with pytest.raises(verifier.VerificationError, match=match):
        verifier.verify_model_provenance(card, directory, model)


def test_failed_verification_writes_hashed_report_and_exits_nonzero(tmp_path, monkeypatch):
    (tmp_path / 'artifact.txt').write_text('immutable result')
    monkeypatch.setattr(verifier, 'verify', lambda *_: verifier.require(False, 'deliberate mismatch'))
    monkeypatch.setattr('sys.argv', ['verify_quinone_repair', '--run-dir', str(tmp_path)])
    with pytest.raises(SystemExit) as result:
        verifier.main()
    assert result.value.code == 1
    report = json.loads((tmp_path / 'independent_verification.json').read_text())
    assert not report['passed'] and report['failure']['message'] == 'deliberate mismatch'
    assert report['result_input_sha256']['artifact.txt'] == verifier.sha(tmp_path / 'artifact.txt')
    assert report['verification_script']['sha256'] == verifier.sha(verifier.__file__)
    with pytest.raises(verifier.VerificationError, match='Preserve the existing report'):
        verifier.main()
