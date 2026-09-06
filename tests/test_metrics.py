import numpy as np
import pytest

from gembench import metrics as M


METRICS = [M.aucpr_bernstein, M.aucpr_standard, M.auroc_standard,
           M.mcc, M.balanced_accuracy, M.accuracy]


@pytest.mark.parametrize("metric", METRICS)
def test_unknown_observations_do_not_become_phenotypes(metric):
    sim = np.array([0.0, 1.0, 0.0, 1.0])
    fit = np.array([-3.0, 0.0, 0.0, -3.0])
    contaminated_sim = np.r_[sim, 0.0, np.nan, np.inf, 1.0]
    contaminated_fit = np.r_[fit, np.nan, 0.0, -3.0, -np.inf]
    assert metric(contaminated_sim, contaminated_fit) == pytest.approx(metric(sim, fit))
    assert sum(M.confusion(contaminated_sim, contaminated_fit).values()) == 4


@pytest.mark.parametrize("metric", METRICS)
def test_no_observations_is_undefined(metric):
    assert np.isnan(metric(np.array([np.nan]), np.array([0.0])))
    assert np.isnan(metric(np.array([]), np.array([])))


def test_threshold_agrees_with_wild_type_protocol():
    assert M.confusion([0.001, 0.0009], [0.0, -3.0]) == dict(tp=1, tn=1, fp=0, fn=0)


def test_shape_mismatch_cannot_silently_realign_observations():
    with pytest.raises(ValueError):
        M.mcc(np.ones((2, 2)), np.ones(4))


def test_bootstrap_empty_gene_intersection():
    assert all(np.isnan(v) for v in M.bootstrap_ci(M.mcc, np.zeros((0, 2)), np.zeros((0, 2)), n_boot=3))
