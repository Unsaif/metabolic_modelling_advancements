import copy
import numpy as np
import pytest
from gembench.comparison import aligned_pairs, compare_runs


def run():
    return dict(benchmark="organism example", params=dict(growth_threshold=0.001, fitness_threshold=-2),
                sim_growth=np.array([[0., 1.], [1., 0.]]), fitness=np.array([[-3., 0.], [0., -3.]]),
                wt_growth=np.ones(2), browser_genes=np.array(["a", "b"]), conditions=np.array(["x", "y"]))


def test_comparison_realigns_order_and_uses_shared_observations():
    a, b = run(), run()
    b["browser_genes"] = b["browser_genes"][::-1]
    b["sim_growth"] = b["sim_growth"][::-1].copy()
    b["fitness"] = b["fitness"][::-1].copy()
    b["sim_growth"][0, 0] = np.nan
    sa, sb, fit, *_ = aligned_pairs(a, b)
    assert np.isfinite(fit).sum() == 3
    assert np.isnan(fit[1, 0])
    report = compare_runs(a, b, n_boot=10)
    assert report["mcc"]["A"] == report["mcc"]["B"]


def test_changed_fitness_is_rejected():
    a, b = run(), run()
    b["fitness"][0, 0] = -1
    with pytest.raises(ValueError, match="fitness differs"):
        aligned_pairs(a, b)


def test_condition_coverage_not_hidden_in_score():
    a, b = run(), run()
    a["wt_growth"][1] = 0
    report = compare_runs(a, b, n_boot=5)
    assert report["wt_grows"] == {"A": 1, "B": 2}
    assert report["n_conditions_both_grow"] == 1


def test_no_overlap_has_no_performance_estimate():
    a, b = run(), run()
    b["browser_genes"] = np.array(["c", "d"])
    report = compare_runs(a, b, n_boot=5)
    assert np.isnan(report["mcc"]["A"])
    assert report["n_shared_finite_pairs"] == 0


def test_different_thresholds_are_rejected():
    a, b = run(), run()
    b["params"]["growth_threshold"] = 0.01
    with pytest.raises(ValueError, match="growth_threshold"):
        aligned_pairs(a, b)
