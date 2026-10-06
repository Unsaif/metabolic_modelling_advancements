"""Synthetic checks of the ranking rules in scripts/iem_disease_ranking.py (no model needed)."""
import importlib.util
import math
import os

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("iem_disease_ranking", os.path.join(ROOT, "scripts", "iem_disease_ranking.py"))
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)


def readout(rxn, pred, h=0.0, d=0.0, st="Optimal"):
    return {"reaction": rxn, "predicted": pred, "healthy": h, "disease": d, "status_healthy": st, "status_disease": st}


def test_rank_stats_ties():
    st = R.rank_stats(np.array([3.0, 1.0, 3.0, 0.0]), 0)
    assert (st["higher"], st["tied"]) == (0, 2)
    assert st["expected_rank"] == 1.5
    assert st["expected_reciprocal_rank"] == pytest.approx(0.75)
    assert st["p_top1"] == 0.5 and st["p_top5"] == 1.0
    st = R.rank_stats(np.array([5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]), 1)
    assert (st["higher"], st["tied"]) == (1, 7)
    assert st["p_top5"] == pytest.approx(4 / 7)


def test_material_call():
    assert R.material_call(0.0, 1.6e-6) == 0          # below the absolute floor
    assert R.material_call(10.0, 10.2) == 0            # 2% change
    assert R.material_call(10.0, 11.0) == 1
    assert R.material_call(11.0, 10.0) == -1
    assert R.material_call(0.0, 0.0) == 0
    assert R.material_call(float("nan"), 1.0) is None


def matrix(spec):
    return [{"iem": name, "call_index": k, "readouts": [readout(r, p) for r, p in rows]}
            for k, (name, rows) in enumerate(spec)]


def test_perfect_specificity_ranks_first():
    n = 10
    names = [f"D{i}" for i in range(n)]
    res = matrix([(nm, [(f"R{j}", "Increased" if j == i else "Unchanged") for j in range(n)]) for i, nm in enumerate(names)])
    order, calls, n_na = R.calls_from_matrix(res, "protocol")
    profiles = {nm: [[f"R{i}", "Increased"]] for i, nm in enumerate(names)}
    out = R.analyse(order, calls, n_na, profiles, draws=2000)
    s = out["summary"]
    assert s["mrr"] == 1.0 and s["expected_top1"] == n
    assert s["null_mrr"] == pytest.approx(sum(1 / i for i in range(1, n + 1)) / n)
    assert s["p_mrr"] < 0.01


def test_identical_predictions_give_the_null():
    n = 8
    names = [f"D{i}" for i in range(n)]
    res = matrix([(nm, [("R0", "Increased"), ("R1", "Decreased")]) for nm in names])
    order, calls, n_na = R.calls_from_matrix(res, "protocol")
    profiles = {nm: [["R0", "Increased"]] for nm in names}
    s = R.analyse(order, calls, n_na, profiles, draws=2000)["summary"]
    assert s["mrr"] == pytest.approx(s["null_mrr"])
    assert s["expected_top1"] == pytest.approx(s["null_top1"])
    assert s["p_mrr"] > 0.5


def test_adjustment_penalises_promiscuous_candidates():
    # D0 changes R0 only; D1 increases everything. Profile of D0: R0 up.
    rows0 = [("R0", "Increased")] + [(f"R{j}", "Unchanged") for j in range(1, 6)]
    rows1 = [(f"R{j}", "Increased") for j in range(6)]
    res = matrix([("D0", rows0), ("D1", rows1)])
    order, calls, n_na = R.calls_from_matrix(res, "protocol")
    profiles = {"D0": [["R0", "Increased"]]}
    plain = R.analyse(order, calls, n_na, profiles, draws=200)["profiles"][0]
    adj = R.analyse(order, calls, n_na, profiles, adjusted=True, panel=[f"R{j}" for j in range(6)], draws=200)["profiles"][0]
    assert plain["tied"] == 2 and plain["expected_rank"] == 1.5
    assert adj["higher"] == 0 and adj["tied"] == 1


def test_na_counts_as_zero_and_is_reported():
    res = matrix([("D0", [("R0", "NA"), ("R1", "Increased")]), ("D1", [("R0", "Increased"), ("R1", "Unchanged")])])
    order, calls, n_na = R.calls_from_matrix(res, "protocol")
    assert n_na == {"D0": 1, "D1": 0}
    names, S = R.score_matrix(order, calls, {"D0": [["R0", "Increased"]]})
    assert list(S[0]) == [0.0, 1.0]


def test_confusers_and_full_lists():
    # D1 always beats D0 on D0's profile; D2 ties with D0 at a positive score.
    res = matrix([("D0", [("R0", "Increased"), ("R1", "Unchanged")]),
                  ("D1", [("R0", "Increased"), ("R1", "Increased")]),
                  ("D2", [("R0", "Increased"), ("R1", "Unchanged")]),
                  ("D3", [("R0", "Unchanged"), ("R1", "Unchanged")])])
    order, calls, n_na = R.calls_from_matrix(res, "protocol")
    out = R.analyse(order, calls, n_na, {"D0": [["R0", "Increased"], ["R1", "Increased"]]}, draws=100)
    prof = out["profiles"][0]
    assert prof["candidates_scoring_higher"] == [["D1", 2.0]] or prof["candidates_scoring_higher"] == [("D1", 2.0)]
    assert prof["candidates_tied"] == ["D2"]
    conf = out["summary"]["confusers_[candidate, profiles_strictly_higher, profiles_tied_with_positive_true_score]"]
    assert conf[0] == ["D1", 1, 0] and ["D2", 0, 1] in conf
