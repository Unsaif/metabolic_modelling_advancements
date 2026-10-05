import importlib.util
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("iem_range_calls", os.path.join(ROOT, "scripts", "iem_range_calls.py"))
rc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rc)


@pytest.mark.parametrize("h_min,h_max,d_min,d_max,expected", [
    (0, 10, 0, 10, "Unchanged"),
    (0, 10, 5, 10, "Increased"),      # equal maxima at a cap, forced minimum rises
    (0, 10, 0, 12, "Increased"),
    (0, 10, 0, 5, "Decreased"),
    (2, 10, 0, 10, "Decreased"),
    (0, 10, 5, 5, "Conflicting"),     # minimum up, maximum down
    (0, 10, 5e-7, 10 + 5e-7, "Unchanged"),
])
def test_full_range_rule(h_min, h_max, d_min, d_max, expected):
    assert rc.range_call(h_min, h_max, d_min, d_max) == expected


@pytest.mark.parametrize("h_min,h_max,d_min,d_max,expected", [
    (0, 10, 5, 10, "Increased"),
    (3, 10, 1, 10, "Decreased"),
    (0, 10, 5, 8, "Decreased"),       # the maximum call is kept whenever it is not Unchanged
    (0, 10, 0, 10, "Unchanged"),
])
def test_tie_break_rule(h_min, h_max, d_min, d_max, expected):
    assert rc.tie_break_call(h_min, h_max, d_min, d_max) == expected


def _record(call_index, biomarkers, with_min):
    out = []
    for rid, exp, hmax, dmax, hmin, dmin in biomarkers:
        b = {"reaction": rid, "expected": exp, "healthy": hmax, "disease": dmax,
             "status_healthy": "Optimal", "status_disease": "Optimal"}
        if with_min:
            b = {"reaction": rid, "expected": exp, "healthy": None, "disease": None,
                 "status_healthy": "not_run", "status_disease": "not_run",
                 "healthy_min": hmin, "disease_min": dmin, "status_healthy_min": "Optimal", "status_disease_min": "Optimal"}
        out.append(b)
    return {"iem": f"I{call_index}", "call_index": call_index, "biomarkers": out}


def test_analysis_counts_transitions_and_excludes_incomplete_pairs():
    data = [
        ("EX_a[u]", "Increased", 10, 10, 0, 5),    # capped, max wrong, R1 and R2 right
        ("DM_b[bc]", "Decreased", 10, 4, 0, 0),    # max right, R1 right
        ("EX_c[u]", "Increased", 4, 8, 3, 1),      # max right, R1 conflicting (right -> wrong), R2 right
        ("EX_d[u]", "Increased", 10, 10, 0, 0),    # capped, all Unchanged
        ("EX_e[u]", "Unchanged", 1, 1, 0, 0),      # not directional: excluded
    ]
    max_run = [_record(1, data, with_min=False)]
    min_run = [_record(1, data, with_min=True)]
    rows = rc.pair_records(max_run, min_run)
    result = rc.analyse(rows)
    assert result["n_scored_four_optima"] == 4
    assert result["correct"] == {"max": 2, "R1": 2, "R2": 3}
    assert result["transitions_vs_max"]["R1"] == {"wrong_to_right": 1, "right_to_wrong": 1, "wrong_to_other_wrong": 0}
    assert result["transitions_vs_max"]["R2"] == {"wrong_to_right": 1, "right_to_wrong": 0, "wrong_to_other_wrong": 0}
    assert result["n_conflicting_R1"] == 1
    assert result["equal_positive_maxima"] == {"n": 2, "R1_correct": 1, "R2_correct": 1, "min_differs": 1}
    assert result["reading"] == {"R1": "no net change", "R2": "improves"}


def test_mismatched_reactions_are_refused():
    max_run = [_record(1, [("EX_a[u]", "Increased", 1, 2, 0, 0)], with_min=False)]
    min_run = [_record(1, [("EX_b[u]", "Increased", 1, 2, 0, 0)], with_min=True)]
    with pytest.raises(ValueError):
        rc.pair_records(max_run, min_run)
