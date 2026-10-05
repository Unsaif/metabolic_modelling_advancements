"""Tests for the MATLAB reference comparison (scripts/compare_matlab_reference.py)."""
import importlib.util
import math
import os

import numpy as np
import scipy.io as sio

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("cmr", os.path.join(ROOT, "scripts", "compare_matlab_reference.py"))
CMR = importlib.util.module_from_spec(spec)
spec.loader.exec_module(CMR)


def test_matlab_call_uses_runiem_rule_and_treats_missing_values_as_unchanged():
    assert CMR.matlab_call(0.0, 5.0) == "Increased"
    assert CMR.matlab_call(5.0, 0.0) == "Decreased"
    assert CMR.matlab_call(1.0, 1.0 + 5e-7) == "Unchanged"
    assert CMR.matlab_call(float("nan"), 1.0) == "Unchanged"
    assert math.isnan(CMR.matlab_value(np.array(["NA"])))


def test_iemsol_parsing_pairs_healthy_and_disease_rows(tmp_path):
    rows = [["IEM Rxns All obj - Healthy", "1", "", ""], ["IEM Rxns All obj - Disease", "0", "", ""],
            ["WB obj - Healthy", "ND", "", ""], ["WB obj - Disease", "1", "1", ""],
            ["Healthy:EX_a[u]", "0", "Disease - Reported:Increased (urine)", ""],
            ["Disease:EX_a[u]", "12.5", "Disease - Reported:Increased (urine)", ""],
            ["Healthy:DM_b[bc]", "3", "Disease - Reported:Decreased (blood)", ""],
            ["Disease:DM_b[bc]", "3", "Disease - Reported:Decreased (blood)", ""]]
    cell = np.empty((len(rows), 4), dtype=object)
    for i, row in enumerate(rows):
        for j, v in enumerate(row):
            cell[i, j] = v
    sio.savemat(tmp_path / "r.mat", {"IEMSol_TOY": cell})
    parsed = CMR.parse_iemsol(str(tmp_path / "r.mat"))
    assert [(r["biomarker"], r["call"]) for r in parsed["TOY"]] == [("EX_a[u]", "Increased"), ("DM_b[bc]", "Unchanged")]
    python = [{"iem": "TOY", "biomarkers": [
        {"reaction": "EX_a[u]", "expected": "Increased", "healthy": 0.0, "disease": 12.5, "predicted": "Increased"},
        {"reaction": "DM_b[bc]", "expected": "Decreased", "healthy": 3.0, "disease": 2.0, "predicted": "Decreased"}]}]
    report = CMR.compare_results(parsed, python)
    assert (report["n_compared"], report["n_same_call"], report["n_matlab_correct"]) == (2, 1, 1)
    assert report["differences"][0]["biomarker"] == "DM_b[bc]"


def test_bound_comparison_counts_identical_and_different(tmp_path):
    rx = ["r1", "r2", "r3"]
    out = CMR.compare_bounds(rx, np.array([0.0, -1.0, -np.inf]), np.array([1.0, 1.0, np.inf]),
                             rx[::-1], np.array([-np.inf, -1.0, 0.0]), np.array([np.inf, 2.0, 1.0]))
    assert out["lb"]["n_different"] == 0
    assert out["ub"]["n_different"] == 1 and out["ub"]["examples"][0]["rxn"] == "r2"
