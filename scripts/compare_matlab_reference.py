"""Compare a MATLAB reference run of runIEM_HH.m with the Python port, bound by bound and biomarker by biomarker.

Inputs: the files written by tools/matlab/run_iem_reference_harvey_1_03d.m (directory --matlab-dir):
  matlab_setup_bounds_Harvey_1_03d.mat   bounds after the diet and physiological setup
  matlab_global_bounds_Harvey_1_03d.mat  bounds after the unified reaction constraints (optional)
  matlab_iem_results_Harvey_1_03d.mat    IEMSol_* cell arrays (optional)
and the Python results file of the same protocol (default: results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json).

MATLAB stores biomarker optima as num2str strings and runIEM_HH compares those strings; the MATLAB call here is
made the same way (healthy minus disease of the printed values, 1e-6 threshold). Python calls use full precision.

Usage: python scripts/compare_matlab_reference.py --matlab-dir <dir> [--python-results f.json] [--out report.json]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np
import scipy.io as sio

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from gembench import wbm as W  # noqa: E402
from gembench import wbm_constraints as C  # noqa: E402
from gembench import wbm_iem as I  # noqa: E402


def _cellstr(a) -> list:
    return [str(np.asarray(x).ravel()[0]) if np.asarray(x).size else "" for x in np.asarray(a, dtype=object).ravel()]


def load_bounds(path: str):
    d = sio.loadmat(path, squeeze_me=False)
    return _cellstr(d["rxns"]), np.asarray(d["lb"], dtype=float).ravel(), np.asarray(d["ub"], dtype=float).ravel()


def compare_bounds(rxns_m, lb_m, ub_m, rxns_p, lb_p, ub_p) -> dict:
    pos = {r: i for i, r in enumerate(rxns_p)}
    missing = [r for r in rxns_m if r not in pos]
    idx = np.array([pos[r] for r in rxns_m if r in pos])
    keep = np.array([r in pos for r in rxns_m])
    out = {"n_matlab": len(rxns_m), "n_python": len(rxns_p), "n_missing_in_python": len(missing),
           "missing_examples": missing[:10]}
    for name, a, b in (("lb", lb_m[keep], lb_p[idx]), ("ub", ub_m[keep], ub_p[idx])):
        same = (a == b) | (np.isinf(a) & np.isinf(b) & (np.sign(a) == np.sign(b)))
        diff = np.abs(a - b)
        rel = diff / np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)
        bad = np.where(~same)[0]
        out[name] = {"n_identical": int(same.sum()), "n_different": int((~same).sum()),
                     "n_relative_diff_over_1e-12": int((rel[~same] > 1e-12).sum()) if len(bad) else 0,
                     "max_abs_diff": float(diff[~same].max()) if len(bad) else 0.0,
                     "examples": [{"rxn": rxns_m[np.where(keep)[0][k]], "matlab": float(a[k]), "python": float(b[k])}
                                  for k in bad[:15]]}
    return out


def matlab_value(cell) -> float:
    s = str(np.asarray(cell).ravel()[0]) if np.asarray(cell).size else ""
    try:
        return float(s)
    except ValueError:
        return float("nan")      # 'NA', 'NaN', '' -> no value (str2num gives [] or NaN)


def matlab_call(h: float, d: float) -> str:
    if not (math.isfinite(h) and math.isfinite(d)):
        return "Unchanged"       # runIEM_HH: comparisons with NaN or [] are false -> 'unchanged'
    hd = h - d
    return "Increased" if hd < -1e-6 else "Decreased" if hd > 1e-6 else "Unchanged"


def parse_iemsol(path: str) -> dict:
    d = sio.loadmat(path, squeeze_me=False)
    out = {}
    for key, value in d.items():
        if not key.startswith("IEMSol_"):
            continue
        cells = np.asarray(value, dtype=object)
        rows = []
        for j in range(4, cells.shape[0] - 1, 2):     # MATLAB j = 5:2:end
            name = str(np.asarray(cells[j, 0]).ravel()[0]) if np.asarray(cells[j, 0]).size else ""
            if not name.startswith("Healthy:"):
                continue
            rxn = name[len("Healthy:"):]
            label = str(np.asarray(cells[j, 2]).ravel()[0]) if cells.shape[1] > 2 and np.asarray(cells[j, 2]).size else ""
            h, dval = matlab_value(cells[j, 1]), matlab_value(cells[j + 1, 1])
            rows.append({"biomarker": rxn, "label": label, "healthy": h, "disease": dval, "call": matlab_call(h, dval)})
        out[key[len("IEMSol_"):]] = rows
    return out


def compare_results(matlab: dict, python_results: list) -> dict:
    py = {(r["iem"], b["reaction"]): b for r in python_results for b in r["biomarkers"]}
    rows, agree, n_both = [], 0, 0
    for iem, items in sorted(matlab.items()):
        for m in items:
            p = py.get((iem, m["biomarker"]))
            entry = {"iem": iem, "biomarker": m["biomarker"], "matlab_healthy": m["healthy"], "matlab_disease": m["disease"],
                     "matlab_call": m["call"]}
            if p is None:
                entry["python_call"] = "absent"
            else:
                entry.update(python_healthy=p["healthy"], python_disease=p["disease"], python_call=p["predicted"])
                if p["predicted"] != "NA":
                    n_both += 1
                    agree += p["predicted"] == m["call"]
            rows.append(entry)
    expected = {(r["iem"], b["reaction"]): b["expected"] for r in python_results for b in r["biomarkers"]}
    m_correct = sum(1 for e in rows if expected.get((e["iem"], e["biomarker"])) == e["matlab_call"])
    return {"n_matlab_biomarkers": len(rows), "n_compared": n_both, "n_same_call": agree,
            "n_matlab_correct": m_correct,
            "differences": [e for e in rows if e.get("python_call") not in (None, "absent", "NA") and e["python_call"] != e["matlab_call"]],
            "rows": rows}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matlab-dir", required=True)
    ap.add_argument("--model", default="Harvey_1_03d")
    ap.add_argument("--python-results", default=os.path.join(ROOT, "results", "wbm_iem", "Harvey_1_03d_iem_results_v0.3.json"))
    ap.add_argument("--out")
    args = ap.parse_args()
    model_file = os.path.join(ROOT, "external", "COBRA.models", "mat", f"{args.model}.mat")
    m = W.load_wbm(model_file)
    lb_p, ub_p, _ = C.runiem_model_setup(m.rxns, m.mets, m.S, m.lb, m.ub, sex=m.meta.get("sex") or "male")
    report = {"matlab_dir": args.matlab_dir}
    f = os.path.join(args.matlab_dir, f"matlab_setup_bounds_{args.model}.mat")
    if os.path.exists(f):
        report["setup_bounds"] = compare_bounds(*load_bounds(f), list(m.rxns), lb_p, ub_p)
    f = os.path.join(args.matlab_dir, f"matlab_global_bounds_{args.model}.mat")
    if os.path.exists(f):
        m.lb, m.ub = lb_p, ub_p
        hw = I.HighsWBM(m)
        I.apply_runiem_global_constraints(hw, bile_duct="toolbox")
        report["global_bounds"] = compare_bounds(*load_bounds(f), list(m.rxns), hw.lb[:m.n_rxns], hw.ub[:m.n_rxns])
    f = os.path.join(args.matlab_dir, f"matlab_iem_results_{args.model}.mat")
    if os.path.exists(f) and os.path.exists(args.python_results):
        report["biomarkers"] = compare_results(parse_iemsol(f), json.load(open(args.python_results)))
    summary = {k: ({b: {kk: vv for kk, vv in v[b].items() if kk != "examples"} for b in ("lb", "ub")}
                   if k.endswith("bounds") else {kk: vv for kk, vv in v.items() if kk not in ("rows", "differences")})
               for k, v in report.items() if k != "matlab_dir"}
    print(json.dumps(summary, indent=1))
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(report, fh, indent=1, default=float)
            fh.write("\n")


if __name__ == "__main__":
    main()
