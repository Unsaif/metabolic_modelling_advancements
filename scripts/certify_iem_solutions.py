"""Certify Python IEM optima independently of the solver, for selected IEMs.

Re-runs the frozen v0.3 protocol (toolbox setup, Toolbox bile-duct list, HiGHS IPM) for the named IEMs and checks
every returned optimal vector against the LP that HiGHS actually held at that moment (all rows, including the
added IEM-sum row and demand columns; all bounds): row and bound violations are recomputed in double precision and
the objective is recomputed from the vector. Used to check the biomarker optima for which the MATLAB/Gurobi
reference run returned no value (NaN): a certified feasible vector with the reported objective shows that the
problem is feasible and the Python value attainable.

Usage: python scripts/certify_iem_solutions.py Harvey_1_03d --iems ASNSD BTD EF GMT PHOX1 SUCLA --out f.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from gembench import wbm as W  # noqa: E402
from gembench import wbm_constraints as C  # noqa: E402
from gembench import wbm_iem as I  # noqa: E402


def certify_current(hw: I.HighsWBM, x: np.ndarray, reported: float) -> dict:
    lp = hw.h.getLp()
    n, m = lp.num_col_, lp.num_row_
    A = sp.csc_matrix((np.asarray(lp.a_matrix_.value_), np.asarray(lp.a_matrix_.index_), np.asarray(lp.a_matrix_.start_)),
                      shape=(m, n))
    x = np.asarray(x, dtype=float)[:n]
    r = A @ x
    rl, ru = np.asarray(lp.row_lower_), np.asarray(lp.row_upper_)
    cl, cu = np.asarray(lp.col_lower_), np.asarray(lp.col_upper_)
    row_viol = np.maximum(np.maximum(rl - r, r - ru), 0)
    col_viol = np.maximum(np.maximum(cl - x, x - cu), 0)
    obj = float(np.asarray(lp.col_cost_) @ x)
    return {"max_row_violation": float(row_viol.max()), "n_rows_over_1e-6": int((row_viol > 1e-6).sum()),
            "max_bound_violation": float(col_viol.max()), "n_bounds_over_1e-6": int((col_viol > 1e-6).sum()),
            "objective_recomputed": obj, "objective_reported": reported,
            "objective_abs_diff": abs(obj - reported) if np.isfinite(reported) else None}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model")
    ap.add_argument("--iems", nargs="+", required=True)
    ap.add_argument("--protocol", default=os.path.join(ROOT, "data", "iem", "iem_protocol_v0.2.json"))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    m = W.load_wbm(os.path.join(ROOT, "external", "COBRA.models", "mat", f"{args.model}.mat"))
    lb, ub, _ = C.runiem_model_setup(m.rxns, m.mets, m.S, m.lb, m.ub, sex=m.meta.get("sex") or "male")
    m.lb, m.ub = lb, ub
    hw = I.HighsWBM(m)
    I.apply_runiem_global_constraints(hw, bile_duct="toolbox")
    records = []
    original = hw.solve

    def hooked():
        status, obj, x, dt = original()
        cert = certify_current(hw, x, obj) if status == "Optimal" and np.all(np.isfinite(x)) else None
        obj_cols = getattr(hw, "_obj_cols", [])
        names = {v: k for k, v in hw.rxn_pos.items()}
        records.append({"status": status, "objective": obj, "seconds": dt,
                        "objective_reactions": [names.get(j, str(j)) for j in obj_cols][:3], "certificate": cert})
        return status, obj, x, dt

    hw.solve = hooked
    protocol = {p["iem"]: p for p in json.load(open(args.protocol))}
    out = {"model": args.model, "iems": {}}
    for iem in args.iems:
        p = protocol[iem]
        start = len(records)
        t = time.time()
        with hw.temporary_state():
            for tw in p["bound_tweaks"]:
                idx = I.match_reactions(m.rxns, [tw["pattern"]])
                hw.set_bounds(idx, **{tw["bound"]: [tw["value"]] * len(idx)})
            r = I.run_iem(hw, iem, p["include_patterns"], p["exclude_patterns"], [tuple(b) for b in p["biomarkers"]],
                          verbose=False, demand_metabolites=p.get("demand_metabolites"))
        solves = records[start:]
        # solve order in run_iem: IEM-flux maximum, disease IEM flux, whole-body objective, then healthy/disease per biomarker
        biomarker_solves = solves[3:]
        rows = []
        for k, b in enumerate(r.biomarkers):
            healthy = biomarker_solves[2 * k] if 2 * k < len(biomarker_solves) else None
            disease = biomarker_solves[2 * k + 1] if 2 * k + 1 < len(biomarker_solves) else None
            rows.append({"biomarker": b.reaction, "healthy": b.healthy, "disease": b.disease, "call": b.predicted,
                         "healthy_solve": healthy, "disease_solve": disease})
        out["iems"][iem] = {"status": r.status, "seconds": round(time.time() - t, 1), "n_solves": len(solves),
                            "biomarkers": rows, "setup_solves": solves[:3]}
        worst = max((s["certificate"]["max_row_violation"] for s in solves if s["certificate"]), default=None)
        print(f"{iem}: {r.status}, {len(solves)} solves, worst row violation {worst}", flush=True)
        with open(args.out, "w") as fh:
            json.dump(out, fh, indent=1, default=float)


if __name__ == "__main__":
    main()
