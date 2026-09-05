"""Solve Harvey/Harvetta with open solvers, certify solutions, and record timings.

Usage: python scripts/run_wbm_solvers.py Harvey_1_03d Harvetta_1_03d --highs simplex ipm --glpk dual
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gembench import wbm as W  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "wbm_solvers")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("models", nargs="+")
    ap.add_argument("--highs", nargs="*", default=["simplex", "ipm"])
    ap.add_argument("--glpk", nargs="*", default=[])
    ap.add_argument("--objectives", nargs="*", default=["feasibility", "EX_phe_L[u]", "EX_glc_D[u]"])
    ap.add_argument("--glpk-time-limit-min", type=float, default=30)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for name in args.models:
        path = os.path.join(ROOT, "external", "COBRA.models", "mat", f"{name}.mat")
        m = W.load_wbm(path)
        print(f"== {name}: {m.S.shape[1]} rxns, {m.S.shape[0]} mets, {m.C.shape[0]} coupling rows; load {m.meta['load_s']}s", flush=True)
        for obj in args.objectives:
            if obj == "feasibility":
                c, sense = np.zeros(m.n_rxns), "min"
            else:
                if obj not in set(m.rxns):
                    print(f"   skip {obj}: not in model"); continue
                c, sense = W.unit_objective(m, obj), "max"
            for meth in args.highs:
                r = W.solve_highs(m, objective=c, sense=sense, method=meth)
                row = {"model": name, "objective": obj, "solver": r.solver, "method": meth, "status": r.status, "value": r.objective,
                       "time_s": round(r.time_s, 1), **{f"cert_{k}": v for k, v in r.certificate.items()}, **{f"info_{k}": v for k, v in r.info.items()}}
                rows.append(row)
                print(f"   HiGHS {meth:8s} {obj:14s} {r.status:10s} value={r.objective:.6g} {r.time_s:6.1f}s cert(maxS={r.certificate.get('max_abs_S_residual', float('nan')):.1e}, "
                      f"maxC={r.certificate.get('max_coupling_violation', float('nan')):.1e})", flush=True)
            for meth in args.glpk:
                r = W.solve_glpk(m, objective=c, sense=sense, method=meth, time_limit_ms=int(args.glpk_time_limit_min * 60_000))
                row = {"model": name, "objective": obj, "solver": r.solver, "method": meth, "status": r.status, "value": r.objective,
                       "time_s": round(r.time_s, 1), **{f"cert_{k}": v for k, v in r.certificate.items()}, **{f"info_{k}": v for k, v in r.info.items()}}
                rows.append(row)
                print(f"   GLPK  {meth:8s} {obj:14s} {r.status:10s} value={r.objective:.6g} {r.time_s:6.1f}s cert(maxS={r.certificate.get('max_abs_S_residual', float('nan')):.1e})", flush=True)
            with open(os.path.join(OUT, "solver_runs.json"), "w") as fh:
                json.dump(rows, fh, indent=2, default=str)


if __name__ == "__main__":
    main()
