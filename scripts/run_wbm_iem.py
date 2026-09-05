"""Run the IEM biomarker protocol (checkIEM_WBM port) on a whole-body model with HiGHS warm starts.

Usage: python scripts/run_wbm_iem.py Harvey_1_03d [--limit N] [--iems HIS AGAT ...]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gembench import wbm as W  # noqa: E402
from gembench import wbm_iem as I  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "wbm_iem")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("model")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--iems", nargs="*", default=None)
    ap.add_argument("--min-flux-healthy", type=float, default=1.0)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    protocol = json.load(open(os.path.join(ROOT, "results", "iem_ground_truth", "iem_protocol_v0.json")))
    if args.iems:
        protocol = [p for p in protocol if p["iem"] in set(args.iems)]
    if args.limit:
        protocol = protocol[: args.limit]

    t0 = time.time()
    m = W.load_wbm(os.path.join(ROOT, "external", "COBRA.models", "mat", f"{args.model}.mat"))
    hw = I.HighsWBM(m)
    g = I.apply_runiem_global_constraints(hw)
    print(f"== {args.model}: {m.n_rxns} rxns; global constraints {g}; setup {time.time()-t0:.0f}s", flush=True)

    results = []
    out_json = os.path.join(OUT, f"{args.model}_iem_results.json")
    if os.path.exists(out_json):          # resume: keep IEMs already computed
        results = json.load(open(out_json))
        done = {(r["iem"], r["call_index"]) for r in results}
        protocol = [p for p in protocol if (p["iem"], p["call_index"]) not in done]
        print(f"resuming: {len(results)} IEMs already done, {len(protocol)} to go", flush=True)
    for p in protocol:
        # per-block bound tweaks (applied before the IEM and reverted after)
        tweaked = []
        for tw in p["bound_tweaks"]:
            idx = I.match_reactions(m.rxns, [tw["pattern"]])
            saved = (hw.lb[idx].copy(), hw.ub[idx].copy())
            if tw["bound"] == "lb":
                hw.set_bounds(idx, lb=[tw["value"]] * len(idx))
            else:
                hw.set_bounds(idx, ub=[tw["value"]] * len(idx))
            tweaked.append((idx, saved))
        r = I.run_iem(hw, p["iem"], p["include_patterns"], p["exclude_patterns"], [tuple(b) for b in p["biomarkers"]],
                      min_flux_healthy=args.min_flux_healthy)
        for idx, (lb, ub) in tweaked:
            hw.set_bounds(idx, lb=lb, ub=ub)
        scored = [b for b in r.biomarkers if b.correct is not None]
        n_ok = sum(1 for b in scored if b.correct)
        print(f"  {p['iem']:8s} rxns={len(r.iem_reactions):3d} vmax={r.vmax_healthy:10.4g} disease_vmax={r.vmax_disease:8.3g} "
              f"wb_feasible={r.wb_objective_disease_feasible} biomarkers={n_ok}/{len(scored)} correct  solves={r.n_solves} {r.time_s:.0f}s {r.notes}", flush=True)
        results.append({"iem": p["iem"], "call_index": p["call_index"], "iem_reactions": r.iem_reactions, "vmax_healthy": r.vmax_healthy,
                        "vmax_disease": r.vmax_disease, "wb_objective_disease_feasible": r.wb_objective_disease_feasible,
                        "n_solves": r.n_solves, "time_s": r.time_s, "notes": r.notes,
                        "biomarkers": [b.__dict__ for b in r.biomarkers]})
        with open(out_json, "w") as fh:
            json.dump(results, fh, indent=2, default=str)
    scored = [b for r in results for b in r["biomarkers"] if b["correct"] is not None]
    n_ok = sum(1 for b in scored if b["correct"])
    up_ok = sum(1 for b in scored if b["correct"] and b["expected"] == "Increased")
    up_all = sum(1 for b in scored if b["expected"] == "Increased")
    dn_ok = sum(1 for b in scored if b["correct"] and b["expected"] == "Decreased")
    dn_all = sum(1 for b in scored if b["expected"] == "Decreased")
    unch = sum(1 for b in scored if b["predicted"] == "Unchanged")
    summary = {"model": args.model, "n_iems_run": len(results), "n_biomarkers_scored": len(scored), "n_correct": n_ok,
               "accuracy": n_ok / len(scored) if scored else float("nan"), "increased_correct": f"{up_ok}/{up_all}",
               "decreased_correct": f"{dn_ok}/{dn_all}", "predicted_unchanged": unch,
               "total_solves": hw.n_solves, "total_solve_time_s": round(hw.solve_time, 1), "wall_s": round(time.time() - t0, 1),
               "caveats": ["shipped bounds only: physiologicalConstraintsHMDBbased, standard physiological parameters and the EU average diet "
                           "(applied in runIEM_HH.m before the IEMs) are NOT yet ported, so this is not the published protocol",
                           "global reaction constraints from the top of runIEM_HH.m and per-block lb/ub tweaks are applied",
                           "solver: HiGHS (IPM+crossover cold start, dual simplex warm starts), feasibility 1e-7"]}
    with open(os.path.join(OUT, f"{args.model}_iem_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
