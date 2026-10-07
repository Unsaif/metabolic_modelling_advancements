"""Re-solve selected LPs of a finished cross-disease matrix from scratch, optionally with another solver.

Usage:
  python scripts/recheck_cross_matrix.py results/wbm_iem/Harvey_1_03d_iem_cross_ranking_v1.json \
         --backend highs --extra-readouts data/iem/iem_ranking_extra_readouts_v0.2.txt

The LPs are rebuilt exactly as scripts/run_wbm_iem_cross.py builds them: the same model set-up, global
constraints and sink columns, each IEM's tweaks, sinks and healthy pin, taken from the matrix record (prepare(...,
stored=record)), so no set-up solve is repeated. Selected are the solves the run marked for a recheck (the hash
selection of --recheck-every), or --every N to select by the same hash anew. Each selected LP is solved from scratch
by interior point with crossover, and the value is compared with the matrix:
  - the absolute and relative difference;
  - whether the call for that IEM and readout changes when the recheck value replaces the matrix value.
Writes results/wbm_iem/ranking/<matrix stem>_recheck_<backend>.json, resumable per readout.

With another solver, a healthy state can be infeasible only because the matrix's pin (that solver's v_max truncated
to six decimals) sits just above this solver's maximum. --retry-infeasible-pins FILE re-solves such rechecks once with
the pin from FILE's v_max for that IEM, computed by this solver; the first attempt is kept under "first_attempt".
Calls are counted as changed only for rechecks that end optimal; the others are listed in the summary.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from gembench import wbm as W  # noqa: E402
from gembench import wbm_iem as I  # noqa: E402
from gembench import wbm_constraints as WC  # noqa: E402

spec = importlib.util.spec_from_file_location("run_wbm_iem_cross", os.path.join(ROOT, "scripts", "run_wbm_iem_cross.py"))
X = importlib.util.module_from_spec(spec)
spec.loader.exec_module(X)


def truncated_pin(v):
    """The protocol's healthy pin from a maximum: truncated to six decimals, towards zero (as in prepare())."""
    return math.floor(v * 1e6) / 1e6 if v > 0 else math.ceil(v * 1e6) / 1e6


def solve_one(hw, mover, s, state, rcol, alt_lo=None):
    """Solve one LP of the matrix from scratch by interior point with crossover.

    If a healthy state is infeasible and alt_lo is given, it is solved once more with alt_lo as the pin. The matrix's
    pin comes from another solver's v_max and, within tolerance, can sit just above this solver's maximum.
    Returns (status, value, seconds, solver info, first attempt or None); s.lo is left as it was."""
    target, row_bounds = X.target_bounds(s, state, rcol, mover.base_lb, mover.base_ub)
    mover.move(target, s.row, row_bounds)
    hw.fresh()
    hw.method = "ipm"
    status, f, _, dt = hw.solve()
    info = dict(getattr(hw, "last_info", {}) or {})
    first_attempt = None
    if status == "Infeasible" and state == "healthy" and alt_lo is not None:
        first_attempt = {"status": status, "pin": s.lo, "time_s": round(dt, 2)}
        lo_matrix, s.lo = s.lo, alt_lo
        try:
            target, row_bounds = X.target_bounds(s, state, rcol, mover.base_lb, mover.base_ub)
            mover.move(target, s.row, row_bounds)
            hw.fresh()
            status, f, _, dt = hw.solve()
            info = dict(getattr(hw, "last_info", {}) or {})
            first_attempt["retry_pin"] = s.lo
        finally:
            s.lo = lo_matrix
    return status, f, dt, info, first_attempt


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("matrix")
    ap.add_argument("--backend", choices=["highs", "gurobi"], default="highs")
    ap.add_argument("--extra-readouts", required=True)
    ap.add_argument("--protocol", default=os.path.join(ROOT, "data", "iem", "iem_protocol_v0.2.json"))
    ap.add_argument("--every", type=int, default=0, help="select by the recheck hash with this N instead of the run's marks")
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--out")
    ap.add_argument("--retry-infeasible-pins", help="results file whose v_max per IEM (computed by this solver) sets the "
                    "healthy pin for a second attempt at rechecks this solver found infeasible in the healthy state")
    args = ap.parse_args()

    with open(args.matrix) as fh:
        recs = json.load(fh)
    prov = recs[0]["provenance"]
    with open(args.protocol) as fh:
        protocol = json.load(fh)
    if X.sha256_file(args.protocol) != prov["protocol_sha256"]:
        sys.exit("protocol differs from the matrix's")
    if X.sha256_file(args.extra_readouts) != prov.get("extra_readouts_sha256"):
        sys.exit("extra readouts differ from the matrix's")
    panel = list(dict.fromkeys(rid for p in protocol for rid, _ in p["biomarkers"]))
    with open(args.extra_readouts) as fh:
        panel = list(dict.fromkeys(panel + [line.strip() for line in fh if line.strip()]))
    by_iem = {r["iem"]: r for r in recs}

    # Select the LPs.
    selected = []
    for r in recs:
        for e in r["readouts"]:
            for st in ("healthy", "disease"):
                info = e[f"solve_{st}"] or {}
                if info.get("method") in (None, "ipm"):
                    continue
                pick = X.rechecked(r["iem"], st, e["reaction"], args.every) if args.every else "recheck" in info
                if pick:
                    selected.append((e["reaction"], r["iem"], st))
    print(f"{len(selected)} LPs selected", flush=True)

    # Build the LP exactly as the runner does.
    m = W.load_wbm(os.path.join(ROOT, "external", "COBRA.models", "mat", "Harvey_1_03d.mat"))
    if X.sha256_file(m.meta["file"]) != prov["model_sha256"]:
        sys.exit("model file differs from the matrix's")
    inputs = WC.load_inputs(WC.DEFAULT_INPUTS)
    lb, ub, _ = WC.runiem_model_setup(m.rxns, m.mets, m.S, m.lb, m.ub, sex=m.meta.get("sex") or "male", inputs=inputs)
    m.lb, m.ub = lb, ub
    hw = X.make_backend(args.backend, m, args.threads)
    I.apply_runiem_global_constraints(hw, bile_duct="toolbox")
    bounds = hashlib.sha256(np.ascontiguousarray(hw.lb).tobytes() + np.ascontiguousarray(hw.ub).tobytes()).hexdigest()
    if bounds != prov["lp_bounds_sha256_before_sinks"]:
        sys.exit("LP bounds differ from the matrix's")
    sink_mets = list(dict.fromkeys([rid[3:] for rid in panel if rid.startswith("DM_")] +
                                   [met for p in protocol for met in (p.get("demand_metabolites") or [])]))
    for met in sink_mets:
        if met in hw.met_pos:
            hw.add_demand(met, lb=0.0, ub=0.0)
    setups = {p["iem"]: X.prepare(hw, p, prov["context"], log=False, stored=by_iem[p["iem"]]) for p in protocol}
    mover = X.Mover(hw)

    stem = os.path.splitext(os.path.basename(args.matrix))[0].replace("_iem_cross_ranking", "_ranking")
    out = args.out or os.path.join(ROOT, "results", "wbm_iem", "ranking", f"{stem}_recheck_{args.backend}.json")
    done = {}
    if os.path.exists(out):
        with open(out) as fh:
            old = json.load(fh)
        done = {(x["reaction"], x["iem"], x["state"]): x for x in old["rechecks"]}
    solver = X.backend_version(args.backend)
    alt_pin = {}
    if args.retry_infeasible_pins:
        with open(args.retry_infeasible_pins) as fh:
            for r in json.load(fh):
                v = r.get("vmax_healthy")
                if v is not None and math.isfinite(v):
                    alt_pin[r["iem"]] = truncated_pin(v)
        for key, x in list(done.items()):
            if x["status"] == "Infeasible" and x["state"] == "healthy" and x["iem"] in alt_pin and "first_attempt" not in x:
                del done[key]     # retried below with the other pin
    by_readout = {}
    for rid, iem, st in selected:
        by_readout.setdefault(rid, []).append((iem, st))
    t0 = time.time()
    for k, (rid, items) in enumerate(by_readout.items()):
        todo = [(iem, st) for iem, st in items if (rid, iem, st) not in done]
        if not todo:
            continue
        rcol = hw.rxn_pos[rid]
        hw.set_objective({rcol: 1.0}, "max")
        for iem, st in todo:
            status, f, dt, info, first_attempt = solve_one(hw, mover, setups[iem], st, rcol, alt_pin.get(iem))
            e = next(x for x in by_iem[iem]["readouts"] if x["reaction"] == rid)
            ok = status == "Optimal" and math.isfinite(f)
            fz = (0.0 if abs(f) <= X.TOL else f) if ok else float("nan")
            stored = e[st]
            other = e["disease" if st == "healthy" else "healthy"]
            h, d = (fz, other) if st == "healthy" else (other, fz)
            new_call = X.call(h if h is not None else float("nan"), d if d is not None else float("nan"))
            diff = abs(fz - stored) if ok and stored is not None else None
            scale = max(abs(fz), abs(stored)) if ok and stored is not None else 0.0
            done[(rid, iem, st)] = {"reaction": rid, "iem": iem, "state": st, "status": status, "value": fz,
                                    "matrix_value": stored, "abs_diff": diff,
                                    "rel_diff": diff / scale if diff is not None and scale > 1e-3 else None,
                                    "matrix_call": e["predicted"], "call_with_recheck": new_call,
                                    "call_changes": ok and new_call != e["predicted"], "time_s": round(dt, 2),
                                    "ipm_iterations": info.get("ipm_iterations")}
            if first_attempt:
                done[(rid, iem, st)]["first_attempt"] = first_attempt
        rows = list(done.values())
        diffs = [x["abs_diff"] for x in rows if x["abs_diff"] is not None]
        summary = {"n": len(rows), "n_selected": len(selected), "non_optimal": sum(x["status"] != "Optimal" for x in rows),
                   "non_optimal_list": [(x["iem"], x["reaction"], x["state"], x["status"]) for x in rows if x["status"] != "Optimal"],
                   "retried_with_own_pin": [(x["iem"], x["reaction"], x["state"], x["status"]) for x in rows if "first_attempt" in x],
                   "calls_changed": sum(bool(x["call_changes"]) for x in rows),
                   "max_abs_diff": max(diffs) if diffs else None,
                   "n_abs_diff_above_1e-6": sum(x > 1e-6 for x in diffs),
                   "max_rel_diff_values_above_1e-3": max([x["rel_diff"] for x in rows if x["rel_diff"] is not None] or [0])}
        X.write_json(out, {"matrix": os.path.relpath(os.path.abspath(args.matrix), ROOT), "matrix_sha256": X.sha256_file(args.matrix),
                           "solver": solver, "method": "interior point with crossover, from scratch",
                           "selection": f"hash every {args.every}" if args.every else "the run's recheck marks",
                           "summary": summary, "rechecks": rows})
        print(f"  readout {k + 1}/{len(by_readout)} {rid}: {len(rows)}/{len(selected)} done, calls changed "
              f"{summary['calls_changed']}, max |diff| {summary['max_abs_diff']}, {time.time() - t0:.0f}s", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
