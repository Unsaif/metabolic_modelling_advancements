"""Run the whole-body IEM protocol of runIEM_HH.m with HiGHS.

Usage: python scripts/run_wbm_iem.py Harvey_1_03d [--model-setup toolbox|shipped] [--bile-duct toolbox|v0.2_all]
                                     [--limit N] [--iems HIS AGAT ...] [--out-suffix _v0.3]

--model-setup toolbox  re-applies physiologicalConstraintsHMDBbased and the EU average diet as runIEM_HH.m
                       does at the pinned Toolbox commit (gembench.wbm_constraints);
--model-setup shipped  keeps the bounds stored in the model file (Harvey 1.03d ships with an earlier
                       version of the same constraints already applied).
--bile-duct v0.2_all   reproduces the v0.2 deviation (ub = 100 on all 261 bile-duct exits instead of 28).
--senses min           computes the minimum biomarker flux in each state instead of the maximum (plan v0.4
                       flux-range analysis; scripts/iem_range_calls.py joins it with a maximum run);
                       --senses max min computes both. The default (max) is the runIEM_HH protocol.
Results are fingerprinted by model, protocol, setup and the resulting LP bounds; a results file is
only resumed under the same fingerprint.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gembench import wbm as W  # noqa: E402
from gembench import wbm_iem as I  # noqa: E402
from gembench import wbm_constraints as WC  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "wbm_iem")
ENGINE_VERSION = "iem-v0.3"
CAVEATS = {
    "toolbox": "Model bounds: the shipped model with physiologicalConstraintsHMDBbased and the EU average diet re-applied "
               "as runIEM_HH.m does at the pinned Toolbox commit (gembench.wbm_constraints).",
    "shipped": "Model bounds as shipped; Harvey 1.03d already carries an earlier version of the physiological and diet "
               "constraints (GFR 129.75 ml/min, CSF export from 0.35 ml/min, older AGORA and diet lists).",
}
TERMINAL_STATUSES = {"complete", "inactive", "no_reactions", "disease_infeasible", "wb_infeasible", "missing_wb_objective"}


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def scored_biomarkers(results):
    """Never trust a stored correctness flag if either optimum is unavailable."""
    return [b for r in results for b in r["biomarkers"]
            if b["correct"] is not None and b["status_healthy"] == b["status_disease"] == "Optimal"
            and finite(b["healthy"]) and finite(b["disease"])]


def load_resume(path, fingerprint):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        results = json.load(fh)
    if any(r.get("run_fingerprint") != fingerprint for r in results):
        raise ValueError("Existing results use a different or unrecorded model/protocol/settings. Choose a new --out-suffix.")
    return results


def json_safe(value):
    """Write unavailable numbers as JSON null, not nonstandard NaN/Infinity."""
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


def write_json(path, value):
    temporary = path + ".tmp"
    with open(temporary, "w") as fh:
        json.dump(json_safe(value), fh, indent=2, allow_nan=False)
    os.replace(temporary, path)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("model")
    ap.add_argument("--model-file", help="MATLAB model file; defaults to external/COBRA.models/mat/<model>.mat")
    ap.add_argument("--protocol", default=os.path.join(ROOT, "data", "iem", "iem_protocol_v0.2.json"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--iems", nargs="*", default=None)
    ap.add_argument("--min-flux-healthy", type=float, default=1.0)
    ap.add_argument("--out-suffix", default="_v0.3", help="Use distinct suffixes for disjoint shards; legacy results cannot be resumed")
    ap.add_argument("--model-setup", choices=["toolbox", "shipped"], default="toolbox")
    ap.add_argument("--bile-duct", choices=["toolbox", "v0.2_all"], default="toolbox")
    ap.add_argument("--constraint-inputs", default=WC.DEFAULT_INPUTS)
    ap.add_argument("--senses", nargs="+", choices=["max", "min"], default=["max"],
                    help="Biomarker optima per state: max (runIEM_HH protocol), min (flux-range analysis), or both")
    args = ap.parse_args()
    if len(set(args.senses)) != len(args.senses):
        ap.error("--senses must not repeat a sense")
    os.makedirs(OUT, exist_ok=True)
    with open(args.protocol) as fh:
        full_protocol = json.load(fh)
    protocol = full_protocol
    if args.iems:
        unknown = set(args.iems) - {p["iem"] for p in protocol}
        if unknown:
            ap.error(f"IEMs not in this protocol: {sorted(unknown)}")
        protocol = [p for p in protocol if p["iem"] in set(args.iems)]
    if args.limit:
        protocol = protocol[:args.limit]
    if not math.isfinite(args.min_flux_healthy) or not 0 <= args.min_flux_healthy <= 1:
        ap.error("--min-flux-healthy must be a finite fraction between 0 and 1")

    t0 = time.time()
    model_file = args.model_file or os.path.join(ROOT, "external", "COBRA.models", "mat", f"{args.model}.mat")
    import highspy
    m = W.load_wbm(model_file)
    setup = {"model_setup": args.model_setup}
    if args.model_setup == "toolbox":
        inputs = WC.load_inputs(args.constraint_inputs)
        lb, ub, report = WC.runiem_model_setup(m.rxns, m.mets, m.S, m.lb, m.ub, sex=m.meta.get("sex") or "male", inputs=inputs)
        setup.update(constraint_inputs_sha256=sha256_file(args.constraint_inputs),
                     toolbox_commit=inputs["provenance"]["toolbox_commit"],
                     n_lb_changed=int((lb != m.lb).sum()), n_ub_changed=int((ub != m.ub).sum()),
                     warnings=report["warnings"], parameters=report["parameters"])
        m.lb, m.ub = lb, ub
    hw = I.HighsWBM(m)
    g = I.apply_runiem_global_constraints(hw, bile_duct=args.bile_duct)
    bounds_sha256 = hashlib.sha256(np.ascontiguousarray(hw.lb).tobytes() + np.ascontiguousarray(hw.ub).tobytes()).hexdigest()
    provenance = {"engine_version": ENGINE_VERSION, "model_sha256": sha256_file(model_file),
                  "protocol_sha256": sha256_file(args.protocol), "min_flux_healthy": args.min_flux_healthy,
                  "solver": f"HiGHS {highspy.Highs().version()}", "method": "ipm", "feas_tol": 1e-7,
                  "opt_tol": 1e-7, "threads": 0, "time_limit_per_solve_s": 1800,
                  "model_setup": setup, "bile_duct": args.bile_duct, "global_constraints": g,
                  "lp_bounds_sha256": bounds_sha256}
    if args.senses != ["max"]:
        # Recorded only when it differs from the protocol default, so max-only runs keep their v0.3 fingerprints.
        provenance["senses"] = list(args.senses)
    fingerprint = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
    out_json = os.path.join(OUT, f"{args.model}_iem_results{args.out_suffix}.json")
    results = load_resume(out_json, fingerprint)
    # Retry interrupted or numerically unavailable IEMs, replacing their previous attempt.
    done = {(r["iem"], r["call_index"]) for r in results if r.get("status") in TERMINAL_STATUSES}
    protocol = [p for p in protocol if (p["iem"], p["call_index"]) not in done]
    print(f"== {args.model}: {m.n_rxns} rxns; setup {args.model_setup}; global constraints {g}; "
          f"bounds {bounds_sha256[:12]}; {len(protocol)} IEMs to run", flush=True)

    for p in protocol:
        # Reset every per-IEM tweak even if a solver raises an exception.
        with hw.temporary_state():
            for tw in p["bound_tweaks"]:
                idx = I.match_reactions(m.rxns, [tw["pattern"]])
                hw.set_bounds(idx, **{tw["bound"]: [tw["value"]] * len(idx)})
            r = I.run_iem(hw, p["iem"], p["include_patterns"], p["exclude_patterns"],
                          [tuple(b) for b in p["biomarkers"]], min_flux_healthy=args.min_flux_healthy,
                          demand_metabolites=p.get("demand_metabolites"), senses=tuple(args.senses))
        record = asdict(r)
        record.update(call_index=p["call_index"], run_fingerprint=fingerprint, provenance=provenance)
        scored = scored_biomarkers([record])
        print(f"  {p['iem']:8s} status={r.status} biomarkers={sum(b['correct'] for b in scored)}/{len(scored)} "
              f"correct; expected={len(p['biomarkers'])}; solves={r.n_solves} {r.time_s:.0f}s {r.notes}", flush=True)
        results = [old for old in results if (old["iem"], old["call_index"]) != (p["iem"], p["call_index"])]
        results.append(record)
        results.sort(key=lambda item: item["call_index"])
        write_json(out_json, results)

    scored = scored_biomarkers(results)
    n_ok = sum(b["correct"] for b in scored)
    by_key = {(p["iem"], p["call_index"]): p for p in full_protocol}
    n_expected = sum(len(by_key[(r["iem"], r["call_index"])]["biomarkers"]) for r in results)
    summary = {"model": args.model, "run_fingerprint": fingerprint, "provenance": provenance,
               "n_iems_attempted": len(results), "n_iems_with_complete_biomarker_solves": sum(r.get("status") == "complete" for r in results),
               "n_iems_in_protocol": len(full_protocol), "n_biomarkers_expected_in_attempted_iems": n_expected,
               "n_biomarkers_scored": len(scored), "n_correct": n_ok,
               "accuracy_among_scored": n_ok / len(scored) if scored else None,
               "scored_coverage_of_attempted_iems": len(scored) / n_expected if n_expected else None,
               "n_biomarkers_expected_in_full_protocol": sum(len(p["biomarkers"]) for p in full_protocol),
               "total_solves_in_recorded_attempts": sum(r["n_solves"] for r in results),
               "session_solves": hw.n_solves, "session_solve_time_s": round(hw.solve_time, 1),
               "session_wall_s": round(time.time() - t0, 1),
               "senses": list(args.senses),
               "n_biomarkers_with_optimal_minima": sum(
                   1 for r in results for b in r["biomarkers"]
                   if b.get("status_healthy_min") == b.get("status_disease_min") == "Optimal"
                   and finite(b.get("healthy_min")) and finite(b.get("disease_min"))),
               "caveats": [CAVEATS[args.model_setup],
                           "Accuracy is conditional on scored, optimal finite solve pairs; report coverage alongside accuracy.",
                           "Demand sinks are isolated per IEM. Solver is HiGHS IPM with crossover for every solve."]}
    write_json(os.path.join(OUT, f"{args.model}_iem_summary{args.out_suffix}.json"), summary)
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
