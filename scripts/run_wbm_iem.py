"""Run the corrected exploratory whole-body IEM protocol with HiGHS.

Usage: python scripts/run_wbm_iem.py Harvey_1_03d [--limit N] [--iems HIS AGAT ...]
The v0.2 output suffix prevents overwriting earlier exploratory results.
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

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gembench import wbm as W  # noqa: E402
from gembench import wbm_iem as I  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "wbm_iem")
ENGINE_VERSION = "iem-v0.2"
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
    ap.add_argument("--out-suffix", default="_v0.2", help="Use distinct suffixes for disjoint shards; legacy results cannot be resumed")
    args = ap.parse_args()
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
    provenance = {"engine_version": ENGINE_VERSION, "model_sha256": sha256_file(model_file),
                  "protocol_sha256": sha256_file(args.protocol), "min_flux_healthy": args.min_flux_healthy,
                  "solver": f"HiGHS {highspy.Highs().version()}", "method": "ipm", "feas_tol": 1e-7,
                  "opt_tol": 1e-7, "threads": 0, "time_limit_per_solve_s": 1800, "physiology_and_diet_ported": False}
    fingerprint = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
    out_json = os.path.join(OUT, f"{args.model}_iem_results{args.out_suffix}.json")
    results = load_resume(out_json, fingerprint)
    # Retry interrupted or numerically unavailable IEMs, replacing their previous attempt.
    done = {(r["iem"], r["call_index"]) for r in results if r.get("status") in TERMINAL_STATUSES}
    protocol = [p for p in protocol if (p["iem"], p["call_index"]) not in done]
    m = W.load_wbm(model_file)
    hw = I.HighsWBM(m)
    g = I.apply_runiem_global_constraints(hw)
    print(f"== {args.model}: {m.n_rxns} rxns; global constraints {g}; {len(protocol)} IEMs to run", flush=True)

    for p in protocol:
        # Reset every per-IEM tweak even if a solver raises an exception.
        with hw.temporary_state():
            for tw in p["bound_tweaks"]:
                idx = I.match_reactions(m.rxns, [tw["pattern"]])
                hw.set_bounds(idx, **{tw["bound"]: [tw["value"]] * len(idx)})
            r = I.run_iem(hw, p["iem"], p["include_patterns"], p["exclude_patterns"],
                          [tuple(b) for b in p["biomarkers"]], min_flux_healthy=args.min_flux_healthy,
                          demand_metabolites=p.get("demand_metabolites"))
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
               "caveats": ["Shipped bounds only: physiological constraints and the EU average diet are not ported; this is not a reproduction of the published results.",
                           "Accuracy is conditional on scored, optimal finite solve pairs; report coverage alongside accuracy.",
                           "Demand sinks are isolated per IEM. Solver is HiGHS IPM with crossover for every solve."]}
    write_json(os.path.join(OUT, f"{args.model}_iem_summary{args.out_suffix}.json"), summary)
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
