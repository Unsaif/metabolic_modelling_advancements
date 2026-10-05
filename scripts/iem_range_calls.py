"""Flux-range calls for whole-body IEM biomarkers (plan v0.4, frozen before any minimum was computed).

The runIEM_HH protocol (checkIEM_WBM) compares only the MAXIMUM biomarker flux of the healthy and the disease
state. When both states reach the same cap (for example a kidney filtration limit), the call is "Unchanged"
whatever happens below the cap. This script adds the MINIMUM flux of each state (a run of
scripts/run_wbm_iem.py with --senses min) and recomputes the calls with two rules fixed in advance:

R1, full range (in the spirit of Shlomi, Cabili and Ruppin 2009): with lo = sign(d_min - h_min) and
    hi = sign(d_max - h_max), each 0 when the absolute difference is at most 1e-6:
      Increased   if at least one of lo, hi is +1 and neither is -1;
      Decreased   if at least one is -1 and neither is +1;
      Unchanged   if both are 0;
      Conflicting if one is +1 and the other -1 (counted as wrong).
R2, tie-break: keep the protocol (maximum) call unless it is Unchanged; then call from the minima alone
    (Increased if d_min - h_min > 1e-6, Decreased if < -1e-6, otherwise Unchanged).

Scored set: biomarkers with an expected direction (Increased or Decreased) whose four optima (healthy and
disease maximum and minimum) are all Optimal and finite. Counts outside that set are reported.

Usage: python scripts/iem_range_calls.py Harvey=max.json,min.json [Harvetta=max.json,min.json ...] [--out f.json]
"""
from __future__ import annotations

import argparse
import json
import math

TOL = 1e-6


def _finite(v) -> bool:
    return isinstance(v, (int, float)) and math.isfinite(v)


def _sign(diff: float) -> int:
    return 1 if diff > TOL else -1 if diff < -TOL else 0


def max_call(h_max: float, d_max: float) -> str:
    s = _sign(d_max - h_max)
    return "Increased" if s > 0 else "Decreased" if s < 0 else "Unchanged"


def range_call(h_min: float, h_max: float, d_min: float, d_max: float) -> str:
    lo, hi = _sign(d_min - h_min), _sign(d_max - h_max)
    if (lo > 0 or hi > 0) and lo >= 0 and hi >= 0:
        return "Increased"
    if (lo < 0 or hi < 0) and lo <= 0 and hi <= 0:
        return "Decreased"
    if lo == 0 and hi == 0:
        return "Unchanged"
    return "Conflicting"


def tie_break_call(h_min: float, h_max: float, d_min: float, d_max: float) -> str:
    call = max_call(h_max, d_max)
    if call != "Unchanged":
        return call
    s = _sign(d_min - h_min)
    return "Increased" if s > 0 else "Decreased" if s < 0 else "Unchanged"


def biofluid(reaction: str) -> str:
    return "urine" if reaction.endswith("[u]") else "csf" if reaction.endswith("[csf]") else \
        "blood" if reaction.endswith("[bc]") else "other"


def pair_records(max_results: list, min_results: list) -> list:
    """Join the two runs biomarker by biomarker (same IEM call, same position, same reaction)."""
    mins = {r["call_index"]: r for r in min_results}
    rows = []
    for r in sorted(max_results, key=lambda x: x["call_index"]):
        m = mins.get(r["call_index"])
        for k, b in enumerate(r["biomarkers"]):
            bm = m["biomarkers"][k] if m is not None and k < len(m["biomarkers"]) else None
            if bm is not None and bm["reaction"] != b["reaction"]:
                raise ValueError(f"{r['iem']} position {k}: {b['reaction']} != {bm['reaction']}")
            rows.append({"iem": r["iem"], "call_index": r["call_index"], "position": k, "reaction": b["reaction"],
                         "expected": b["expected"], "biofluid": biofluid(b["reaction"]),
                         "h_max": b.get("healthy"), "d_max": b.get("disease"),
                         "status_max": (b.get("status_healthy"), b.get("status_disease")),
                         "h_min": bm.get("healthy_min") if bm else None, "d_min": bm.get("disease_min") if bm else None,
                         "status_min": (bm.get("status_healthy_min"), bm.get("status_disease_min")) if bm else (None, None)})
    return rows


def analyse(rows: list) -> dict:
    directional = [r for r in rows if r["expected"] in ("Increased", "Decreased")]
    complete = [r for r in directional
                if r["status_max"] == ("Optimal", "Optimal") and r["status_min"] == ("Optimal", "Optimal")
                and all(_finite(r[k]) for k in ("h_max", "d_max", "h_min", "d_min"))]
    out_rows = []
    counts = {"max": 0, "R1": 0, "R2": 0}
    transitions = {rule: {"wrong_to_right": 0, "right_to_wrong": 0, "wrong_to_other_wrong": 0} for rule in ("R1", "R2")}
    for r in complete:
        calls = {"max": max_call(r["h_max"], r["d_max"]),
                 "R1": range_call(r["h_min"], r["h_max"], r["d_min"], r["d_max"]),
                 "R2": tie_break_call(r["h_min"], r["h_max"], r["d_min"], r["d_max"])}
        ok = {k: v == r["expected"] for k, v in calls.items()}
        for k in counts:
            counts[k] += ok[k]
        for rule in ("R1", "R2"):
            if not ok["max"] and ok[rule]:
                transitions[rule]["wrong_to_right"] += 1
            elif ok["max"] and not ok[rule]:
                transitions[rule]["right_to_wrong"] += 1
            elif not ok["max"] and not ok[rule] and calls[rule] != calls["max"]:
                transitions[rule]["wrong_to_other_wrong"] += 1
        capped = calls["max"] == "Unchanged" and r["h_max"] > TOL and abs(r["d_max"] - r["h_max"]) <= TOL
        out_rows.append({**{k: r[k] for k in ("iem", "call_index", "position", "reaction", "expected", "biofluid",
                                               "h_min", "h_max", "d_min", "d_max")},
                         "call_max": calls["max"], "call_R1": calls["R1"], "call_R2": calls["R2"],
                         "equal_positive_maxima": capped})
    n = len(complete)
    by_fluid = {}
    for row in out_rows:
        f = by_fluid.setdefault(row["biofluid"], {"n": 0, "max": 0, "R1": 0, "R2": 0})
        f["n"] += 1
        for k in ("max", "R1", "R2"):
            f[k] += row[f"call_{k}"] == row["expected"]
    capped_rows = [row for row in out_rows if row["equal_positive_maxima"]]
    reading = {}
    for rule in ("R1", "R2"):
        gain = counts[rule] - counts["max"]
        reading[rule] = "improves" if gain > 0 else "worsens" if gain < 0 else "no net change"
    return {
        "n_biomarkers": len(rows), "n_directional": len(directional), "n_scored_four_optima": n,
        "n_directional_excluded_missing_or_nonoptimal": len(directional) - n,
        "correct": counts, "accuracy": {k: (v / n if n else None) for k, v in counts.items()},
        "net_gain_vs_max": {rule: counts[rule] - counts["max"] for rule in ("R1", "R2")},
        "reading": reading, "transitions_vs_max": transitions,
        "n_conflicting_R1": sum(row["call_R1"] == "Conflicting" for row in out_rows),
        "equal_positive_maxima": {"n": len(capped_rows),
                                  "R1_correct": sum(row["call_R1"] == row["expected"] for row in capped_rows),
                                  "R2_correct": sum(row["call_R2"] == row["expected"] for row in capped_rows),
                                  "min_differs": sum(_sign(row["d_min"] - row["h_min"]) != 0 for row in capped_rows)},
        "minima_differ": sum(_sign(row["d_min"] - row["h_min"]) != 0 for row in out_rows),
        "by_biofluid": by_fluid,
        "changed_calls": [row for row in out_rows if row["call_R1"] != row["call_max"] or row["call_R2"] != row["call_max"]],
        "rows": out_rows,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="label=max_results.json,min_results.json")
    ap.add_argument("--out")
    args = ap.parse_args()
    report = {"rules": {"tolerance": TOL, "R1": "full range", "R2": "tie-break on equal maxima"}, "models": {}}
    for spec in args.runs:
        label, files = spec.split("=", 1)
        max_file, min_file = files.split(",")
        with open(max_file) as fh:
            max_results = json.load(fh)
        with open(min_file) as fh:
            min_results = json.load(fh)
        result = analyse(pair_records(max_results, min_results))
        result["files"] = {"max": max_file, "min": min_file}
        report["models"][label] = result
        print(label, json.dumps({k: result[k] for k in ("n_scored_four_optima", "correct", "net_gain_vs_max", "reading",
                                                        "transitions_vs_max", "n_conflicting_R1", "equal_positive_maxima")}))
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(report, fh, indent=1)
            fh.write("\n")


if __name__ == "__main__":
    main()
