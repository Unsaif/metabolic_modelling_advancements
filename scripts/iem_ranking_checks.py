"""Checks on the disease-ranking prediction matrix before any ranking (plan: "Checks before any ranking").

Usage: python scripts/iem_ranking_checks.py [results/wbm_iem/Harvey_1_03d_iem_cross_ranking_v1.json]

Reports, without printing any cross-disease value:
  1. completeness: IEM statuses, non-optimal solves, fallbacks, NA and absent readouts;
  2. own biomarkers against v0.3: call agreement (gate: more than 5 differing calls stops the analysis) and value
     differences;
  3. rechecks: barrier re-solves from scratch against the warm solves;
  4. the HIS readouts against the one-IEM-at-a-time Gurobi feasibility run (call agreement, maximum difference);
  5. solve times.
Writes results/wbm_iem/ranking/<matrix stem>_checks.json.
"""
from __future__ import annotations

import json
import math
import os
import statistics
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
V03 = os.path.join(ROOT, "results", "wbm_iem", "Harvey_1_03d_iem_results_v0.3.json")
FEAS = os.path.join(ROOT, "results", "wbm_iem", "feasibility", "Harvey_1_03d_iem_cross_feasibility_gurobi_primal.json")
GATE = 5


def finite(x):
    return x is not None and isinstance(x, (int, float)) and math.isfinite(x)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results", "wbm_iem", "Harvey_1_03d_iem_cross_ranking_v1.json")
    with open(path) as fh:
        recs = json.load(fh)
    with open(V03) as fh:
        v03 = {r["iem"]: {b["reaction"]: b for b in r["biomarkers"]} for r in json.load(fh)}
    out = {"matrix": os.path.relpath(path, ROOT)}

    # 1. completeness
    statuses = {}
    non_optimal, fallbacks, na, absent = [], 0, {}, set()
    for r in recs:
        statuses[r["status"]] = statuses.get(r["status"], 0) + 1
        for e in r["readouts"]:
            for st in ("healthy", "disease"):
                s = e[f"status_{st}"]
                if s == "absent":
                    absent.add(e["reaction"])
                elif s != "Optimal":
                    non_optimal.append((r["iem"], e["reaction"], st, s))
                fallbacks += bool(e[f"fallback_{st}"])
            if e["predicted"] == "NA" and e["status_healthy"] != "absent":
                na.setdefault(e["reaction"], []).append(r["iem"])
    n_readouts = {len(r["readouts"]) for r in recs}
    out["completeness"] = {"n_iems": len(recs), "statuses": statuses, "readouts_per_iem": sorted(n_readouts),
                           "absent_readouts": sorted(absent), "non_optimal_solves": non_optimal,
                           "fallbacks": fallbacks, "na_readouts_other_than_absent": na}

    # 2. own biomarkers against v0.3
    rows = []
    for r in recs:
        for e in r["readouts"]:
            if not e["own"]:
                continue
            b = v03[r["iem"]].get(e["reaction"])
            if b is None:
                continue
            d = [abs(x - y) for x, y in ((b["healthy"], e["healthy"]), (b["disease"], e["disease"])) if finite(x) and finite(y)]
            rel = [abs(x - y) / max(abs(x), abs(y)) for x, y in ((b["healthy"], e["healthy"]), (b["disease"], e["disease"]))
                   if finite(x) and finite(y) and max(abs(x), abs(y)) > 1e-3]
            rows.append({"iem": r["iem"], "reaction": e["reaction"], "v03_call": b["predicted"], "call": e["predicted"],
                         "max_abs_diff": max(d) if d else None, "max_rel_diff_above_1e-3": max(rel) if rel else None})
    differing = [x for x in rows if x["v03_call"] != x["call"]]
    big = sorted((x for x in rows if (x["max_rel_diff_above_1e-3"] or 0) > 1e-4), key=lambda x: -x["max_rel_diff_above_1e-3"])
    out["own_biomarkers_vs_v03"] = {
        "n": len(rows), "n_same_call": len(rows) - len(differing), "differing_calls": differing,
        "gate_passed": len(differing) <= GATE,
        "max_abs_diff": max((x["max_abs_diff"] or 0) for x in rows),
        "max_rel_diff_values_above_1e-3": max((x["max_rel_diff_above_1e-3"] or 0) for x in rows),
        "values_with_rel_diff_above_1e-4": [(x["iem"], x["reaction"], x["max_rel_diff_above_1e-3"]) for x in big]}

    # 3. rechecks
    checks = [(r["iem"], e["reaction"], st, e[f"solve_{st}"]["recheck"]) for r in recs for e in r["readouts"]
              for st in ("healthy", "disease") if "recheck" in (e[f"solve_{st}"] or {})]
    out["rechecks"] = {"n": len(checks), "non_optimal": sum(c["status"] != "Optimal" for *_, c in checks),
                       "max_abs_diff": max([c["abs_diff"] or 0 for *_, c in checks] + [0]),
                       "n_abs_diff_above_1e-6": sum((c["abs_diff"] or 0) > 1e-6 for *_, c in checks)}

    # 4. HIS against the one-IEM-at-a-time feasibility run (values compared by program only)
    if os.path.exists(FEAS):
        with open(FEAS) as fh:
            feas = {e["reaction"]: e for e in json.load(fh)[0]["readouts"]}
        his = next(r for r in recs if r["iem"] == "HIS")
        same = n = 0
        worst = 0.0
        for e in his["readouts"]:
            f = feas.get(e["reaction"])
            if f is None or e["status_healthy"] == "absent":
                continue
            n += 1
            same += f["predicted"] == e["predicted"]
            for x, y in ((f["healthy"], e["healthy"]), (f["disease"], e["disease"])):
                if finite(x) and finite(y):
                    worst = max(worst, abs(x - y))
        out["his_vs_feasibility_run"] = {"n_readouts_compared": n, "n_same_call": same, "max_abs_diff": worst}

    # 5. timing
    times = {}
    for r in recs:
        for e in r["readouts"]:
            for st in ("healthy", "disease"):
                m = (e[f"solve_{st}"] or {}).get("method")
                if m:
                    times.setdefault(m, []).append(e[f"time_{st}_s"])
    out["timing"] = {m: {"n": len(t), "median_s": statistics.median(t), "total_h": sum(t) / 3600} for m, t in times.items()}

    stem = os.path.splitext(os.path.basename(path))[0].replace("_iem_cross_ranking", "_ranking")
    dest = os.path.join(ROOT, "results", "wbm_iem", "ranking", f"{stem}_checks.json")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=1)
    c, o, rc = out["completeness"], out["own_biomarkers_vs_v03"], out["rechecks"]
    print(f"completeness: {c['n_iems']} IEMs, statuses {c['statuses']}, readouts per IEM {c['readouts_per_iem']}, "
          f"absent {c['absent_readouts']}, non-optimal {len(c['non_optimal_solves'])}, fallbacks {c['fallbacks']}, "
          f"other NA readouts {len(c['na_readouts_other_than_absent'])}")
    print(f"own biomarkers vs v0.3: {o['n_same_call']}/{o['n']} same call (gate {'passed' if o['gate_passed'] else 'FAILED'}); "
          f"max abs diff {o['max_abs_diff']:.3g}; max rel diff (values > 1e-3) {o['max_rel_diff_values_above_1e-3']:.3g}; "
          f"{len(o['values_with_rel_diff_above_1e-4'])} values differ by more than 1e-4 relative")
    print(f"rechecks: {rc['n']}, non-optimal {rc['non_optimal']}, max abs diff {rc['max_abs_diff']:.3g}")
    if "his_vs_feasibility_run" in out:
        h = out["his_vs_feasibility_run"]
        print(f"HIS vs feasibility run: {h['n_same_call']}/{h['n_readouts_compared']} same call, max abs diff {h['max_abs_diff']:.3g}")
    print("timing:", {m: (v["n"], round(v["median_s"], 2), round(v["total_h"], 1)) for m, v in out["timing"].items()})
    print("written", os.path.relpath(dest, ROOT))


if __name__ == "__main__":
    main()
