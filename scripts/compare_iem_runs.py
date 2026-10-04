"""Compare whole-body IEM runs biomarker by biomarker.

Usage:
  python scripts/compare_iem_runs.py --runs v0.2=results/wbm_iem/Harvey_1_03d_iem_results_v0.2.json \
      v0.2b=results/wbm_iem/Harvey_1_03d_iem_results_v0.2b.json v0.3=results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json \
      --protocol data/iem/iem_protocol_v0.2.json --out results/wbm_iem/Harvey_1_03d_iem_comparison_v0.3

Writes <out>.json (per-run summaries and pairwise transitions between consecutive runs) and <out>.tsv
(one row per protocol biomarker with each run's healthy and disease values and call).
A biomarker is scored when its expected direction is Increased or Decreased and both maximisations
returned an optimal, finite value (as in run_wbm_iem.scored_biomarkers).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from collections import Counter


def finite(v):
    return isinstance(v, (int, float)) and math.isfinite(v)


def index_run(results):
    out = {}
    for r in results:
        for b in r["biomarkers"]:
            key = (r["iem"], r["call_index"], b["reaction"])
            scored = (b["correct"] is not None and b["status_healthy"] == b["status_disease"] == "Optimal"
                      and finite(b["healthy"]) and finite(b["disease"]))
            out[key] = {"expected": b["expected"], "predicted": b["predicted"], "healthy": b["healthy"],
                        "disease": b["disease"], "scored": scored, "correct": b["correct"] if scored else None,
                        "iem_status": r["status"]}
    return out


def error_kind(rec):
    if not rec["scored"] or rec["correct"]:
        return None
    return "no change predicted" if rec["predicted"] == "Unchanged" else "opposite direction"


def summarise(protocol, results, idx):
    statuses = Counter(r["status"] for r in results)
    scored = [v for v in idx.values() if v["scored"]]
    correct = sum(v["correct"] for v in scored)
    expected = sum(len(p["biomarkers"]) for p in protocol)
    errors = Counter(error_kind(v) for v in scored if not v["correct"])
    per_iem = {}
    for p in protocol:
        recs = [idx.get((p["iem"], p["call_index"], b[0])) for b in p["biomarkers"]]
        s = [r for r in recs if r and r["scored"]]
        per_iem[f'{p["iem"]}#{p["call_index"]}'] = {"scored": len(s), "correct": sum(r["correct"] for r in s),
                                                     "expected": len(p["biomarkers"])}
    return {"n_iems": len(results), "iem_status": dict(statuses), "n_biomarkers_in_protocol": expected,
            "n_scored": len(scored), "n_correct": correct, "accuracy_among_scored": correct / len(scored) if scored else None,
            "coverage": len(scored) / expected if expected else None, "errors": dict(errors),
            "n_iems_all_scored_correct": sum(1 for v in per_iem.values() if v["scored"] == v["expected"] and v["correct"] == v["scored"]),
            "n_iems_all_correct_among_scored": sum(1 for v in per_iem.values() if v["scored"] and v["correct"] == v["scored"]),
            "per_iem": per_iem}


def transitions(protocol, a, b):
    t = Counter()
    changed = []
    for p in protocol:
        for bm in p["biomarkers"]:
            key = (p["iem"], p["call_index"], bm[0])
            ra, rb = a.get(key), b.get(key)
            sa = "unscored" if not ra or not ra["scored"] else ("right" if ra["correct"] else "wrong")
            sb = "unscored" if not rb or not rb["scored"] else ("right" if rb["correct"] else "wrong")
            t[f"{sa}->{sb}"] += 1
            pa = ra["predicted"] if ra else "absent"
            pb = rb["predicted"] if rb else "absent"
            if pa != pb:
                changed.append({"iem": p["iem"], "call_index": p["call_index"], "biomarker": bm[0], "expected": bm[1],
                                "from": pa, "to": pb, "from_status": sa, "to_status": sb})
    return {"counts": dict(t), "n_prediction_changes": len(changed), "changed": changed}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", required=True, help="label=path, in comparison order")
    ap.add_argument("--protocol", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    protocol = json.load(open(args.protocol))
    runs = []
    for item in args.runs:
        label, path = item.split("=", 1)
        results = json.load(open(path))
        fps = {r.get("run_fingerprint") for r in results}
        if len(fps) != 1:
            raise SystemExit(f"{path}: results from more than one run fingerprint")
        runs.append((label, path, results, index_run(results)))
    report = {"protocol": args.protocol, "runs": {}, "transitions": {}}
    for label, path, results, idx in runs:
        report["runs"][label] = {"path": path, "run_fingerprint": results[0]["run_fingerprint"],
                                 **summarise(protocol, results, idx)}
    for (la, _, _, ia), (lb, _, _, ib) in zip(runs, runs[1:]):
        report["transitions"][f"{la}->{lb}"] = transitions(protocol, ia, ib)
    if len(runs) > 2:
        la, _, _, ia = runs[0]
        lb, _, _, ib = runs[-1]
        report["transitions"][f"{la}->{lb}"] = transitions(protocol, ia, ib)
    with open(args.out + ".json", "w") as fh:
        json.dump(report, fh, indent=1, allow_nan=False)
        fh.write("\n")
    with open(args.out + ".tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        head = ["iem", "call_index", "biomarker", "expected"]
        for label, *_ in runs:
            head += [f"{label}_healthy", f"{label}_disease", f"{label}_call", f"{label}_correct"]
        w.writerow(head)
        for p in protocol:
            for bm in p["biomarkers"]:
                key = (p["iem"], p["call_index"], bm[0])
                row = [p["iem"], p["call_index"], bm[0], bm[1]]
                for _, _, _, idx in runs:
                    r = idx.get(key)
                    if r is None:
                        row += ["", "", "absent", ""]
                    else:
                        row += [r["healthy"], r["disease"], r["predicted"], "" if r["correct"] is None else int(r["correct"])]
                w.writerow(row)
    for label, s in report["runs"].items():
        print(f"{label:6s} IEMs {s['n_iems']:2d} scored {s['n_scored']}/{s['n_biomarkers_in_protocol']} correct {s['n_correct']} "
              f"acc {s['accuracy_among_scored']:.3f} errors {s['errors']} all-correct IEMs {s['n_iems_all_correct_among_scored']}")
    for k, t in report["transitions"].items():
        print(f"{k}: {t['counts']} prediction changes {t['n_prediction_changes']}")


if __name__ == "__main__":
    main()
