"""Per-disease breakdown of the own-biomarker calls (post hoc, descriptive; no new solves).

Usage: python scripts/iem_per_disease_calls.py

For each IEM and model run, every biomarker with an expected direction is counted as
  correct     predicted change in the expected direction
  missed      predicted Unchanged
  opposite    predicted change in the opposite direction
  unavailable no optimum, or the reaction is absent from the model
Per disease: recall = correct / (correct + missed + opposite); precision = correct / (correct + opposite), i.e. how
often a predicted change has the right direction. Unavailable biomarkers are reported, not scored. Changes the model
predicts on readouts outside a disease's biomarker list are not judged here: absence from the list is not evidence of
no change (the disease-ranking study looks at them).
"""
from __future__ import annotations

import csv
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = {"Harvey_v0.3": "results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json",
        "Harvetta_v0.4": "results/wbm_iem/Harvetta_1_03d_iem_results_v0.4.json"}
OUT = os.path.join(ROOT, "results", "wbm_iem", "per_disease_calls_v0.3_v0.4")


def classify(b):
    if b["expected"] not in ("Increased", "Decreased"):
        return None
    if b["status_healthy"] != "Optimal" or b["status_disease"] != "Optimal" or b["predicted"] == "NA":
        return "unavailable"
    if b["predicted"] == b["expected"]:
        return "correct"
    if b["predicted"] == "Unchanged":
        return "missed"
    return "opposite"


def main():
    os.makedirs(OUT, exist_ok=True)
    report, rows = {"note": "Post hoc, descriptive. See the module docstring for definitions.", "runs": {}}, []
    for run, path in RUNS.items():
        with open(os.path.join(ROOT, path)) as fh:
            results = sorted(json.load(fh), key=lambda r: r["call_index"])
        per = []
        for r in results:
            counts = {"correct": 0, "missed": 0, "opposite": 0, "unavailable": 0}
            for b in r["biomarkers"]:
                c = classify(b)
                if c:
                    counts[c] += 1
            scored = counts["correct"] + counts["missed"] + counts["opposite"]
            changed = counts["correct"] + counts["opposite"]
            rec = {"iem": r["iem"], "call_index": r["call_index"], "n_biomarkers": len(r["biomarkers"]), **counts,
                   "recall": counts["correct"] / scored if scored else None,
                   "precision": counts["correct"] / changed if changed else None}
            per.append(rec)
            rows.append({"run": run, **rec})
        def n(cond):
            return sum(1 for x in per if cond(x))
        totals = {k: sum(x[k] for x in per) for k in ("correct", "missed", "opposite", "unavailable")}
        report["runs"][run] = {
            "file": path, "n_iems": len(per), "totals": totals,
            "pooled_recall": totals["correct"] / (totals["correct"] + totals["missed"] + totals["opposite"]),
            "pooled_precision": totals["correct"] / (totals["correct"] + totals["opposite"]),
            "n_iems_all_scored_correct": n(lambda x: x["recall"] == 1.0),
            "n_iems_none_correct": n(lambda x: x["recall"] == 0.0),
            "n_iems_with_any_opposite": n(lambda x: x["opposite"] > 0),
            "n_iems_with_any_missed": n(lambda x: x["missed"] > 0),
            "median_recall": sorted(x["recall"] for x in per if x["recall"] is not None)[len(per) // 2],
            "per_disease": per}
        s = report["runs"][run]
        print(f"{run}: pooled recall {s['pooled_recall']:.3f}, precision {s['pooled_precision']:.3f}; "
              f"all correct in {s['n_iems_all_scored_correct']}/{len(per)} IEMs, none correct in {s['n_iems_none_correct']}; "
              f"any opposite call in {s['n_iems_with_any_opposite']}, any missed change in {s['n_iems_with_any_missed']}; "
              f"totals {totals}")
    with open(os.path.join(OUT, "per_disease_calls.json"), "w") as fh:
        json.dump(report, fh, indent=1)
    with open(os.path.join(OUT, "per_disease_calls.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
