"""Post hoc: how much of the IEM accuracy rests on very small healthy-versus-disease differences?

The protocol (checkIEM_WBM) calls a biomarker Increased or Decreased when the disease maximum differs from the
healthy maximum by more than 1e-6 (absolute). This descriptive analysis, which is not part of the frozen protocol,
recomputes the calls with a relative threshold tau: Increased if (d - h) / max(|h|, |d|) > tau, Decreased if
< -tau, otherwise Unchanged, for tau in {0 (protocol), 0.001, 0.01, 0.05, 0.10}. Only biomarkers scored by the
protocol (optimal, finite pairs with an expected direction) are used, so the denominators do not change.

Usage: python scripts/iem_effect_size_sensitivity.py label=results.json [label=results.json ...] [--out f.json]
"""
from __future__ import annotations

import argparse
import json
import math

TAUS = [0.0, 0.001, 0.01, 0.05, 0.10]


def call(h: float, d: float, tau: float, tol: float = 1e-6) -> str:
    diff = d - h
    if tau == 0.0:
        return "Increased" if diff > tol else "Decreased" if diff < -tol else "Unchanged"
    scale = max(abs(h), abs(d))
    rel = diff / scale if scale > 0 else 0.0
    return "Increased" if rel > tau else "Decreased" if rel < -tau else "Unchanged"


def scored(results):
    for r in results:
        for b in r["biomarkers"]:
            if (b["correct"] is not None and b["status_healthy"] == b["status_disease"] == "Optimal"
                    and isinstance(b["healthy"], (int, float)) and isinstance(b["disease"], (int, float))
                    and math.isfinite(b["healthy"]) and math.isfinite(b["disease"])):
                yield r["iem"], b


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out")
    args = ap.parse_args()
    report = {"note": "Post hoc, descriptive; not part of the frozen protocol.", "taus": TAUS, "runs": {}}
    for item in args.runs:
        label, path = item.split("=", 1)
        rows = list(scored(json.load(open(path))))
        res = {"n_scored": len(rows)}
        for tau in TAUS:
            correct = sum(call(b["healthy"], b["disease"], tau) == b["expected"] for _, b in rows)
            res[f"tau_{tau}"] = {"correct": correct, "accuracy": correct / len(rows) if rows else None}
        protocol_correct = [(i, b) for i, b in rows if b["correct"]]
        small = []
        for i, b in protocol_correct:
            scale = max(abs(b["healthy"]), abs(b["disease"]))
            rel = abs(b["disease"] - b["healthy"]) / scale if scale > 0 else 0.0
            if rel < 0.05:
                small.append({"iem": i, "biomarker": b["reaction"], "healthy": b["healthy"], "disease": b["disease"],
                              "relative_change": rel})
        res["protocol_correct_with_relative_change_below_5pct"] = sorted(small, key=lambda x: x["relative_change"])
        report["runs"][label] = res
        print(f"{label}: scored {len(rows)}; " + "; ".join(
            f"tau {tau:g}: {res[f'tau_{tau}']['correct']} ({res[f'tau_{tau}']['accuracy']:.3f})" for tau in TAUS)
              + f"; correct calls with relative change < 5%: {len(small)}")
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(report, fh, indent=1)
            fh.write("\n")


if __name__ == "__main__":
    main()
