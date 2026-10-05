"""Error anatomy of an IEM results file (plan v0.4, analysis A3).

Among scored biomarkers (expected direction, both maxima optimal and finite), splits the errors into "opposite
direction" and "no change predicted", and counts the no-change errors whose healthy and disease maxima are equal
and positive (|d - h| <= 1e-6 and h > 1e-6): values at a cap, which a maximum-only comparison cannot separate.

Usage: python scripts/iem_error_anatomy.py label=results.json [label=results.json ...] [--out f.json]
"""
from __future__ import annotations

import argparse
import json
import math

TOL = 1e-6


def finite(v) -> bool:
    return isinstance(v, (int, float)) and math.isfinite(v)


def anatomy(results: list) -> dict:
    scored = [(r["iem"], b) for r in results for b in r["biomarkers"]
              if b.get("correct") is not None and b.get("status_healthy") == b.get("status_disease") == "Optimal"
              and finite(b.get("healthy")) and finite(b.get("disease"))]
    wrong = [(iem, b) for iem, b in scored if not b["correct"]]
    opposite = [(iem, b) for iem, b in wrong if b["predicted"] != "Unchanged"]
    no_change = [(iem, b) for iem, b in wrong if b["predicted"] == "Unchanged"]
    capped = [(iem, b) for iem, b in no_change if abs(b["disease"] - b["healthy"]) <= TOL and b["healthy"] > TOL]
    zero = [(iem, b) for iem, b in no_change if abs(b["healthy"]) <= TOL and abs(b["disease"]) <= TOL]
    iems = {}
    for r in results:
        bs = r["biomarkers"]
        iems[(r["iem"], r["call_index"])] = bool(bs) and all(
            b.get("correct") is True and b.get("status_healthy") == b.get("status_disease") == "Optimal" for b in bs)
    return {"n_scored": len(scored), "n_correct": len(scored) - len(wrong), "n_errors": len(wrong),
            "n_opposite": len(opposite), "n_no_change": len(no_change),
            "n_no_change_equal_positive_maxima": len(capped), "n_no_change_zero_both": len(zero),
            "n_iems_fully_correct": sum(iems.values()),
            "opposite": [{"iem": i, "reaction": b["reaction"], "expected": b["expected"], "healthy": b["healthy"],
                          "disease": b["disease"]} for i, b in opposite],
            "no_change": [{"iem": i, "reaction": b["reaction"], "expected": b["expected"], "healthy": b["healthy"],
                           "disease": b["disease"]} for i, b in no_change]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out")
    args = ap.parse_args()
    report = {}
    for spec in args.runs:
        label, path = spec.split("=", 1)
        with open(path) as fh:
            report[label] = anatomy(json.load(fh))
        print(label, {k: v for k, v in report[label].items() if k.startswith("n_")})
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(report, fh, indent=1)
            fh.write("\n")


if __name__ == "__main__":
    main()
