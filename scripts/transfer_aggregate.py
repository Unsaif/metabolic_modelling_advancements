"""Organism-level aggregate of paired arm comparisons (transfer study v1, docs/studies/transfer-method-v1.md section 5).

Input: the JSON written by scripts/transfer_compare.py --out (one row per organism x pair). For each pair A:B the
script reports every organism's paired MCC difference (or why it is not evaluable), the unweighted mean over evaluable
organisms with a percentile organism-level bootstrap interval (10,000 resamples, seed 0), the numbers of organisms
improved / worsened / unchanged (|difference| < 0.001), and the exact two-sided sign-test p-value over organisms
that changed. Each organism counts once, however many genes and conditions it contributes.

The primary metric is the union-of-genes paired difference (see scripts/transfer_compare.py); --metric common gives the
common-gene difference.

Usage: python scripts/transfer_aggregate.py results/transfer_v1/<role>/paired.json [--metric union|common]
       [--exclude MycoTube] [--out f.json]
"""
from __future__ import annotations

import argparse
import json
import math
from collections import OrderedDict

import numpy as np

UNCHANGED = 1e-3


def sign_test_p(n_pos: int, n_neg: int) -> float:
    n = n_pos + n_neg
    if n == 0:
        return float("nan")
    k = min(n_pos, n_neg)
    p = sum(math.comb(n, i) for i in range(0, k + 1)) / 2 ** n
    return min(1.0, 2 * p)


METRICS = {"union": ("union_delta", "union_delta_ci95"), "common": ("delta", "delta_ci95")}


def aggregate(rows, n_boot=10000, seed=0, exclude=(), metric="union"):
    key, ci_key = METRICS[metric]
    pairs = OrderedDict()
    for r in rows:
        pairs.setdefault((r["A"], r["B"]), []).append(r)
    out = []
    for (a, b), rs in pairs.items():
        per_org, deltas = [], []
        for r in rs:
            if r["org"] in exclude:
                continue
            d = r.get(key)
            ok = (not r.get("missing")) and d is not None and isinstance(d, (int, float)) and math.isfinite(d)
            per_org.append({"org": r["org"], "delta": d if ok else None,
                            "ci95": r.get(ci_key) if ok else None,
                            "wt_grows": [r.get("wt_grows_A"), r.get("wt_grows_B")],
                            "n_mapped_conditions": r.get("n_mapped_conditions"),
                            "n_pairs": r.get("n_pairs"),
                            "status": "evaluable" if ok else ("missing run" if r.get("missing") else "undefined")})
            if ok:
                deltas.append(d)
        x = np.array(deltas, dtype=float)
        res = {"A": a, "B": b, "n_organisms": len(per_org), "n_evaluable": int(len(x)), "per_organism": per_org}
        if len(x):
            rng = np.random.default_rng(seed)
            boots = x[rng.integers(0, len(x), size=(n_boot, len(x)))].mean(axis=1)
            pos, neg = int((x >= UNCHANGED).sum()), int((x <= -UNCHANGED).sum())
            res.update({"mean_delta": float(x.mean()), "median_delta": float(np.median(x)),
                        "organism_bootstrap_ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
                        "n_improved": pos, "n_worsened": neg, "n_unchanged": int(len(x) - pos - neg),
                        "sign_test_p_two_sided": sign_test_p(pos, neg)})
        out.append(res)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("paired_json")
    ap.add_argument("--exclude", default="")
    ap.add_argument("--metric", choices=sorted(METRICS), default="union",
                    help="union of both arms' genes (primary) or genes common to both (secondary)")
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--out")
    a = ap.parse_args()
    rows = json.load(open(a.paired_json))
    excl = tuple(x for x in a.exclude.split(",") if x)
    res = aggregate(rows, n_boot=a.n_boot, exclude=excl, metric=a.metric)
    print(f"metric: {a.metric}")
    for r in res:
        if "mean_delta" in r:
            lo, hi = r["organism_bootstrap_ci95"]
            print(f"{r['A']:>10s} -> {r['B']:<10s} evaluable {r['n_evaluable']}/{r['n_organisms']}: mean {r['mean_delta']:+.4f} "
                  f"[{lo:+.4f}, {hi:+.4f}] improved {r['n_improved']} worsened {r['n_worsened']} unchanged {r['n_unchanged']} "
                  f"sign p={r['sign_test_p_two_sided']:.3f}  " +
                  " ".join(f"{o['org']}={o['delta']:+.4f}" if o["delta"] is not None else f"{o['org']}=n/a" for o in r["per_organism"]))
        else:
            print(f"{r['A']:>10s} -> {r['B']:<10s} no evaluable organism")
    if a.out:
        json.dump({"metric": a.metric, "excluded": list(excl), "n_boot": a.n_boot, "seed": 0, "unchanged_tolerance": UNCHANGED, "pairs": res},
                  open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
