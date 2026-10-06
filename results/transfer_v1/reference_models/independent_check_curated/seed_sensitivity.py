"""Seed sensitivity of the gene-bootstrap 95% intervals (curated - arm), as scored and with the corrected gene map.

Usage: python -I seed_sensitivity.py   (writes seed_sensitivity.json)
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import OUT, dump  # noqa: E402
from recompute_comparison import DEV, REF, counts, grows, load, mcc, table  # noqa: E402


def intervals(a, b, mode, seeds=range(50), n=1000):
    genes = sorted(set(a["browser_genes"]) & set(b["browser_genes"])) if mode == "common" else \
        sorted(set(a["browser_genes"]) | set(b["browser_genes"]))
    conds = sorted(grows(a) & grows(b))
    (sa, sb), fit = table([a, b], genes, conds, absent_as_wt=(mode == "union"))
    ca, cb = counts(sa, fit), counts(sb, fit)
    lows, highs = [], []
    for s in seeds:
        rng = np.random.default_rng(1000 + s)
        d = []
        for _ in range(n):
            w = np.bincount(rng.integers(0, len(ca), len(ca)), minlength=len(ca))
            d.append(mcc(w @ cb) - mcc(w @ ca))
        lo, hi = np.percentile(d, [2.5, 97.5])
        lows.append(float(lo)); highs.append(float(hi))
    lows, highs = np.array(lows), np.array(highs)
    return {"point": mcc(cb.sum(0)) - mcc(ca.sum(0)), "n_seeds": len(lows),
            "lower_min_max": [float(lows.min()), float(lows.max())], "upper_min_max": [float(highs.min()), float(highs.max())],
            "seeds_excluding_zero": int(((lows > 0) | (highs < 0)).sum())}


def main():
    arms = {"B0": load(os.path.join(DEV, "B0", "matrices.npz")), "U'": load(os.path.join(DEV, "UNQ", "matrices.npz")),
            "M": load(os.path.join(DEV, "M", "matrices.npz"))}
    out = {}
    for tag, path in (("as_scored", REF), ("corrected_gene_map", os.path.join(OUT, "matrices_iSO783_corrected_gene_map.npz"))):
        cur = load(path)
        out[tag] = {n: {m: intervals(arms[n], cur, m) for m in ("common", "union")} for n in arms}
        for n in arms:
            for m in ("common", "union"):
                x = out[tag][n][m]
                print(tag, n, m, round(x["point"], 4), "lower", [round(v, 4) for v in x["lower_min_max"]],
                      "upper", [round(v, 4) for v in x["upper_min_max"]], "seeds excluding 0:", x["seeds_excluding_zero"], "/", x["n_seeds"], flush=True)
    dump(out, "seed_sensitivity.json")


if __name__ == "__main__":
    main()
