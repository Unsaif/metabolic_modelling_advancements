"""Paired comparison of two result sets for the same organism on the intersection of genes and of
conditions where both wild types grow. Pooled metrics over different gene sets are not comparable
(the extra genes one model covers are not a random sample), so this is the comparison to quote.

Usage: python scripts/compare_models_paired.py <org> <results_dir_A> <results_dir_B> [label_A label_B]
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from gembench import metrics as M  # noqa: E402


def load(d):
    z = np.load(os.path.join(d, "matrices.npz"), allow_pickle=True)
    return z["sim_growth"], z["fitness"], z["wt_growth"], list(z["browser_genes"]), list(z["conditions"])


def main() -> None:
    org, da, db = sys.argv[1:4]
    la, lb = (sys.argv[4], sys.argv[5]) if len(sys.argv) > 5 else (os.path.basename(da), os.path.basename(db))
    sa, fa, wa, ga, ca = load(da); sb, fb, wb, gb, cb = load(db)
    genes = [g for g in ga if g in set(gb)]
    conds = [c for c in ca if c in set(cb)]
    ia = [ga.index(g) for g in genes]; ib = [gb.index(g) for g in genes]
    ja, jb = [], []
    for c in conds:
        x, y = ca.index(c), cb.index(c)
        if wa[x] >= 1e-3 and wb[y] >= 1e-3:
            ja.append(x); jb.append(y)
    out = {"org": org, "n_common_genes": len(genes), "n_common_conditions": len(conds), "n_conditions_both_grow": len(ja),
           "wt_grows": {la: int((wa >= 1e-3).sum()), lb: int((wb >= 1e-3).sum())}, "models": {}}
    for label, s, f in [(la, sa[np.ix_(ia, ja)], fa[np.ix_(ia, ja)]), (lb, sb[np.ix_(ib, jb)], fb[np.ix_(ib, jb)])]:
        res = {}
        for name, fn in [("aucpr_bernstein", M.aucpr_bernstein), ("auroc_standard", M.auroc_standard), ("mcc", M.mcc),
                         ("balanced_accuracy", M.balanced_accuracy)]:
            pt, lo, hi = M.bootstrap_ci(fn, s, f, n_boot=500)
            res[name] = {"point": pt, "ci95": [lo, hi]}
        res["confusion"] = M.confusion(s, f)
        out["models"][label] = res
    # paired bootstrap on the difference in MCC and AUC-PR (same gene resamples for both models)
    rng = np.random.default_rng(0)
    A = (sa[np.ix_(ia, ja)], fa[np.ix_(ia, ja)]); B = (sb[np.ix_(ib, jb)], fb[np.ix_(ib, jb)])
    diffs = {"mcc": [], "aucpr_bernstein": []}
    n = len(genes)
    for _ in range(500):
        idx = rng.integers(0, n, n)
        for k, fn in [("mcc", M.mcc), ("aucpr_bernstein", M.aucpr_bernstein)]:
            try:
                diffs[k].append(fn(B[0][idx], B[1][idx]) - fn(A[0][idx], A[1][idx]))
            except Exception:  # noqa: BLE001
                pass
    out["paired_difference_B_minus_A"] = {k: {"point": out["models"][lb][k]["point"] - out["models"][la][k]["point"],
                                             "ci95": [float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))]}
                                         for k, v in diffs.items()}
    print(json.dumps(out, indent=2))
    with open(os.path.join(ROOT, "results", "carbon_fitness_multi", f"{org}_paired_{la}_vs_{lb}.json"), "w") as fh:
        json.dump(out, fh, indent=2)


if __name__ == "__main__":
    main()
