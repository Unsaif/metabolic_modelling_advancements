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
from gembench.comparison import load_run, aligned_pairs  # noqa: E402


def load(d):
    z = np.load(os.path.join(d, "matrices.npz"), allow_pickle=True)
    return z["sim_growth"], z["fitness"], z["wt_growth"], list(z["browser_genes"]), list(z["conditions"])


def main() -> None:
    org, da, db = sys.argv[1:4]
    la, lb = (sys.argv[4], sys.argv[5]) if len(sys.argv) > 5 else (os.path.basename(da), os.path.basename(db))
    if la == lb:
        raise ValueError("Comparison labels must be distinct")
    a, b = load_run(da), load_run(db)
    sa, sb, fa, genes, grows, conds = aligned_pairs(a, b)
    fb = fa  # identical observations and shared finite mask, checked by aligned_pairs
    wa, wb = a['wt_growth'], b['wt_growth']
    st, ft = a['params']['growth_threshold'], a['params']['fitness_threshold']
    ia = ib = list(range(len(genes)))
    ja = jb = list(range(len(grows)))
    metrics = [("aucpr_bernstein", lambda s,f: M.aucpr_bernstein(s,f,st)),
               ("auroc_standard", lambda s,f: M.auroc_standard(s,f,ft)),
               ("mcc", lambda s,f: M.mcc(s,f,st,ft)),
               ("balanced_accuracy", lambda s,f: M.balanced_accuracy(s,f,st,ft))]
    out = {"org": org, "n_common_genes": len(genes), "n_common_conditions": len(conds), "n_conditions_both_grow": len(ja),
           "n_shared_finite_pairs": int(np.isfinite(fa).sum()), "evaluation_role": "retrospective_development",
           "wt_grows": {la: int((np.isfinite(wa) & (wa >= st)).sum()), lb: int((np.isfinite(wb) & (wb >= st)).sum())}, "models": {}}
    for label, s, f in [(la, sa[np.ix_(ia, ja)], fa[np.ix_(ia, ja)]), (lb, sb[np.ix_(ib, jb)], fb[np.ix_(ib, jb)])]:
        res = {}
        for name, fn in metrics:
            pt, lo, hi = M.bootstrap_ci(fn, s, f, n_boot=500)
            res[name] = {"point": pt, "ci95": [lo, hi]}
        res["confusion"] = M.confusion(s, f, st, ft)
        out["models"][label] = res
    # paired bootstrap on the difference in MCC and AUC-PR (same gene resamples for both models)
    rng = np.random.default_rng(0)
    A = (sa[np.ix_(ia, ja)], fa[np.ix_(ia, ja)]); B = (sb[np.ix_(ib, jb)], fb[np.ix_(ib, jb)])
    diffs = {"mcc": [], "aucpr_bernstein": []}
    n = len(genes)
    for _ in range(500 if n else 0):
        idx = rng.integers(0, n, n)
        for k, fn in [(k, fn) for k, fn in metrics if k in diffs]:
            try:
                diffs[k].append(fn(B[0][idx], B[1][idx]) - fn(A[0][idx], A[1][idx]))
            except Exception:  # noqa: BLE001
                pass
    out["paired_difference_B_minus_A"] = {k: {"point": out["models"][lb][k]["point"] - out["models"][la][k]["point"],
                                             "ci95": [float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))] if np.isfinite(v).any() else [float("nan"), float("nan")]}
                                         for k, v in diffs.items()}
    print(json.dumps(out, indent=2))
    with open(os.path.join(ROOT, "results", "carbon_fitness_multi", f"{org}_paired_{la}_vs_{lb}.json"), "w") as fh:
        json.dump(out, fh, indent=2)


if __name__ == "__main__":
    main()
