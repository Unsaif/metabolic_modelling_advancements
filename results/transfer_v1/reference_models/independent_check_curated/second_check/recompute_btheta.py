"""Check 3b: recompute the exploratory iAH991 entry of comparison.json from the saved matrices (own numpy code,
the same definitions as ../recompute_comparison.py), plus two sensitivity views:

  * conditions where iAH991's wild type grows only just above the 1e-3 threshold (L-fucose, L-rhamnose: 0.0012),
    where almost any knockout that lowers growth by a fifth becomes an "important" call;
  * the same comparison with the drafts re-scored with the same B12 and sulfide supplement
    (matrices from rescore_drafts_supplemented.py, when present).

Usage: python -I recompute_btheta.py [--drafts-dir <dir with B0.npz, UNQ.npz, M.npz>] [--tag name]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from common import ROOT  # noqa: E402
from recompute_comparison import FT, ST, boot, counts, grows, load, mcc, pair, table  # noqa: E402

DEV = os.path.join(ROOT, "results", "transfer_v1", "development", "Btheta")
CUR = os.environ.get("CUR_NPZ") or os.path.join(ROOT, "results", "transfer_v1", "reference_models", "Btheta", "iAH991_exploratory_B12_H2S", "matrices.npz")


def analyse(cur, arms, drop_conditions=()):
    if drop_conditions:
        def drop(r):
            keep = [j for j, c in enumerate(r["conditions"]) if c not in drop_conditions]
            return {**r, "sim_growth": r["sim_growth"][:, keep], "fitness": r["fitness"][:, keep],
                    "wt_growth": r["wt_growth"][keep], "conditions": [r["conditions"][j] for j in keep]}
        cur, arms = drop(cur), {k: drop(v) for k, v in arms.items()}
    rep = {}
    cols = [j for j, c in enumerate(cur["conditions"]) if c in grows(cur)]
    c = counts(cur["sim_growth"][:, cols], cur["fitness"][:, cols]).sum(0)
    rep["own_iAH991"] = {"n_genes": len(cur["browser_genes"]), "n_conditions": len(cols), "mcc": mcc(c),
                         "confusion_important_positive": dict(zip(["tp", "tn", "fp", "fn"], map(int, c)))}
    rep["own_arms"] = {}
    for k, r in arms.items():
        cc = [j for j, x in enumerate(r["conditions"]) if x in grows(r)]
        rep["own_arms"][k] = mcc(counts(r["sim_growth"][:, cc], r["fitness"][:, cc]).sum(0))
    rep["coverage"] = {k: len(grows(r)) for k, r in {"iAH991": cur, **arms}.items()}
    runs = [arms["B0"], arms["U'"], arms["M"], cur]
    genes = sorted(set().union(*[set(r["browser_genes"]) for r in runs]))
    conds = sorted(set.intersection(*[grows(r) for r in runs]))
    sims, fit = table(runs, genes, conds, absent_as_wt=True)
    rep["four_way"] = {"union_genes": len(genes), "conditions_all_grow": len(conds), "n_observations": int(np.isfinite(fit).sum()),
                       "mcc": {n: mcc(counts(s, fit).sum(0)) for n, s in zip(["B0", "U'", "M", "iAH991"], sims)}}
    rep["pairwise"] = {n: {"common": pair(arms[n], cur, "common"), "union": pair(arms[n], cur, "union")} for n in ("B0", "U'", "M")}
    u = arms["U'"]
    cg = sorted(set(u["browser_genes"]) & set(cur["browser_genes"]))
    cc = sorted(grows(u) & grows(cur))
    (su, sc), fc = table([u, cur], cg, cc, absent_as_wt=False)
    ok = np.isfinite(su) & np.isfinite(sc) & np.isfinite(fc)
    obs = fc <= FT
    wu, wc = ok & ((su < ST) != obs), ok & ((sc < ST) != obs)
    n = int(ok.sum())
    rep["error_overlap"] = {"n_genes": len(cg), "n_conditions": len(cc), "n_observations": n, "wrong_Uprime": int(wu.sum()),
                            "wrong_curated": int(wc.sum()), "wrong_both": int((wu & wc).sum()),
                            "expected_both_if_independent": float(wu.sum() * wc.sum() / n) if n else None,
                            "ratio_to_expected": float((wu & wc).sum() / (wu.sum() * wc.sum() / n)) if n else None,
                            "share_of_Uprime_errors_shared": float((wu & wc).sum() / wu.sum()) if wu.sum() else None,
                            "both_missed_important": int((wu & wc & (su >= ST) & (sc >= ST)).sum()),
                            "both_false_important": int((wu & wc & (su < ST) & (sc < ST)).sum()),
                            "Uprime_false_important": int((ok & (su < ST) & ~obs).sum()), "Uprime_missed": int((ok & (su >= ST) & obs).sum()),
                            "curated_false_important": int((ok & (sc < ST) & ~obs).sum()), "curated_missed": int((ok & (sc >= ST) & obs).sum())}
    draft = set().union(*[set(r["browser_genes"]) for r in arms.values()])
    only = [g for g in cur["browser_genes"] if g not in draft]
    gi = {g: i for i, g in enumerate(cur["browser_genes"])}
    cols = [j for j, x in enumerate(cur["conditions"]) if x in grows(cur) and x in set(cc)]
    so = cur["sim_growth"][np.ix_([gi[g] for g in only], cols)]
    fo = cur["fitness"][np.ix_([gi[g] for g in only], cols)]
    imp = np.isfinite(so) & np.isfinite(fo) & (so < ST)
    per_cond = {cur["conditions"][j].split(" |")[0]: int(imp[:, k].sum()) for k, j in enumerate(cols)}
    rep["curated_only"] = {"n_genes": len(only), "n_conditions": len(cols), "important_calls": int(imp.sum()),
                           "important_calls_confirmed": int((imp & (fo <= FT)).sum()),
                           "observed_important_cells": int((np.isfinite(fo) & (fo <= FT)).sum()),
                           "important_calls_by_condition": per_cond}
    return rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drafts-dir", default="")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    cur = load(CUR)
    if a.drafts_dir:
        arms = {"B0": load(os.path.join(a.drafts_dir, "B0.npz")), "U'": load(os.path.join(a.drafts_dir, "UNQ.npz")),
                "M": load(os.path.join(a.drafts_dir, "M.npz"))}
    else:
        arms = {"B0": load(os.path.join(DEV, "B0", "matrices.npz")), "U'": load(os.path.join(DEV, "UNQ", "matrices.npz")),
                "M": load(os.path.join(DEV, "M", "matrices.npz"))}
    rep = {"curated": os.path.relpath(CUR, ROOT), "drafts": a.drafts_dir or os.path.relpath(DEV, ROOT)}
    rep["all"] = analyse(cur, arms)
    near = [c for c, w in zip(cur["conditions"], cur["wt_growth"]) if 1e-3 <= w < 2e-3]
    rep["near_threshold_conditions"] = {c: float(w) for c, w in zip(cur["conditions"], cur["wt_growth"]) if 1e-3 <= w < 2e-3}
    rep["without_near_threshold_conditions"] = analyse(cur, arms, drop_conditions=set(near))
    name = f"recompute_btheta{'_' + a.tag if a.tag else ''}.json"
    json.dump(rep, open(os.path.join(HERE, name), "w"), indent=1, default=float)
    for k in ("all", "without_near_threshold_conditions"):
        r = rep[k]
        print("==", k)
        print("  own iAH991", round(r["own_iAH991"]["mcc"], 4), r["own_iAH991"]["n_genes"], r["own_iAH991"]["n_conditions"],
              "| own arms", {x: round(y, 4) for x, y in r["own_arms"].items()}, "| coverage", r["coverage"])
        print("  four-way", r["four_way"]["union_genes"], r["four_way"]["conditions_all_grow"], r["four_way"]["n_observations"],
              {x: round(y, 4) for x, y in r["four_way"]["mcc"].items()})
        for n in ("B0", "U'", "M"):
            c, u = r["pairwise"][n]["common"], r["pairwise"][n]["union"]
            print(f"  {n}: common {c['n_genes']}x{c['n_conditions']} arm {c['mcc_arm']:.4f} cur {c['mcc_curated']:.4f} d={c['curated_minus_arm']:+.4f} "
                  f"ci0={[round(x, 4) for x in c['ci95_by_seed'][0]]} | union {u['n_genes']} d={u['curated_minus_arm']:+.4f} ci0={[round(x, 4) for x in u['ci95_by_seed'][0]]}")
        print("  overlap", r["error_overlap"])
        print("  curated-only", r["curated_only"])
    print("near-threshold:", rep["near_threshold_conditions"])


if __name__ == "__main__":
    main()
