"""Check 5: recompute the MR1/iSO783 entry of comparison.json from the saved matrices with independent code.

Definitions (as stated in the plan and the transfer method; written here from the definitions, not from the code):
  * observation = gene x condition cell with finite fitness and finite knockout growth;
  * observed important = replicate-averaged fitness <= -2 (no value equals -2 exactly, so < and <= agree);
  * predicted important = knockout growth < 1e-3;
  * own MCC: the model's genes, conditions where its wild type grows (>= 1e-3);
  * common genes: genes in both runs, conditions where both wild types grow;
  * union: genes in any of the runs, conditions where every wild type grows; a gene absent from a run is predicted
    to have no effect in that run (its knockout growth = that run's wild-type growth);
  * error overlap: cells of the common-gene table where U' and the curated model are wrong;
  * curated-only genes: genes of the curated run in none of B0, U', M; conditions where the curated model and U'
    both grow.
Bootstrap: genes resampled with replacement (1,000 draws), several seeds, percentile 95% intervals.

Usage: python -I recompute_comparison.py [--curated <npz>] [--tag <name>]
Writes recompute_comparison[_<tag>].json.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import OUT, ROOT, dump  # noqa: E402

ST, FT = 1e-3, -2.0
DEV = os.path.join(ROOT, "results", "transfer_v1", "development", "MR1")
REF = os.path.join(ROOT, "results", "transfer_v1", "reference_models", "MR1", "iSO783", "matrices.npz")


def load(path):
    z = np.load(path, allow_pickle=False)
    r = {k: z[k] for k in ("sim_growth", "wt_growth", "fitness", "browser_genes", "conditions")}
    r["browser_genes"] = [str(x) for x in r["browser_genes"]]
    r["conditions"] = [str(x) for x in r["conditions"]]
    assert len(set(r["browser_genes"])) == len(r["browser_genes"])
    return r


def counts(sim, fit):
    ok = np.isfinite(sim) & np.isfinite(fit)
    pi = sim < ST
    oi = fit <= FT
    tp = (ok & pi & oi).sum(axis=-1)
    tn = (ok & ~pi & ~oi).sum(axis=-1)
    fp = (ok & pi & ~oi).sum(axis=-1)
    fn = (ok & ~pi & oi).sum(axis=-1)
    return np.stack([tp, tn, fp, fn], axis=-1).astype(float)   # per gene (rows)


def mcc(c):
    tp, tn, fp, fn = c
    den = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    return float((tp * tn - fp * fn) / np.sqrt(den)) if den > 0 else 0.0


def grows(r):
    return {c for c, w in zip(r["conditions"], r["wt_growth"]) if np.isfinite(w) and w >= ST}


def own(r):
    cols = [j for j, c in enumerate(r["conditions"]) if c in grows(r)]
    c = counts(r["sim_growth"][:, cols], r["fitness"][:, cols]).sum(0)
    return {"n_genes": len(r["browser_genes"]), "n_conditions": len(cols), "n_cells": int(c.sum()),
            "confusion_important_positive": dict(zip(["tp", "tn", "fp", "fn"], map(int, c))), "mcc": mcc(c)}


def table(runs, genes, conds, absent_as_wt):
    """Stack runs on given genes x conditions. Returns sims (list) and the common fitness matrix."""
    fit = np.full((len(genes), len(conds)), np.nan)
    sims = []
    for r in runs:
        gi = {g: i for i, g in enumerate(r["browser_genes"])}
        ci = {c: j for j, c in enumerate(r["conditions"])}
        cols = [ci[c] for c in conds]
        wt = np.array([r["wt_growth"][j] for j in cols])
        sim = np.tile(wt, (len(genes), 1)) if absent_as_wt else np.full((len(genes), len(conds)), np.nan)
        for i, g in enumerate(genes):
            if g in gi:
                sim[i] = r["sim_growth"][gi[g], cols]
                f = r["fitness"][gi[g], cols]
                if np.any(np.isfinite(fit[i]) & np.isfinite(f) & (fit[i] != f)):
                    raise ValueError(f"fitness differs between runs for {g}")
                fit[i] = np.where(np.isfinite(f), f, fit[i])
        sims.append(sim)
    return sims, fit


def boot(ca, cb, seeds=(0, 1, 2, 3, 4), n=1000):
    out = []
    for s in seeds:
        rng = np.random.default_rng(s)
        d = []
        for _ in range(n):
            w = np.bincount(rng.integers(0, len(ca), len(ca)), minlength=len(ca))
            d.append(mcc(w @ cb) - mcc(w @ ca))
        out.append([float(x) for x in np.percentile(d, [2.5, 97.5])])
    return out


def pair(a, b, mode):
    if mode == "common":
        genes = sorted(set(a["browser_genes"]) & set(b["browser_genes"]))
    else:
        genes = sorted(set(a["browser_genes"]) | set(b["browser_genes"]))
    conds = sorted(grows(a) & grows(b))
    (sa, sb), fit = table([a, b], genes, conds, absent_as_wt=(mode == "union"))
    if mode == "common":
        ok = np.isfinite(sa) & np.isfinite(sb) & np.isfinite(fit)
        sa, sb, fit = np.where(ok, sa, np.nan), np.where(ok, sb, np.nan), np.where(ok, fit, np.nan)
    ca, cb = counts(sa, fit), counts(sb, fit)
    return {"n_genes": len(genes), "n_conditions": len(conds), "mcc_arm": mcc(ca.sum(0)), "mcc_curated": mcc(cb.sum(0)),
            "curated_minus_arm": mcc(cb.sum(0)) - mcc(ca.sum(0)), "ci95_by_seed": boot(ca, cb)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curated", default=REF)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    cur = load(a.curated)
    arms = {"B0": load(os.path.join(DEV, "B0", "matrices.npz")), "U'": load(os.path.join(DEV, "UNQ", "matrices.npz")),
            "M": load(os.path.join(DEV, "M", "matrices.npz"))}
    rep = {"curated_matrices": os.path.relpath(a.curated, ROOT)}
    rep["exact_minus_two_fitness_values"] = int(sum((r["fitness"] == -2).sum() for r in [cur, *arms.values()]))
    rep["own"] = {"iSO783": own(cur), **{k: own(v) for k, v in arms.items()}}
    rep["coverage"] = {k: {"mapped": len(r["conditions"]), "grows": len(grows(r)), "grows_in": sorted(grows(r))}
                       for k, r in {"iSO783": cur, **arms}.items()}
    # four-way union
    runs = [arms["B0"], arms["U'"], arms["M"], cur]
    names = ["B0", "U'", "M", "iSO783"]
    genes = sorted(set().union(*[set(r["browser_genes"]) for r in runs]))
    conds = sorted(set.intersection(*[grows(r) for r in runs]))
    sims, fit = table(runs, genes, conds, absent_as_wt=True)
    m4 = {n: mcc(counts(s, fit).sum(0)) for n, s in zip(names, sims)}
    gap = m4["iSO783"] - m4["B0"]
    rep["four_way"] = {"union_genes": len(genes), "conditions_all_grow": len(conds), "conditions": conds,
                       "n_observations": int(np.isfinite(fit).sum()), "mcc": m4,
                       "share_of_gap_closed": {"U'": (m4["U'"] - m4["B0"]) / gap, "M": (m4["M"] - m4["B0"]) / gap} if gap > 0 else None}
    rep["pairwise_curated_minus_arm"] = {n: {"common": pair(arms[n], cur, "common"), "union": pair(arms[n], cur, "union")}
                                         for n in ("B0", "U'", "M")}
    # error overlap, common genes, U' vs curated
    u = arms["U'"]
    cg = sorted(set(u["browser_genes"]) & set(cur["browser_genes"]))
    cc = sorted(grows(u) & grows(cur))
    (su, sc), fc = table([u, cur], cg, cc, absent_as_wt=False)
    ok = np.isfinite(su) & np.isfinite(sc) & np.isfinite(fc)
    obs = fc <= FT
    wu = ok & ((su < ST) != obs)
    wc = ok & ((sc < ST) != obs)
    n = int(ok.sum())
    rep["error_overlap_common_genes_Uprime_vs_curated"] = {
        "n_genes": len(cg), "n_conditions": len(cc), "n_observations": n, "wrong_Uprime": int(wu.sum()),
        "wrong_curated": int(wc.sum()), "wrong_both": int((wu & wc).sum()),
        "expected_both_if_independent": float(wu.sum() * wc.sum() / n),
        "wrong_both_by_type": {"both_false_important": int((wu & wc & (su < ST) & (sc < ST)).sum()),
                               "both_false_no_effect": int((wu & wc & (su >= ST) & (sc >= ST)).sum())},
        "share_of_Uprime_errors_shared": float((wu & wc).sum() / wu.sum())}
    # curated-only genes
    draft = set().union(*[set(r["browser_genes"]) for r in arms.values()])
    only = [g for g in cur["browser_genes"] if g not in draft]
    gi = {g: i for i, g in enumerate(cur["browser_genes"])}
    cols = [j for j, c in enumerate(cur["conditions"]) if c in grows(cur) and c in set(cc)]
    so = cur["sim_growth"][np.ix_([gi[g] for g in only], cols)]
    fo = cur["fitness"][np.ix_([gi[g] for g in only], cols)]
    oko = np.isfinite(so) & np.isfinite(fo)
    imp = oko & (so < ST)
    rep["curated_only_genes"] = {"n_genes": len(only), "n_conditions": len(cols), "important_calls": int(imp.sum()),
                                 "important_calls_confirmed": int((imp & (fo <= FT)).sum()),
                                 "observed_important_cells": int((oko & (fo <= FT)).sum()),
                                 "genes_with_any_important_call": int(imp.any(axis=1).sum()),
                                 "conditions_used": [cur["conditions"][j] for j in cols]}
    # curated-only genes over all conditions where the curated model grows (for the wording)
    cols_all = [j for j, c in enumerate(cur["conditions"]) if c in grows(cur)]
    so2 = cur["sim_growth"][np.ix_([gi[g] for g in only], cols_all)]
    fo2 = cur["fitness"][np.ix_([gi[g] for g in only], cols_all)]
    imp2 = np.isfinite(so2) & np.isfinite(fo2) & (so2 < ST)
    rep["curated_only_genes_all_curated_conditions"] = {"n_conditions": len(cols_all), "important_calls": int(imp2.sum()),
                                                       "important_calls_confirmed": int((imp2 & (fo2 <= FT)).sum())}
    name = f"recompute_comparison{'_' + a.tag if a.tag else ''}.json"
    dump(rep, name)
    print(json.dumps({"own": {k: round(v["mcc"], 4) for k, v in rep["own"].items()},
                      "four_way": {k: (round(v, 4) if isinstance(v, float) else v) for k, v in rep["four_way"]["mcc"].items()},
                      "union_genes": rep["four_way"]["union_genes"], "conds": rep["four_way"]["conditions_all_grow"],
                      "common": {k: (v["common"]["n_genes"], v["common"]["n_conditions"], round(v["common"]["mcc_arm"], 4),
                                     round(v["common"]["mcc_curated"], 4), round(v["common"]["curated_minus_arm"], 4),
                                     [[round(x, 3) for x in ci] for ci in v["common"]["ci95_by_seed"]])
                                 for k, v in rep["pairwise_curated_minus_arm"].items()},
                      "union": {k: (v["union"]["n_genes"], round(v["union"]["curated_minus_arm"], 4),
                                    [[round(x, 3) for x in ci] for ci in v["union"]["ci95_by_seed"]])
                                for k, v in rep["pairwise_curated_minus_arm"].items()},
                      "overlap": rep["error_overlap_common_genes_Uprime_vs_curated"],
                      "curated_only": rep["curated_only_genes"]}, indent=1, default=float))


if __name__ == "__main__":
    main()
