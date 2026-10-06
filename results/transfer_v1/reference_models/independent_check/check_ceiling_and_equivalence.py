#!/usr/bin/env python3
"""Context for two wording claims in Paper 1 (no gembench/scripts imports, no LP solve).

(1) "on this benchmark a realistic target for a draft is nearer 0.6 than 1" / "in practice the scale does not run
    to 1": how well do two replicate screens of the same condition agree with each other under the same call rule
    (fitness <= -2)? Computed on the reference models' scored genes and conditions from the raw Fitness Browser
    tables (carbon-source experiments, condition_2 empty or DMSO, grouped by compound x medium). Single-replicate
    agreement is a lower bound on the reliability of the replicate-averaged calls that the benchmark scores.
(2) "as good as" / "matches": paired gene-bootstrap 90% and 95% intervals for curated - arm on the common genes, the
    smallest symmetric margin d for which the 90% interval lies inside (-d, d) (a two one-sided tests reading), and,
    for contrast, the within-P. putida gains U' - B0 and M - B0 on the same genes.
Writes check_ceiling_and_equivalence.json next to this file.
"""
from __future__ import annotations

import itertools
import json
import os

import numpy as np
import pandas as pd

ROOT = "/home/claude/mma"
RES = os.path.join(ROOT, "results", "transfer_v1")
OUT = os.path.join(RES, "reference_models", "independent_check")
FT, GT = -2.0, 1e-3
N_BOOT, SEED = 2000, 7


def mcc_bool(pred, obs):
    tp = float(np.sum(pred & obs)); tn = float(np.sum(~pred & ~obs))
    fp = float(np.sum(pred & ~obs)); fn = float(np.sum(~pred & obs))
    den = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    return (tp * tn - fp * fn) / np.sqrt(den) if den > 0 else 0.0


def load(d):
    with np.load(os.path.join(d, "matrices.npz"), allow_pickle=False) as z:
        return {k: z[k].copy() for k in z.files}


# ------------------------------------------------------------------------------------------- (1) replicate ceiling
def replicate_agreement(org, genes, conds):
    od = os.path.join(ROOT, "data", "fitness_browser", org)
    exp = pd.read_csv(os.path.join(od, "experiments.tsv"), sep="\t", dtype=str, keep_default_na=False)
    fit = pd.read_csv(os.path.join(od, "fit_logratios.tsv"), sep="\t", dtype=str, keep_default_na=False)
    col_of = {c.split(" ", 1)[0]: c for c in fit.columns if c not in ("orgId", "locusId", "sysName", "geneName", "desc")}
    fit.index = fit["sysName"].where(fit["sysName"] != "", fit["locusId"])
    sel = exp[(exp["expGroup"] == "carbon source") & (exp["condition_2"].isin(["", "Dimethyl Sulfoxide"]))]
    groups, set_of = {}, {}
    for _, r in sel.iterrows():
        k = f"{r['condition_1']} | {r['media']}"
        if k in conds and r["expName"] in col_of:
            groups.setdefault(k, []).append(r["expName"])
            set_of[r["expName"]] = r["setName"]
    g = [x for x in genes if x in fit.index]
    a_all, b_all, n_groups = [], [], 0
    first_pair_a, first_pair_b = [], []
    by_set = {"same_set": ([], []), "different_set": ([], [])}
    n_pairs = {"same_set": 0, "different_set": 0}
    for k, exps in groups.items():
        if len(exps) < 2:
            continue
        n_groups += 1
        vals = {e: fit.loc[g, col_of[e]].replace("", np.nan).astype(float).to_numpy() for e in exps}
        for i, (e1, e2) in enumerate(itertools.combinations(exps, 2)):
            x, y = vals[e1], vals[e2]
            ok = np.isfinite(x) & np.isfinite(y)
            a_all.append(x[ok] <= FT); b_all.append(y[ok] <= FT)
            key = "same_set" if set_of[e1] == set_of[e2] else "different_set"
            by_set[key][0].append(x[ok] <= FT); by_set[key][1].append(y[ok] <= FT)
            n_pairs[key] += 1
            if i == 0:
                first_pair_a.append(x[ok] <= FT); first_pair_b.append(y[ok] <= FT)
    a, b = np.concatenate(a_all), np.concatenate(b_all)
    fa, fb = np.concatenate(first_pair_a), np.concatenate(first_pair_b)
    split = {k: {"pairs": n_pairs[k], "mcc": (mcc_bool(np.concatenate(v[0]), np.concatenate(v[1])) if v[0] else None),
                 "cells": int(sum(len(x) for x in v[0]))} for k, v in by_set.items()}
    return {"conditions_with_replicates": n_groups, "conditions_scored": len(conds),
            "mcc_all_replicate_pairs": mcc_bool(a, b), "cells_all_pairs": int(len(a)),
            "mcc_first_pair_per_condition": mcc_bool(fa, fb), "cells_first_pair": int(len(fa)),
            "by_set": split, "replicate_counts": sorted(len(v) for v in groups.values())}


# --------------------------------------------------------------------------------------- (2) equivalence reading
def common(a, b):
    ga = {str(x): i for i, x in enumerate(a["browser_genes"])}
    gb = {str(x): i for i, x in enumerate(b["browser_genes"])}
    ca = {str(x): j for j, x in enumerate(a["conditions"])}
    cb = {str(x): j for j, x in enumerate(b["conditions"])}
    genes = sorted(set(ga) & set(gb))
    conds = [c for c in ca if c in cb and a["wt_growth"][ca[c]] >= GT and b["wt_growth"][cb[c]] >= GT]
    ia, ib = [ga[x] for x in genes], [gb[x] for x in genes]
    ja, jb = [ca[c] for c in conds], [cb[c] for c in conds]
    fa = a["fitness"][np.ix_(ia, ja)]
    assert np.array_equal(fa, b["fitness"][np.ix_(ib, jb)], equal_nan=True)
    return a["sim_growth"][np.ix_(ia, ja)], b["sim_growth"][np.ix_(ib, jb)], fa


def mcc_mat(s, f):
    ok = np.isfinite(s) & np.isfinite(f)
    return mcc_bool((s < GT)[ok], (f <= FT)[ok])


def paired(sa, sb, f):
    rng = np.random.default_rng(SEED)
    n = sa.shape[0]
    d = np.array([mcc_mat(sb[ix], f[ix]) - mcc_mat(sa[ix], f[ix]) for ix in (rng.integers(0, n, n) for _ in range(N_BOOT))])
    pt = mcc_mat(sb, f) - mcc_mat(sa, f)
    ci90, ci95 = np.percentile(d, [5, 95]).tolist(), np.percentile(d, [2.5, 97.5]).tolist()
    return {"point": pt, "ci90": ci90, "ci95": ci95, "smallest_equivalence_margin_90": float(max(abs(ci90[0]), abs(ci90[1]))),
            "share_of_resamples_B_better": float(np.mean(d > 0))}


def main():
    arms = {n: load(os.path.join(RES, "development", "Putida", a)) for n, a in (("B0", "B0"), ("U'", "UNQ"), ("M", "M"))}
    ijn = load(os.path.join(RES, "reference_models", "Putida", "iJN1463"))
    iml = load(os.path.join(RES, "reference_models", "Keio", "iML1515"))
    out = {"n_boot": N_BOOT, "seed": SEED}

    # (1) replicate agreement on the scored genes and growing conditions
    keio_conds = {str(c) for c, w in zip(iml["conditions"], iml["wt_growth"]) if w >= GT}
    out["replicate_agreement_Keio_iML1515_scope"] = replicate_agreement("Keio", [str(x) for x in iml["browser_genes"]], keio_conds)
    put_conds = {str(c) for c, w in zip(arms["B0"]["conditions"], arms["B0"]["wt_growth"]) if w >= GT}
    common_genes = sorted({str(x) for x in arms["B0"]["browser_genes"]} & {str(x) for x in ijn["browser_genes"]})
    out["replicate_agreement_Putida_764_genes_28_conditions"] = replicate_agreement("Putida", common_genes, put_conds)

    # (2) equivalence readings on the common genes
    eq = {}
    for n in ("B0", "U'", "M"):
        sa, sb, f = common(arms[n], ijn)
        eq[f"iJN1463_minus_{n}"] = paired(sa, sb, f)
    for x, y in (("B0", "U'"), ("B0", "M")):
        # within the drafts, restricted to genes shared with the curated model (764), conditions where all grow
        sx, sc, f = common(arms[x], ijn)
        sy, _, _ = common(arms[y], ijn) if y != "M" else (None, None, None)
        if y == "M":
            # M has one extra common gene; align B0 and M on B0's common genes
            gm = {str(g): i for i, g in enumerate(arms["M"]["browser_genes"])}
            gb = sorted({str(g) for g in arms["B0"]["browser_genes"]} & {str(g) for g in ijn["browser_genes"]})
            cm = {str(c): j for j, c in enumerate(arms["M"]["conditions"])}
            conds = [str(c) for c, w in zip(arms["B0"]["conditions"], arms["B0"]["wt_growth"]) if w >= GT]
            sy = arms["M"]["sim_growth"][np.ix_([gm[g] for g in gb], [cm[c] for c in conds])]
        eq[f"{y}_minus_{x}_on_764_shared_genes"] = paired(sx, sy, f)
    out["equivalence"] = eq
    with open(os.path.join(OUT, "check_ceiling_and_equivalence.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
        fh.write("\n")
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
