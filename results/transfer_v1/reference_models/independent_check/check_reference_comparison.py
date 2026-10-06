#!/usr/bin/env python3
"""Independent check of the hand-curated reference comparison in Paper 1 (docs/paper/paper1-draft.md).

Written from the task's definitions only. It reads the saved run matrices (matrices.npz) and cards (card.json)
with numpy/json/pandas; it does not import gembench or anything under scripts/, and it runs no LP solve.

Definitions
  * growth threshold and fitness threshold are read from each card (protocol.params);
  * a condition counts for a model if its wild type grows (wt_growth >= growth threshold, finite);
  * predicted important: knockout growth < growth threshold; observed important: fitness <= fitness threshold
    (the task's definition; the code under test uses fitness < threshold, so both are computed and the number of
    cells with fitness exactly equal to the threshold is reported);
  * only cells with finite fitness (and finite simulation) count;
  * common-gene comparison: genes in both models, conditions where both wild types grow;
  * union comparison: genes of either model; a gene absent from a model gets that model's wild-type growth
    (predicted "no effect"); conditions where both (or all) wild types grow;
  * gene bootstrap: resample gene rows with replacement and recompute the MCC difference (row resampling of the
    actual arrays, not per-gene count weights).

Writes check_reference_comparison.json next to this file.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = "/home/claude/mma"
RES = os.path.join(ROOT, "results", "transfer_v1")
REF = os.path.join(RES, "reference_models")
OUT = os.path.join(REF, "independent_check")

N_BOOT = 1000
SEEDS = [20261006, 1, 2, 3, 4]          # first seed is the reported one; the others show seed sensitivity


# ----------------------------------------------------------------------------------------------------------------
# loading

def load(d):
    with np.load(os.path.join(d, "matrices.npz"), allow_pickle=False) as z:
        r = {k: z[k].copy() for k in z.files}
    card = json.load(open(os.path.join(d, "card.json")))
    p = card["protocol"]["params"]
    r["gt"] = float(p["growth_threshold"])
    r["ft"] = float(p["fitness_threshold"])
    r["card"] = card
    r["dir"] = d
    genes = [str(g) for g in r["browser_genes"]]
    conds = [str(c) for c in r["conditions"]]
    assert len(set(genes)) == len(genes), f"duplicate genes in {d}"
    assert len(set(conds)) == len(conds), f"duplicate conditions in {d}"
    assert r["sim_growth"].shape == (len(genes), len(conds)) == r["fitness"].shape
    assert r["wt_growth"].shape == (len(conds),)
    r["genes"], r["conds"] = genes, conds
    r["gi"] = {g: i for i, g in enumerate(genes)}
    r["ci"] = {c: j for j, c in enumerate(conds)}
    return r


def grows(r):
    w = r["wt_growth"]
    return np.isfinite(w) & (w >= r["gt"])


# ----------------------------------------------------------------------------------------------------------------
# metric

def confusion(sim, fit, gt, ft, le=True):
    """Counts with 'important' as the positive class. le=True: observed important if fit <= ft (task definition);
    le=False: fit < ft (definition in the code under test)."""
    ok = np.isfinite(sim) & np.isfinite(fit)
    pred = sim < gt
    obs = (fit <= ft) if le else (fit < ft)
    tp = int(np.sum(ok & pred & obs))
    tn = int(np.sum(ok & ~pred & ~obs))
    fp = int(np.sum(ok & pred & ~obs))
    fn = int(np.sum(ok & ~pred & obs))
    return {"tp": tp, "tn": tn, "fp": fp, "fn": fn, "n": tp + tn + fp + fn}


def mcc_c(c):
    tp, tn, fp, fn = (float(c[k]) for k in ("tp", "tn", "fp", "fn"))
    den = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    if tp + tn + fp + fn == 0:
        return float("nan")
    return (tp * tn - fp * fn) / np.sqrt(den) if den > 0 else 0.0


def mcc(sim, fit, gt, ft, le=True):
    return mcc_c(confusion(sim, fit, gt, ft, le))


# ----------------------------------------------------------------------------------------------------------------
# alignments

def own_scope(r, cond_mask=None):
    m = grows(r) if cond_mask is None else (grows(r) & cond_mask)
    return r["sim_growth"][:, m], r["fitness"][:, m]


def merged_fitness(runs, genes, conds):
    """genes x conds fitness from whichever run has the gene; raises if two runs disagree on a finite value."""
    fit = np.full((len(genes), len(conds)), np.nan)
    n_compared = 0
    for r in runs:
        rows = [r["gi"].get(g) for g in genes]
        cols = [r["ci"][c] for c in conds]
        for i, ri in enumerate(rows):
            if ri is None:
                continue
            f = r["fitness"][ri, cols]
            both = np.isfinite(fit[i]) & np.isfinite(f)
            n_compared += int(both.sum())
            if np.any(both & (fit[i] != f)):
                raise ValueError(f"fitness differs between runs for gene {genes[i]}")
            # a value finite in one run and NaN in another would also be a data mismatch
            if np.any(np.isfinite(fit[i]) != np.isfinite(f)) and np.any(np.isfinite(fit[i])):
                raise ValueError(f"finite/NaN mismatch between runs for gene {genes[i]}")
            fit[i] = np.where(np.isfinite(f), f, fit[i])
    return fit, n_compared


def conds_all_grow(runs):
    common = [c for c in runs[0]["conds"] if all(c in r["ci"] for r in runs)]
    return [c for c in common if all(np.isfinite(r["wt_growth"][r["ci"][c]]) and r["wt_growth"][r["ci"][c]] >= r["gt"]
                                     for r in runs)]


def sim_on(r, genes, conds):
    """Knockout growth on the given genes/conditions; absent genes get the wild-type growth (no effect)."""
    cols = [r["ci"][c] for c in conds]
    wt = r["wt_growth"][cols]
    s = np.tile(wt, (len(genes), 1)).astype(float)
    present = np.zeros(len(genes), bool)
    for i, g in enumerate(genes):
        ri = r["gi"].get(g)
        if ri is not None:
            s[i] = r["sim_growth"][ri, cols]
            present[i] = True
    return s, present


def common_view(runs):
    genes = sorted(set.intersection(*[set(r["genes"]) for r in runs]))
    conds = conds_all_grow(runs)
    fit, ncomp = merged_fitness(runs, genes, conds)
    sims = [sim_on(r, genes, conds)[0] for r in runs]
    return genes, conds, sims, fit, ncomp


def union_view(runs):
    genes = sorted(set().union(*[set(r["genes"]) for r in runs]))
    conds = conds_all_grow(runs)
    fit, ncomp = merged_fitness(runs, genes, conds)
    sims, pres = zip(*[sim_on(r, genes, conds) for r in runs])
    return genes, conds, list(sims), fit, ncomp, list(pres)


# ----------------------------------------------------------------------------------------------------------------
# bootstrap

def boot_diff(sa, sb, fit, gt, ft, n_boot, seed, le=True):
    """Percentile interval of MCC(b) - MCC(a), resampling gene rows with replacement."""
    rng = np.random.default_rng(seed)
    n = sa.shape[0]
    d = np.empty(n_boot)
    for k in range(n_boot):
        ix = rng.integers(0, n, n)
        d[k] = mcc(sb[ix], fit[ix], gt, ft, le) - mcc(sa[ix], fit[ix], gt, ft, le)
    d = d[np.isfinite(d)]
    return [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))], int(len(d))


def boot_single(s, fit, gt, ft, n_boot, seed, le=True):
    rng = np.random.default_rng(seed)
    n = s.shape[0]
    v = np.array([mcc(s[ix], fit[ix], gt, ft, le) for ix in (rng.integers(0, n, n) for _ in range(n_boot))])
    v = v[np.isfinite(v)]
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


# ----------------------------------------------------------------------------------------------------------------

def errors(sim, fit, gt, ft, le=True):
    ok = np.isfinite(sim) & np.isfinite(fit)
    pred = sim < gt
    obs = (fit <= ft) if le else (fit < ft)
    return ok & (pred != obs), ok & pred & ~obs, ok & ~pred & obs, ok


def main():
    out = {"definitions": __doc__.strip().splitlines()[0], "n_boot": N_BOOT, "seeds": SEEDS}

    putida = {n: load(os.path.join(RES, "development", "Putida", a)) for n, a in (("B0", "B0"), ("U'", "UNQ"), ("M", "M"))}
    ijn = load(os.path.join(REF, "Putida", "iJN1463"))
    iml = load(os.path.join(REF, "Keio", "iML1515"))
    allp = dict(putida, iJN1463=ijn)

    thr = {n: (r["gt"], r["ft"]) for n, r in list(allp.items()) + [("iML1515", iml)]}
    assert len(set(thr.values())) == 1, thr
    gt, ft = thr["B0"]
    out["thresholds"] = {"growth": gt, "fitness": ft}

    # boundary cells: fitness exactly at the threshold
    out["cells_with_fitness_exactly_at_threshold"] = {
        n: int(np.sum(own_scope(r)[1] == ft)) for n, r in list(allp.items()) + [("iML1515", iml)]}

    # ---------------- 1. iML1515 ----------------
    s, f = own_scope(iml)
    c_le, c_lt = confusion(s, f, gt, ft, True), confusion(s, f, gt, ft, False)
    rows_any = np.isfinite(f).any(axis=1)
    out["check1_iML1515"] = {
        "n_conditions_mapped": len(iml["conds"]), "n_conditions_wt_grows": int(grows(iml).sum()),
        "n_genes": int(rows_any.sum()), "n_cells": c_le["n"],
        "mcc_le": mcc_c(c_le), "mcc_lt": mcc_c(c_lt), "confusion_le": c_le,
        "mcc_ci95_gene_bootstrap": boot_single(s, f, gt, ft, N_BOOT, SEEDS[0]),
        "card_mcc": iml["card"]["results"]["gene_level_conditions_where_wt_grows"]["mcc"]["point"],
        "card_confusion_(positive=no defect)": iml["card"]["results"]["gene_level_conditions_where_wt_grows"]["confusion"],
        "wt_growth_min": float(np.min(iml["wt_growth"]))}

    # ---------------- own-scope MCCs of all P. putida models (context) ----------------
    own = {}
    for n, r in allp.items():
        s, f = own_scope(r)
        own[n] = {"n_genes": int(np.isfinite(f).any(axis=1).sum()), "n_conditions": int(grows(r).sum()),
                  "mcc_le": mcc(s, f, gt, ft, True), "mcc_lt": mcc(s, f, gt, ft, False),
                  "card_mcc": r["card"]["results"]["gene_level_conditions_where_wt_grows"]["mcc"]["point"]}
    out["own_scope_mcc_putida"] = own

    # ---------------- 4. coverage ----------------
    cov = {n: {"mapped": len(r["conds"]), "grows": int(grows(r).sum()),
               "card": r["card"]["results"]["condition_level"]} for n, r in allp.items()}
    assert all(r["conds"] == ijn["conds"] for r in allp.values()), "condition lists differ"
    gsets = {n: {c for c, g in zip(r["conds"], grows(r)) if g} for n, r in allp.items()}
    draft_sets_equal = gsets["B0"] == gsets["U'"] == gsets["M"]
    subset = gsets["B0"] <= gsets["iJN1463"]
    ctab = {}
    for n, r in allp.items():
        t = pd.read_csv(os.path.join(r["dir"], "conditions.tsv"), sep="\t", keep_default_na=False)
        ctab[n] = t
    cond_rows = []
    for j, c in enumerate(ijn["conds"]):
        row = {"condition": c}
        for n in allp:
            row[f"grows_{n}"] = c in gsets[n]
            t = ctab[n]
            k = t.index[(t["condition"] + " | " + t["media"]) == c]
            row[f"absent_exchange_{n}"] = bool(len(k) and str(t.loc[k[0], "missing_exchanges"]).strip() != "")
        cond_rows.append(row)
    pd.DataFrame(cond_rows).to_csv(os.path.join(OUT, "putida_condition_coverage.tsv"), sep="\t", index=False)
    out["check4_coverage"] = {
        "per_model": cov, "draft_growth_sets_identical_B0_U_M": draft_sets_equal,
        "draft_28_subset_of_curated_36": subset,
        "curated_only_conditions": sorted(gsets["iJN1463"] - gsets["B0"]),
        "draft_only_conditions": sorted(gsets["B0"] - gsets["iJN1463"]),
        "neither_grows": sorted(set(ijn["conds"]) - gsets["iJN1463"] - gsets["B0"]),
        "draft_nongrowth_with_absent_exchange": int(sum(1 for r in cond_rows if not r["grows_B0"] and r["absent_exchange_B0"])),
        "curated_nongrowth_with_absent_exchange": int(sum(1 for r in cond_rows if not r["grows_iJN1463"] and r["absent_exchange_iJN1463"])),
        "conditions_total_in_card_counts": {n: r["card"]["results"]["counts"]["conditions_total"] for n, r in allp.items()}}

    # ---------------- 2. common-gene comparison ----------------
    common = {}
    for n in ("B0", "U'", "M"):
        genes, conds, (sa, sb), fit, ncomp = common_view([putida[n], ijn])
        rec = {"n_common_genes": len(genes), "n_conditions_both_grow": len(conds),
               "n_cells": int(np.isfinite(fit).sum()), "fitness_cells_cross_checked": ncomp,
               "mcc_arm_le": mcc(sa, fit, gt, ft, True), "mcc_curated_le": mcc(sb, fit, gt, ft, True),
               "mcc_arm_lt": mcc(sa, fit, gt, ft, False), "mcc_curated_lt": mcc(sb, fit, gt, ft, False)}
        rec["curated_minus_arm"] = rec["mcc_curated_le"] - rec["mcc_arm_le"]
        rec["ci95_by_seed"] = {str(sd): boot_diff(sa, sb, fit, gt, ft, N_BOOT, sd)[0] for sd in SEEDS}
        rec["ci95"] = rec["ci95_by_seed"][str(SEEDS[0])]
        common[n] = rec
    # B0 vs M directly on the same common set (for "brought ... to the level")
    out["check2_common"] = common

    # all four on one common set (genes present in all four models, conditions where all grow)
    genes4, conds4, sims4, fit4, _ = common_view([putida["B0"], putida["U'"], putida["M"], ijn])
    names4 = ["B0", "U'", "M", "iJN1463"]
    out["common_all_four"] = {"n_genes": len(genes4), "n_conditions": len(conds4),
                              "mcc": {n: mcc(s, fit4, gt, ft) for n, s in zip(names4, sims4)}}

    # ---------------- 3. union ----------------
    genesu, condsu, simsu, fitu, ncu, presu = union_view([putida["B0"], putida["U'"], putida["M"], ijn])
    four = {"union_genes": len(genesu), "conditions_all_grow": len(condsu), "n_cells": int(np.isfinite(fitu).sum()),
            "mcc": {n: mcc(s, fitu, gt, ft) for n, s in zip(names4, simsu)},
            "mcc_lt": {n: mcc(s, fitu, gt, ft, False) for n, s in zip(names4, simsu)}}
    pair_union = {}
    for n in ("B0", "U'", "M"):
        g, c, (sa, sb), fit, _, _ = union_view([putida[n], ijn])
        ci, _ = boot_diff(sa, sb, fit, gt, ft, N_BOOT, SEEDS[0])
        pair_union[n] = {"union_genes": len(g), "conditions": len(c), "mcc_arm": mcc(sa, fit, gt, ft),
                         "mcc_curated": mcc(sb, fit, gt, ft), "curated_minus_arm": mcc(sb, fit, gt, ft) - mcc(sa, fit, gt, ft),
                         "ci95": ci}
    out["check3_union"] = {"four_way": four, "pairwise": pair_union}

    # decomposition of the union: genes only in one side (U' vs curated, 28 conditions)
    g, c, (sa, sb), fit, _, (pa, pb) = union_view([putida["U'"], ijn])
    dec = {}
    for label, mask, s in (("curated_only_genes_curated_calls", pb & ~pa, sb), ("draft_only_genes_draft_calls", pa & ~pb, sa),
                           ("shared_genes_curated_calls", pa & pb, sb), ("shared_genes_draft_calls", pa & pb, sa)):
        cc = confusion(s[mask], fit[mask], gt, ft)
        dec[label] = {"n_genes": int(mask.sum()), "predicted_important": cc["tp"] + cc["fp"], "confirmed": cc["tp"],
                      "observed_important_cells": cc["tp"] + cc["fn"], "missed": cc["fn"]}
    out["union_decomposition_Uprime_vs_iJN1463"] = dec

    # ---------------- curated model on its 8 extra conditions and on the 28 shared ones ----------------
    shared_mask = np.array([c in gsets["B0"] for c in ijn["conds"]])
    s28, f28 = own_scope(ijn, shared_mask)
    s8, f8 = own_scope(ijn, ~shared_mask)
    out["curated_by_condition_set"] = {
        "own_genes_28_shared_conditions": {"mcc": mcc(s28, f28, gt, ft), "confusion": confusion(s28, f28, gt, ft)},
        "own_genes_8_extra_conditions": {"n_conditions": int((grows(ijn) & ~shared_mask).sum()), "mcc": mcc(s8, f8, gt, ft),
                                         "confusion": confusion(s8, f8, gt, ft),
                                         "ci95": boot_single(s8, f8, gt, ft, N_BOOT, SEEDS[0])}}

    # ---------------- error overlap on the common observations ----------------
    ov = {}
    for n in ("B0", "U'", "M"):
        genes, conds, (sa, sb), fit, _ = common_view([putida[n], ijn])
        ea, fpa, fna, ok = errors(sa, fit, gt, ft)
        eb, fpb, fnb, _ = errors(sb, fit, gt, ft)
        N = int(ok.sum())
        both = int((ea & eb).sum())
        exp_indep = float(ea.sum()) * float(eb.sum()) / N
        obs_imp = ok & (fit <= ft)
        ov[n] = {"n_cells": N, "errors_arm": int(ea.sum()), "errors_curated": int(eb.sum()), "errors_both": both,
                 "share_of_arm_errors_also_curated_errors": both / max(1, int(ea.sum())),
                 "share_of_curated_errors_also_arm_errors": both / max(1, int(eb.sum())),
                 "jaccard": both / max(1, int((ea | eb).sum())),
                 "expected_both_if_independent": exp_indep,
                 "arm_false_important": int(fpa.sum()), "arm_missed_important": int(fna.sum()),
                 "curated_false_important": int(fpb.sum()), "curated_missed_important": int(fnb.sum()),
                 "missed_by_both": int((fna & fnb).sum()), "false_important_by_both": int((fpa & fpb).sum()),
                 "observed_important_cells": int(obs_imp.sum()),
                 "genes_with_any_error_arm": int(ea.any(axis=1).sum()), "genes_with_any_error_curated": int(eb.any(axis=1).sum()),
                 "genes_with_errors_in_both": int((ea.any(axis=1) & eb.any(axis=1)).sum())}
    out["error_overlap_common"] = ov

    # ---------------- B0 -> U' -> M changes on common genes relative to curated (where do they move?) -------------
    genes, conds, (s_b0, s_u, s_m, s_c), fit, _ = common_view([putida["B0"], putida["U'"], putida["M"], ijn])
    ok = np.isfinite(fit)
    pb0, pu, pm, pc = (s < gt for s in (s_b0, s_u, s_m, s_c))
    chg = ok & (pb0 != pu)
    out["B0_to_Uprime_changes_on_four_way_common"] = {
        "n_cells_changed": int(chg.sum()), "now_agree_with_curated": int((chg & (pu == pc)).sum()),
        "now_agree_with_data": int((chg & (pu == (fit <= ft))).sum())}

    with open(os.path.join(OUT, "check_reference_comparison.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
        fh.write("\n")
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    sys.exit(main())
