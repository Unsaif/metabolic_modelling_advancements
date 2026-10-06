"""How close do the corrected drafts get to a hand-curated model? (Paper 1, post hoc and descriptive.)

For each organism with a curated reference scored by scripts/score_reference_model.py, this compares the transfer
study's draft arms (B0 = gap-filled draft, U' = UNQ = B0 + the five automatic rules, M = U' + blind AI curation)
with the curated model, using the study's own definitions:

  * four-way table: all four models on the same observations (union of their genes, conditions where all four wild
    types grow; a gene absent from a model is predicted to have no effect there, as in the transfer study's primary
    analysis), with the share of the draft-to-curated gap that U' and M close;
  * pairwise curated - arm differences on the union and on common genes, with gene-bootstrap 95% intervals
    (scripts/transfer_compare.py);
  * coverage: conditions in which each model grows, out of the conditions mapped (the screen grew the organism in
    every one of them).

Organisms without draft arms (E. coli, the benchmark's best-curated control) get coverage and their own MCC only.

Usage: python scripts/compare_reference_models.py [--out results/transfer_v1/reference_models/comparison.json]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gembench.comparison import load_run  # noqa: E402
import transfer_compare as TC  # noqa: E402

RES = os.path.join(ROOT, "results", "transfer_v1")
REF = os.path.join(RES, "reference_models")
ARMS = (("B0", "B0"), ("U'", "UNQ"), ("M", "M"))


def union_n(runs):
    """Union of all runs' genes, conditions where every wild type grows (generalises transfer_compare.union_pairs)."""
    st = runs[0]["params"]["growth_threshold"]
    genes = sorted(set().union(*[set(r["browser_genes"]) for r in runs]))
    conds = sorted(set.intersection(*[set(r["conditions"]) for r in runs]))
    gi = [{g: i for i, g in enumerate(r["browser_genes"])} for r in runs]
    ci = [{c: i for i, c in enumerate(r["conditions"])} for r in runs]
    grows = [c for c in conds if all(np.isfinite(r["wt_growth"][ix[c]]) and r["wt_growth"][ix[c]] >= st
                                     for r, ix in zip(runs, ci))]
    fit = np.full((len(genes), len(grows)), np.nan)
    sims = []
    for r, g_ix, c_ix in zip(runs, gi, ci):
        cols = [c_ix[c] for c in grows]
        sim = np.tile(np.array([r["wt_growth"][j] for j in cols]), (len(genes), 1))
        for i, g in enumerate(genes):
            if g in g_ix:
                sim[i, :] = r["sim_growth"][g_ix[g], cols]
                f = r["fitness"][g_ix[g], cols]
                known = np.isfinite(fit[i, :])
                if np.any(known & np.isfinite(f) & (fit[i, :] != f)):
                    raise ValueError(f"fitness differs between runs for {g}")
                fit[i, :] = np.where(np.isfinite(f), f, fit[i, :])
        sims.append(sim)
    return sims, fit, genes, grows


def coverage(run):
    st = run["params"]["growth_threshold"]
    w = run["wt_growth"]
    return {"n_conditions_mapped": int(len(w)), "n_conditions_grows": int((np.isfinite(w) & (w >= st)).sum())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(REF, "comparison.json"))
    ap.add_argument("--n-boot", type=int, default=1000)
    args = ap.parse_args()
    report = {"note": "post hoc and descriptive; development organisms and the E. coli control are not held out",
              "notes": ["share_of_gap_closed is exploratory (not in docs/studies/transfer-v1-curated-references-plan.md): "
                        "a ratio of differences whose intervals include zero",
                        "entries whose card says exploratory=true used a medium supplement (protocol deviation)"],
              "organisms": {}}
    for org_dir in sorted(glob.glob(os.path.join(REF, "*"))):
        if not os.path.isdir(org_dir):
            continue
        org = os.path.basename(org_dir)
        for ref_dir in sorted(glob.glob(os.path.join(org_dir, "*"))):
            if not os.path.exists(os.path.join(ref_dir, "card.json")):
                continue
            label = os.path.basename(ref_dir)
            ref = load_run(ref_dir)
            card = json.load(open(os.path.join(ref_dir, "card.json")))
            own = card["results"]["gene_level_conditions_where_wt_grows"]
            entry = {"reference": label, "coverage": {label: coverage(ref)},
                     "own_mcc": {label: own.get("mcc", {}).get("point")},
                     "exploratory": bool(card["protocol"].get("exploratory", False)),
                     "medium_supplement": card["protocol"].get("medium_supplement", {}),
                     "conditions_reference_grows": [c for c, w in zip(ref["conditions"], ref["wt_growth"])
                                                    if np.isfinite(w) and w >= ref["params"]["growth_threshold"]]}
            arm_dirs = {name: os.path.join(RES, "development", org, arm) for name, arm in ARMS}
            have_arms = all(os.path.exists(os.path.join(d, "card.json")) for d in arm_dirs.values())
            if have_arms and coverage(ref)["n_conditions_grows"] == 0:
                for name, d in arm_dirs.items():
                    r = load_run(d)
                    entry["coverage"][name] = coverage(r)
                    c = json.load(open(os.path.join(d, "card.json")))
                    entry["own_mcc"][name] = c["results"]["gene_level_conditions_where_wt_grows"].get("mcc", {}).get("point")
                entry["gene_level"] = "not computed: the reference model grows in none of the mapped conditions"
            elif have_arms:
                arms = {name: load_run(d) for name, d in arm_dirs.items()}
                for name, r in arms.items():
                    entry["coverage"][name] = coverage(r)
                    c = json.load(open(os.path.join(arm_dirs[name], "card.json")))
                    entry["own_mcc"][name] = c["results"]["gene_level_conditions_where_wt_grows"].get("mcc", {}).get("point")
                runs = [arms[n] for n, _ in ARMS] + [ref]
                names = [n for n, _ in ARMS] + [label]
                sims, fit, genes, grows = union_n(runs)
                st, ft = ref["params"]["growth_threshold"], ref["params"]["fitness_threshold"]
                m = {n: TC.mcc_from_counts(TC.gene_counts(s, fit, st, ft).sum(0)) for n, s in zip(names, sims)}
                gap = m[label] - m["B0"]
                entry["four_way"] = {"union_genes": len(genes), "conditions_all_grow": len(grows),
                                     "n_observations": int(np.isfinite(fit).sum()), "mcc": m,
                                     "share_of_gap_closed": ({"U'": (m["U'"] - m["B0"]) / gap, "M": (m["M"] - m["B0"]) / gap}
                                                             if gap > 0 else None)}
                # Error overlap on common genes (U' vs curated) and the curated model's calls on genes only it contains.
                from gembench.comparison import aligned_pairs
                su, sc, fitc, cgenes, cgrows, _ = aligned_pairs(arms["U'"], ref)
                ok = np.isfinite(su) & np.isfinite(sc) & np.isfinite(fitc)
                obs = np.where(ok, fitc <= ft, False)
                wrong_u = ok & ((su < st) != obs)
                wrong_c = ok & ((sc < st) != obs)
                n_ok = int(ok.sum())
                entry["error_overlap_common_genes_Uprime_vs_curated"] = {
                    "n_genes": len(cgenes), "n_conditions": len(cgrows), "n_observations": n_ok,
                    "wrong_Uprime": int(wrong_u.sum()), "wrong_curated": int(wrong_c.sum()), "wrong_both": int((wrong_u & wrong_c).sum()),
                    "expected_both_if_independent": float(wrong_u.sum() * wrong_c.sum() / n_ok) if n_ok else None}
                draft_genes = set().union(*[set(arms[n]["browser_genes"]) for n, _ in ARMS])
                only = [g for g in ref["browser_genes"] if g not in draft_genes]
                gi = {g: i for i, g in enumerate(ref["browser_genes"])}
                cols = [j for j, w in enumerate(ref["wt_growth"]) if np.isfinite(w) and w >= st
                        and ref["conditions"][j] in set(cgrows)]
                sim_o = ref["sim_growth"][np.ix_([gi[g] for g in only], cols)]
                fit_o = ref["fitness"][np.ix_([gi[g] for g in only], cols)]
                ok_o = np.isfinite(sim_o) & np.isfinite(fit_o)
                imp = ok_o & (sim_o < st)
                entry["curated_only_genes"] = {"n_genes": len(only), "n_conditions": len(cols),
                                               "conditions": [ref["conditions"][j] for j in cols],
                                               "important_calls": int(imp.sum()),
                                               "important_calls_confirmed": int((imp & (fit_o <= ft)).sum())}
                # Why models do not grow: conditions without growth that lack an exchange reaction for the carbon source.
                import pandas as pd
                gaps = {}
                for name, d in list(arm_dirs.items()) + [(label, ref_dir)]:
                    ct = pd.read_table(os.path.join(d, "conditions.tsv"), keep_default_na=False)
                    no = ct[~ct["wt_grows"].astype(str).str.lower().eq("true")]
                    gaps[name] = {"no_growth": int(len(no)), "no_growth_missing_exchange": int((no["missing_exchanges"].astype(str) != "").sum())}
                entry["no_growth_conditions"] = gaps
                entry["pairwise_curated_minus_arm"] = {}
                for n in names[:-1]:
                    un = TC.union_delta(arms[n], ref, n_boot=args.n_boot)
                    co = TC.common_delta(arms[n], ref, n_boot=args.n_boot)
                    entry["pairwise_curated_minus_arm"][n] = {
                        "union": un,
                        "common": {"n_common_genes": co["n_common_genes"], "n_conditions_both_grow": co["n_conditions_both_grow"],
                                   "mcc_arm": co["mcc"]["A"], "mcc_curated": co["mcc"]["B"],
                                   "difference": co["paired_mcc_difference_B_minus_A"]}}
            report["organisms"][f"{org}/{label}"] = entry
            print(json.dumps({k: entry[k] for k in entry if k in ("coverage", "own_mcc", "four_way")}, indent=1, default=float))
    with open(args.out, "w") as fh:
        json.dump(report, fh, indent=1, default=float)
        fh.write("\n")


if __name__ == "__main__":
    main()
