#!/usr/bin/env python3
"""Independent checks of (a) coverage across the transfer study, (b) the reference runs' inputs and protocol, and
(c) the 94% figure. Reads cards, matrices and raw Fitness Browser tables with json/numpy/pandas only; no gembench or
scripts/ imports and no LP solve.

(a) For the twelve organisms with a working draft (4 development + 8 new: all evaluable organisms except
    M. tuberculosis and R. palustris), coverage of B0 and U' (=UNQ) from each arm's card (results.condition_level,
    results.counts.conditions_total) and, independently, from the arm's matrices (wt_growth >= threshold).
(b) Protocol equivalence of the reference runs with the P. putida development arms; model-file hashes; identity gene
    mapping; and the fitness matrices of both reference runs rebuilt from the raw Fitness Browser tables
    (carbon-source experiments with condition_2 empty or DMSO, media present in the media table, grouped by
    compound x medium, NaN-aware mean over replicates, carbon source mapped to BiGG ids).
(c) The 94% figure from paper1_derived_numbers.json, and the underlying counts recomputed from the B0 and UNQ
    matrices of the ten evaluable new organisms (common genes, conditions where both grow).
Writes check_coverage_and_inputs.json next to this file.
"""
from __future__ import annotations

import hashlib
import json
import os
import statistics

import numpy as np
import pandas as pd

ROOT = "/home/claude/mma"
RES = os.path.join(ROOT, "results", "transfer_v1")
OUT = os.path.join(RES, "reference_models", "independent_check")

TWELVE = [("development", "Btheta"), ("development", "Putida"), ("development", "MR1"), ("development", "Smeli"),
          ("evaluation_panel_A", "Cola"), ("evaluation_panel_A", "Dino"), ("evaluation_panel_A", "Dyella79"),
          ("evaluation_panel_B", "Caulo"), ("evaluation_panel_B", "Cup4G11"), ("evaluation_panel_B", "Marino"),
          ("evaluation_panel_B", "PV4"), ("evaluation_panel_B", "SB2B")]
EXCLUDED = [("evaluation_panel_A", "MycoTube"), ("evaluation_panel_A", "RPal_CGA009")]
TEN_NEW = [x for x in TWELVE if x[0] != "development"] + EXCLUDED


def card(phase, org, arm):
    return json.load(open(os.path.join(RES, phase, org, arm, "card.json")))


def npz(d):
    with np.load(os.path.join(d, "matrices.npz"), allow_pickle=False) as z:
        return {k: z[k].copy() for k in z.files}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------- (a) coverage
def coverage():
    rows = []
    for phase, org in TWELVE + EXCLUDED:
        rec = {"phase": phase, "org": org}
        for arm, lab in (("B0", "B0"), ("UNQ", "U")):
            c = card(phase, org, arm)
            m = npz(os.path.join(RES, phase, org, arm))
            gt = c["protocol"]["params"]["growth_threshold"]
            w = m["wt_growth"]
            rec[f"mapped_{lab}"] = c["results"]["condition_level"]["n_conditions_mapped"]
            rec[f"grows_{lab}_card"] = c["results"]["condition_level"]["n_conditions_wt_grows"]
            rec[f"grows_{lab}_matrix"] = int((np.isfinite(w) & (w >= gt)).sum())
            rec[f"mapped_{lab}_matrix"] = int(len(w))
            rec[f"conditions_total_{lab}"] = c["results"]["counts"]["conditions_total"]
        rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "coverage_by_organism.tsv"), sep="\t", index=False)
    t = df[df.apply(lambda r: (r["phase"], r["org"]) in TWELVE, axis=1)]
    assert (t["grows_B0_card"] == t["grows_B0_matrix"]).all() and (t["grows_U_card"] == t["grows_U_matrix"]).all()
    assert (t["mapped_B0"] == t["mapped_U"]).all() and (t["mapped_B0"] == t["mapped_B0_matrix"]).all()
    frac = (t["grows_B0_card"] / t["mapped_B0"]).tolist()
    frac_u = (t["grows_U_card"] / t["mapped_U"]).tolist()
    total_screened = int(t["conditions_total_B0"].sum())
    q = np.percentile(frac, [25, 75])
    # table1 cross-check
    t1 = pd.read_csv(os.path.join(RES, "table1_organisms.tsv"), sep="\t", dtype=str)
    t1m = {r["org"]: r for _, r in t1.iterrows()}
    t1_agree = all(str(int(r["mapped_B0"])) == t1m[r["org"]]["mapped_conditions"] and
                   str(int(r["grows_B0_card"])) == t1m[r["org"]]["grows_B0"] and
                   str(int(r["grows_U_card"])) == t1m[r["org"]]["grows_U"] for _, r in t.iterrows())
    all14 = df
    return {
        "n_organisms": int(len(t)),
        "B0_grows_total": int(t["grows_B0_card"].sum()), "U_grows_total": int(t["grows_U_card"].sum()),
        "mapped_total": int(t["mapped_B0"].sum()),
        "conditions_total_screened_in_mapped_media": total_screened,
        "B0_share_of_mapped": float(t["grows_B0_card"].sum() / t["mapped_B0"].sum()),
        "B0_share_of_all_screened_in_mapped_media": float(t["grows_B0_card"].sum() / total_screened),
        "B0_fraction_min": float(min(frac)), "B0_fraction_max": float(max(frac)), "B0_fraction_median": float(statistics.median(frac)),
        "B0_fraction_iqr": [float(q[0]), float(q[1])],
        "U_fraction_median": float(statistics.median(frac_u)),
        "per_organism_B0_fraction": {r["org"]: round(r["grows_B0_card"] / r["mapped_B0"], 4) for _, r in t.iterrows()},
        "per_organism_mapped_vs_total": {r["org"]: [int(r["mapped_B0"]), int(r["conditions_total_B0"])] for _, r in t.iterrows()},
        "organisms_where_U_grows_more": {r["org"]: int(r["grows_U_card"] - r["grows_B0_card"]) for _, r in t.iterrows()
                                         if r["grows_U_card"] != r["grows_B0_card"]},
        "table1_matches_cards": bool(t1_agree),
        "with_MycoTube_and_RPal": {"B0_grows": int(all14["grows_B0_card"].sum()), "mapped": int(all14["mapped_B0"].sum())},
    }


# --------------------------------------------------------------------------------------- (b) inputs and protocol
def fb_conditions(org_dir, media_table, carbon_table):
    exp = pd.read_csv(os.path.join(org_dir, "experiments.tsv"), sep="\t", dtype=str, keep_default_na=False)
    fit = pd.read_csv(os.path.join(org_dir, "fit_logratios.tsv"), sep="\t", dtype=str, keep_default_na=False)
    meta_cols = ["orgId", "locusId", "sysName", "geneName", "desc"]
    exp_cols = [c for c in fit.columns if c not in meta_cols]
    # fitness columns are "<expName> <description>"; map expName -> column
    col_of = {c.split(" ", 1)[0]: c for c in exp_cols}
    media = set(pd.read_csv(media_table, sep="\t", dtype=str, keep_default_na=False)["media"])
    cmap = pd.read_csv(carbon_table, sep="\t", dtype=str, keep_default_na=False).set_index("condition")
    sel = exp[(exp["expGroup"] == "carbon source") & (exp["condition_2"].isin(["", "Dimethyl Sulfoxide"]))]
    groups, all_cs = {}, set()
    for _, r in sel.iterrows():
        all_cs.add((r["condition_1"], r["media"]))
        if r["media"] not in media or r["expName"] not in col_of:
            continue
        groups.setdefault(f"{r['condition_1']} | {r['media']}", []).append(r["expName"])
    mapped = {}
    for key, exps in groups.items():
        name = key.split(" | ")[0]
        ids = [x for x in cmap.loc[name, "bigg_ids"].split(";") if x] if name in cmap.index else []
        if ids:
            mapped[key] = exps
    vals = fit[[col_of[e] for exps in groups.values() for e in exps]].replace("", np.nan).astype(float)
    vals.index = fit["sysName"].where(fit["sysName"] != "", fit["locusId"])
    mat = {key: vals[[col_of[e] for e in exps]].mean(axis=1, skipna=True) for key, exps in mapped.items()}
    return groups, mapped, pd.DataFrame(mat), all_cs, sel


def check_reference_inputs():
    out = {}
    ref_dirs = {"iJN1463": ("Putida", os.path.join(RES, "reference_models", "Putida", "iJN1463")),
                "iML1515": ("Keio", os.path.join(RES, "reference_models", "Keio", "iML1515"))}
    arm_params = {a: card("development", "Putida", a)["protocol"] for a in ("B0", "UNQ", "M")}
    for label, (org, d) in ref_dirs.items():
        c = json.load(open(os.path.join(d, "card.json")))
        p = c["protocol"]
        m = npz(d)
        rec = {"params_equal_to_putida_arms": all(p["params"] == ap["params"] for ap in arm_params.values()),
               "media_table_equal": all(p["media_mapping"] == ap["media_mapping"] for ap in arm_params.values()),
               "carbon_table_equal": all(p["carbon_source_mapping"] == ap["carbon_source_mapping"] for ap in arm_params.values()),
               "applied_transforms": p["applied"], "gene_mapping": p["gene_mapping"],
               "identity_mapping_in_matrices": bool(np.array_equal(m["model_genes"], m["browser_genes"])),
               "model_file": c["model"]["file"], "card_sha256": c["model"]["sha256"],
               "file_sha256": sha256(os.path.join(ROOT, c["model"]["file"])),
               "software": c.get("software"), "model_source": c["model"]["source"]}
        rec["sha256_matches"] = rec["card_sha256"] == rec["file_sha256"]
        # fitness rebuilt from raw Fitness Browser tables
        groups, mapped, mat, all_cs, sel = fb_conditions(os.path.join(ROOT, "data", "fitness_browser", org),
                                                        os.path.join(ROOT, p["media_mapping"]),
                                                        os.path.join(ROOT, p["carbon_source_mapping"]))
        conds = [str(x) for x in m["conditions"]]
        genes = [str(x) for x in m["browser_genes"]]
        rec["raw_condition_groups_in_mapped_media"] = len(groups)
        rec["raw_mapped_conditions"] = len(mapped)
        rec["raw_carbon_source_conditions_any_medium"] = len(all_cs)
        rec["run_conditions_equal_raw_mapped"] = sorted(conds) == sorted(mapped)
        rec["unmapped_conditions"] = sorted(set(groups) - set(mapped))
        sub = mat.loc[genes, conds].to_numpy(dtype=float)
        f = m["fitness"]
        same_nan = np.array_equal(np.isfinite(sub), np.isfinite(f))
        diff = np.nanmax(np.abs(sub - f)) if np.isfinite(f).any() else float("nan")
        rec["fitness_rebuilt_from_raw"] = {"cells": int(np.isfinite(f).sum()), "nan_pattern_equal": bool(same_nan),
                                           "max_abs_diff": float(diff)}
        out[label] = rec
    # Putida draft arms: same check, so the four runs' fitness is known to come from the same table
    groups, mapped, mat, _, _ = fb_conditions(os.path.join(ROOT, "data", "fitness_browser", "Putida"),
                                              os.path.join(ROOT, arm_params["B0"]["media_mapping"]),
                                              os.path.join(ROOT, arm_params["B0"]["carbon_source_mapping"]))
    arms = {}
    for a in ("B0", "UNQ", "M"):
        m = npz(os.path.join(RES, "development", "Putida", a))
        sub = mat.loc[[str(g) for g in m["browser_genes"]], [str(c) for c in m["conditions"]]].to_numpy(dtype=float)
        arms[a] = {"nan_pattern_equal": bool(np.array_equal(np.isfinite(sub), np.isfinite(m["fitness"]))),
                   "max_abs_diff": float(np.nanmax(np.abs(sub - m["fitness"])))}
    out["putida_arms_fitness_rebuilt_from_raw"] = arms
    # software of the transfer arms vs the reference runs
    out["software_putida_arms"] = card("development", "Putida", "B0").get("software")
    return out


# ------------------------------------------------------------------------------------------------ (c) 94 percent
def changed_predictions():
    d = json.load(open(os.path.join(RES, "paper1_derived_numbers.json")))["U_prime_vs_B0"]
    tot = {"to_important": 0, "to_important_agree": 0, "to_unimportant": 0, "to_unimportant_agree": 0, "genes": 0}
    per = {}
    for phase, org in TEN_NEW:
        a, b = (npz(os.path.join(RES, phase, org, x)) for x in ("B0", "UNQ"))
        ca = card(phase, org, "B0")["protocol"]["params"]
        gt, ft = ca["growth_threshold"], ca["fitness_threshold"]
        ga = {str(g): i for i, g in enumerate(a["browser_genes"])}
        gb = {str(g): i for i, g in enumerate(b["browser_genes"])}
        cb = {str(c): j for j, c in enumerate(b["conditions"])}
        genes = sorted(set(ga) & set(gb))
        conds = [str(c) for j, c in enumerate(a["conditions"]) if str(c) in cb and a["wt_growth"][j] >= gt
                 and b["wt_growth"][cb[str(c)]] >= gt]
        ia, ib = [ga[g] for g in genes], [gb[g] for g in genes]
        ja = [{str(c): j for j, c in enumerate(a["conditions"])}[c] for c in conds]
        jb = [cb[c] for c in conds]
        sa, sb = a["sim_growth"][np.ix_(ia, ja)], b["sim_growth"][np.ix_(ib, jb)]
        fa, fb = a["fitness"][np.ix_(ia, ja)], b["fitness"][np.ix_(ib, jb)]
        assert np.array_equal(fa, fb, equal_nan=True)
        ok = np.isfinite(fa) & np.isfinite(sa) & np.isfinite(sb)
        imp = fa <= ft
        to_imp = ok & (sa >= gt) & (sb < gt)
        to_un = ok & (sa < gt) & (sb >= gt)
        r = {"to_important": int(to_imp.sum()), "to_important_agree": int((to_imp & imp).sum()),
             "to_unimportant": int(to_un.sum()), "to_unimportant_agree": int((to_un & ~imp).sum()),
             "genes": int((to_imp | to_un).any(axis=1).sum())}
        per[org] = r
        for k in tot:
            tot[k] += r[k]
    share = (tot["to_unimportant_agree"] + tot["to_important_agree"]) / (tot["to_unimportant"] + tot["to_important"])
    return {"derived_file": {k: d[k] for k in ("genes_changed_common", "to_important", "to_important_agree",
                                                "to_unimportant", "to_unimportant_agree")},
            "derived_file_share": (d["to_unimportant_agree"] + d["to_important_agree"]) / (d["to_unimportant"] + d["to_important"]),
            "recomputed": tot, "recomputed_share": share, "per_organism": per}


def main():
    out = {"coverage": coverage(), "reference_inputs": check_reference_inputs(), "changed_predictions": changed_predictions()}
    with open(os.path.join(OUT, "check_coverage_and_inputs.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
        fh.write("\n")
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
