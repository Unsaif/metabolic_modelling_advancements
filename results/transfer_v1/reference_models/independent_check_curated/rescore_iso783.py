"""Re-run the iSO783 scoring with independent code, then again with the gene map corrected.

Part A reproduces results/transfer_v1/reference_models/MR1/iSO783/matrices.npz from the BiGG view with this
script's own implementation of the protocol (no gembench import): exchanges closed, medium completion (exchange +
gene-less uptake for medium components present in the cytosol, pantothenate/folate/bicarbonate excluded), the
Fitness Browser medium at the study's bounds, each carbon source at -10, wild-type growth, every mapped gene knocked
out where the wild type grows; fitness averaged over the condition's experiments.

Part B repeats it with the identity rule extended by one fallback: a model locus tag SOnnnn that is not a Fitness
Browser sysName is mapped to SO_nnnn when that is one (the Fitness Browser lists 124 MR-1 loci in the RefSeq
form). The corrected matrices are saved as matrices_iSO783_corrected_gene_map.npz for recompute_comparison.py.

Usage: python -I rescore_iso783.py
"""
from __future__ import annotations

import gzip
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import MODELS, OUT, ROOT, dump, fb_conditions, medium_bounds, p, read_model  # noqa: E402
from check_energy_gate import complete_medium, egc_values, CURRENCIES_BIGG  # noqa: E402

ST = 1e-3


def fitness_table(org):
    fit = pd.read_table(p(f"data/fitness_browser/{org}/fit_logratios.tsv"), dtype=str, keep_default_na=False)
    meta = [c for c in ["orgId", "locusId", "sysName", "geneName", "desc"] if c in fit.columns]
    cols = [c for c in fit.columns if c not in meta]
    f = fit[cols].apply(pd.to_numeric, errors="coerce")
    f.columns = [c.split(" ")[0] for c in cols]
    f.index = fit["sysName"].where(fit["sysName"] != "", fit["locusId"]).values
    return f


def growth(model):
    v = model.slim_optimize(error_value=float("nan"))
    st = model.solver.status
    if st == "infeasible":
        return 0.0
    if st != "optimal" or not np.isfinite(v):
        raise RuntimeError(st)
    return float(v)


def simulate(model, genes, conds, medium):
    sim = np.zeros((len(genes), len(conds)))
    wt = np.zeros(len(conds))
    mb = medium_bounds(medium)
    for j, c in enumerate(conds):
        with model:
            for ex in model.exchanges:
                ex.bounds = (0.0, 1000.0)
            for rid, lb in mb.items():
                if rid in model.reactions:
                    model.reactions.get_by_id(rid).lower_bound = lb
            for b in c["bigg_ids"]:
                if f"EX_{b}_e" in model.reactions:
                    model.reactions.get_by_id(f"EX_{b}_e").lower_bound = -10.0
            wt[j] = growth(model)
            if wt[j] >= ST:
                for i, g in enumerate(genes):
                    with model:
                        model.genes.get_by_id(g).knock_out()
                        sim[i, j] = growth(model)
    return sim, wt


def main():
    cfg = MODELS["iSO783"]
    org, medium = cfg["org"], cfg["medium"]
    fit = fitness_table(org)
    sysnames = set(pd.read_table(p(f"data/fitness_browser/{org}/genes.tsv"), dtype=str, keep_default_na=False)["sysName"])
    conds = [c for c in fb_conditions(org) if c["bigg_ids"]]
    fmat = pd.DataFrame({f"{c['name']} | {c['media']}": fit[c["experiments"]].mean(axis=1, skipna=True) for c in conds})
    model = read_model(cfg["view"])
    for ex in model.exchanges:
        ex.bounds = (0.0, 1000.0)
    added = complete_medium(model, sorted({c["media"] for c in conds}))
    for ex in model.exchanges:
        ex.bounds = (0.0, 1000.0)
    gate = egc_values(model, CURRENCIES_BIGG)

    def gene_pairs(rule):
        pairs, seen = [], set()
        for g in model.genes:
            s = rule(g.id)
            if s is not None and s in fit.index and s not in seen:
                pairs.append((g.id, s)); seen.add(s)
        return pairs

    identity = lambda g: g if g in sysnames else None                                     # noqa: E731

    def corrected(g):
        if g in sysnames:
            return g
        m = re.match(r"^SO(\d{4}[A-Za-z]?)$", g)
        return f"SO_{m.group(1)}" if m and f"SO_{m.group(1)}" in sysnames else None

    rep = {"medium_completion_added": added, "energy_check": gate, "conditions": [f"{c['name']} | {c['media']}" for c in conds]}
    saved = np.load(p("results/transfer_v1/reference_models/MR1/iSO783/matrices.npz"))
    # Part A
    pa = gene_pairs(identity)
    sim, wt = simulate(model, [g for g, _ in pa], conds, medium)
    bg = [s for _, s in pa]
    order = {g: i for i, g in enumerate(str(x) for x in saved["browser_genes"])}
    same_genes = sorted(bg) == sorted(order)
    idx = [order[g] for g in bg]
    s_saved, f_saved = saved["sim_growth"][idx], saved["fitness"][idx]
    fm = fmat.loc[bg].to_numpy()
    rep["part_A_reproduction"] = {
        "n_genes": len(bg), "same_gene_set_as_saved": same_genes,
        "conditions_same_order": [str(x) for x in saved["conditions"]] == rep["conditions"],
        "max_abs_diff_wt": float(np.max(np.abs(wt - saved["wt_growth"]))),
        "max_abs_diff_sim": float(np.max(np.abs(sim - s_saved))),
        "binary_calls_differ": int(((sim < ST) != (s_saved < ST)).sum()),
        "max_abs_diff_fitness": float(np.nanmax(np.abs(fm - f_saved))),
    }
    print("A:", rep["part_A_reproduction"], flush=True)
    # Part B
    pb = gene_pairs(corrected)
    extra = [(g, s) for g, s in pb if s not in set(bg)]
    sim_x, _ = simulate(model, [g for g, _ in extra], conds, medium)
    genes_b = [g for g, _ in pa] + [g for g, _ in extra]
    bg_b = bg + [s for _, s in extra]
    sim_b = np.vstack([sim, sim_x])
    fit_b = fmat.loc[bg_b].to_numpy()
    out = os.path.join(OUT, "matrices_iSO783_corrected_gene_map.npz")
    np.savez_compressed(out, sim_growth=sim_b, wt_growth=wt, fitness=fit_b, model_genes=np.array(genes_b),
                        browser_genes=np.array(bg_b), conditions=np.array(rep["conditions"]))
    grows = wt >= ST
    rep["part_B_extra_genes"] = [{"model_gene": g, "sysName": s,
                                  "ko_no_growth_in_growing_conditions": int((sim_x[k, grows] < ST).sum()),
                                  "fitness_le_-2_in_growing_conditions": int((fmat.loc[s].to_numpy()[grows] <= -2).sum()),
                                  "n_growing_conditions": int(grows.sum())}
                                 for k, (g, s) in enumerate(extra)]
    rep["part_B_matrices"] = os.path.relpath(out, ROOT)
    print("B:", rep["part_B_extra_genes"], flush=True)
    dump(rep, "rescore_iso783.json")


if __name__ == "__main__":
    main()
