"""Check 3e: iAH991's gene map (BT_0554 -> BT0554) misses loci that the Fitness Browser lists with the underscore
(BT_0823, BT_2070 with fitness data), the mirror image of the iSO783 problem. This script scores the two missing genes
in the exploratory medium with this check's protocol implementation, appends them to the saved exploratory matrix
(matrices_iAH991_exploratory_gene_fallback.npz) and recomputes the B. theta comparison.

Usage: python -I rescore_iah991_gene_fallback.py
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
from common import ROOT, p, read_model  # noqa: E402
from protocol import ST, complete_medium, fitness_matrix, fitness_table, mapped_conditions, simulate  # noqa: E402

SAVED = os.path.join(ROOT, "results/transfer_v1/reference_models/Btheta/iAH991_exploratory_B12_H2S/matrices.npz")


def main():
    conds = mapped_conditions("Btheta")
    sysn = set(pd.read_table(p("data/fitness_browser/Btheta/genes.tsv"), dtype=str, keep_default_na=False)["sysName"])
    fit_index = set(fitness_table("Btheta").index)
    z = np.load(SAVED)
    have = set(str(x) for x in z["browser_genes"])
    m = read_model("models/curated/iAH991/iAH991_bigg_view.xml.gz")
    extra = []
    for g in m.genes:
        s = g.id.replace("_", "")
        if s in sysn:
            continue
        if g.id in sysn and g.id in fit_index and g.id not in have:
            extra.append((g.id, g.id))
    sup = {"cbl1": -0.001, "h2s": -1000.0}
    complete_medium(m, ["Varel_Bryant_medium"], sup)
    sim, wt = simulate(m, [a for a, _ in extra], conds, sup)
    fit = fitness_matrix("Btheta", conds, [b for _, b in extra])
    grows = wt >= ST
    rep = {"extra_genes": [{"model_gene": a, "sysName": b, "ko_no_growth_conditions": int((sim[i, grows] < ST).sum()),
                            "fitness_le_-2_conditions": int((fit[i, grows] <= -2).sum()), "n_growing": int(grows.sum())}
                           for i, (a, b) in enumerate(extra)],
           "wt_matches_saved": bool(np.allclose(wt, z["wt_growth"], atol=1e-9))}
    out = os.path.join(HERE, "matrices_iAH991_exploratory_gene_fallback.npz")
    np.savez_compressed(out, sim_growth=np.vstack([z["sim_growth"], sim]), wt_growth=z["wt_growth"],
                        fitness=np.vstack([z["fitness"], fit]), model_genes=np.concatenate([z["model_genes"], [a for a, _ in extra]]),
                        browser_genes=np.concatenate([z["browser_genes"], [b for _, b in extra]]), conditions=z["conditions"])
    json.dump(rep, open(os.path.join(HERE, "rescore_iah991_gene_fallback.json"), "w"), indent=1)
    print(rep, flush=True)


if __name__ == "__main__":
    main()
