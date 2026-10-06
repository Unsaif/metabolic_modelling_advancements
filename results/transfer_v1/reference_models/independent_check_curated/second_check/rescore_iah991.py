"""Check 3d: re-run the iAH991 scoring with this check's own protocol implementation.

  * strict protocol: wild-type growth in the 25 mapped conditions (expected: none);
  * exploratory (cob(I)alamin -0.001, sulfide -1000): full knockout matrix, compared with the saved
    results/transfer_v1/reference_models/Btheta/iAH991_exploratory_B12_H2S/matrices.npz;
  * sensitivity, B12 only (cob(I)alamin -0.001, no sulfide): full knockout matrix, saved as
    iah991_B12_only.npz, to see how much the sulfide choice moves iAH991's calls.

Gene map: BT_0554 -> BT0554, genes with fitness data, model gene order.

Usage: python -I rescore_iah991.py   (writes rescore_iah991.json, iah991_B12_only.npz)
"""
from __future__ import annotations

import json
import os
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
    keys = [f"{c['name']} | {c['media']}" for c in conds]
    sysn = set(pd.read_table(p("data/fitness_browser/Btheta/genes.tsv"), dtype=str, keep_default_na=False)["sysName"])
    fit_index = set(fitness_table("Btheta").index)
    rep = {}
    m = read_model("models/curated/iAH991/iAH991_bigg_view.xml.gz")
    pairs, seen = [], set()
    for g in m.genes:
        s = g.id.replace("_", "")
        if s in sysn and s in fit_index and s not in seen:
            pairs.append((g.id, s)); seen.add(s)
    genes, bgenes = [a for a, _ in pairs], [b for _, b in pairs]
    rep["n_genes"] = len(genes)
    saved = np.load(SAVED)
    assert [str(x) for x in saved["conditions"]] == keys
    order = {g: i for i, g in enumerate(str(x) for x in saved["browser_genes"])}
    rep["same_gene_set_as_saved"] = sorted(order) == sorted(bgenes)
    fit = fitness_matrix("Btheta", conds, bgenes)
    idx = [order[g] for g in bgenes]
    rep["max_abs_diff_fitness_vs_saved"] = float(np.nanmax(np.abs(fit - saved["fitness"][idx])))
    # strict
    ms = read_model("models/curated/iAH991/iAH991_bigg_view.xml.gz")
    rep["strict_completion_added"] = complete_medium(ms, ["Varel_Bryant_medium"])
    _, wt0 = simulate(ms, genes, conds, None, knockouts=False)
    rep["strict_wt_grows"] = int((wt0 >= ST).sum())
    rep["strict_max_wt"] = float(wt0.max())
    print("strict:", rep["strict_wt_grows"], rep["strict_max_wt"], flush=True)
    # exploratory
    sup = {"cbl1": -0.001, "h2s": -1000.0}
    me = read_model("models/curated/iAH991/iAH991_bigg_view.xml.gz")
    rep["exploratory_completion_added"] = complete_medium(me, ["Varel_Bryant_medium"], sup)
    sim, wt = simulate(me, genes, conds, sup)
    rep["exploratory_reproduction"] = {"max_abs_diff_wt": float(np.max(np.abs(wt - saved["wt_growth"]))),
                                       "max_abs_diff_sim": float(np.max(np.abs(sim - saved["sim_growth"][idx]))),
                                       "binary_calls_differ": int(((sim < ST) != (saved["sim_growth"][idx] < ST)).sum())}
    print("exploratory reproduction:", rep["exploratory_reproduction"], flush=True)
    # B12 only
    sup2 = {"cbl1": -0.001}
    mb = read_model("models/curated/iAH991/iAH991_bigg_view.xml.gz")
    complete_medium(mb, ["Varel_Bryant_medium"], sup2)
    sim2, wt2 = simulate(mb, genes, conds, sup2)
    both = (wt >= ST) & (wt2 >= ST)
    rep["B12_only"] = {"wt_grows": int((wt2 >= ST).sum()), "wt_range": [float(wt2[wt2 >= ST].min()), float(wt2.max())],
                       "important_calls_exploratory": int((sim[:, both] < ST).sum()),
                       "important_calls_B12_only": int((sim2[:, both] < ST).sum()),
                       "calls_changed": int(((sim < ST) != (sim2 < ST))[:, both].sum())}
    print("B12 only:", rep["B12_only"], flush=True)
    np.savez_compressed(os.path.join(HERE, "iah991_B12_only.npz"), sim_growth=sim2, wt_growth=wt2, fitness=fit,
                        model_genes=np.array(genes), browser_genes=np.array(bgenes), conditions=np.array(keys))
    json.dump(rep, open(os.path.join(HERE, "rescore_iah991.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
