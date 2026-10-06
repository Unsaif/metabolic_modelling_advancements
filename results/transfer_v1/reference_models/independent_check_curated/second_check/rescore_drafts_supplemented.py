"""Check 3c: would the B. theta drafts' calls change if they got the same exploratory supplement as iAH991
(cob(I)alamin -0.001, sulfide -1000)?

For each arm (B0, UNQ = U', M): the saved arm model (results/transfer_v1/development/Btheta/<arm>/model.xml.gz) is
scored with this check's protocol implementation, first as the study did (to show the implementation reproduces the
saved matrix), then with the supplement. The gene list and its Fitness Browser mapping are taken from the saved
matrices. Writes drafts_supplemented/<arm>.npz (same layout as matrices.npz) and rescore_drafts_supplemented_<arm>.json.

Usage: python -I rescore_drafts_supplemented.py <arm>
"""
from __future__ import annotations

import gzip
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
from common import ROOT  # noqa: E402
from protocol import ST, complete_medium, fitness_matrix, mapped_conditions, simulate  # noqa: E402

import cobra  # noqa: E402

SUP = {"cbl1": -0.001, "h2s": -1000.0}


def main():
    arm = sys.argv[1]
    d = os.path.join(ROOT, "results", "transfer_v1", "development", "Btheta", arm)
    z = np.load(os.path.join(d, "matrices.npz"))
    genes = [str(x) for x in z["model_genes"]]
    bgenes = [str(x) for x in z["browser_genes"]]
    conds = mapped_conditions("Btheta")
    keys = [f"{c['name']} | {c['media']}" for c in conds]
    assert keys == [str(x) for x in z["conditions"]], "condition order differs"
    fit = fitness_matrix("Btheta", conds, bgenes)
    rep = {"arm": arm, "n_genes": len(genes), "max_abs_diff_fitness_vs_saved": float(np.nanmax(np.abs(fit - z["fitness"])))}
    with gzip.open(os.path.join(d, "model.xml.gz"), "rt") as fh:
        model = cobra.io.read_sbml_model(fh)
    model.solver = "glpk"
    for ex in model.exchanges:
        ex.bounds = (0.0, 1000.0)
    rep["completion_added_strict"] = complete_medium(model, ["Varel_Bryant_medium"])
    sim, wt = simulate(model, genes, conds)
    rep["reproduction"] = {"max_abs_diff_wt": float(np.max(np.abs(wt - z["wt_growth"]))),
                           "max_abs_diff_sim": float(np.max(np.abs(sim - z["sim_growth"]))),
                           "binary_calls_differ": int(((sim < ST) != (z["sim_growth"] < ST)).sum()),
                           "wt_grows_differ": int(((wt >= ST) != (z["wt_growth"] >= ST)).sum())}
    print(arm, "reproduction", rep["reproduction"], flush=True)
    with gzip.open(os.path.join(d, "model.xml.gz"), "rt") as fh:
        model = cobra.io.read_sbml_model(fh)
    model.solver = "glpk"
    for ex in model.exchanges:
        ex.bounds = (0.0, 1000.0)
    rep["completion_added_supplemented"] = complete_medium(model, ["Varel_Bryant_medium"], SUP)
    sim2, wt2 = simulate(model, genes, conds, SUP)
    grow_both = (wt >= ST) & (wt2 >= ST)
    changed = (sim < ST) != (sim2 < ST)
    rep["supplemented"] = {"wt_grows": int((wt2 >= ST).sum()), "wt_grows_strict": int((wt >= ST).sum()),
                           "conditions_gained": [k for k, a, b in zip(keys, wt, wt2) if a < ST <= b],
                           "conditions_lost": [k for k, a, b in zip(keys, wt, wt2) if b < ST <= a],
                           "calls_changed_in_conditions_both_grow": int((changed[:, grow_both]).sum()),
                           "genes_with_changed_calls": sorted({bgenes[i] for i in np.where(changed[:, grow_both].any(axis=1))[0]}),
                           "to_not_important": int(((sim < ST) & (sim2 >= ST))[:, grow_both].sum()),
                           "to_important": int(((sim >= ST) & (sim2 < ST))[:, grow_both].sum())}
    print(arm, "supplemented", {k: v for k, v in rep["supplemented"].items() if k != "genes_with_changed_calls"},
          "genes:", rep["supplemented"]["genes_with_changed_calls"][:30], flush=True)
    os.makedirs(os.path.join(HERE, "drafts_supplemented"), exist_ok=True)
    np.savez_compressed(os.path.join(HERE, "drafts_supplemented", f"{arm}.npz"), sim_growth=sim2, wt_growth=wt2,
                        fitness=z["fitness"], model_genes=z["model_genes"], browser_genes=z["browser_genes"],
                        conditions=z["conditions"])
    json.dump(rep, open(os.path.join(HERE, f"rescore_drafts_supplemented_{arm}.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
