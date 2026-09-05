"""Evaluate a universe patch set: for every organism/model with an unpatched and a patched result directory,
compare predictions gene by gene (same genes, same conditions) and report what changed and whether the
changes agree with the fitness data. This is the verifier step of the fix-at-source loop.

Usage: python scripts/evaluate_patch.py v0.1
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from gembench import metrics as M  # noqa: E402
from gembench.fitness_browser import load_organism  # noqa: E402

RES = os.path.join(ROOT, "results", "carbon_fitness_multi")


def load(d):
    z = np.load(os.path.join(d, "matrices.npz"), allow_pickle=True)
    return z["sim_growth"], z["fitness"], z["wt_growth"], list(z["browser_genes"]), list(z["conditions"])


def main() -> None:
    ver = sys.argv[1] if len(sys.argv) > 1 else "v0.1"
    rows, changes = [], []
    for pd_dir in sorted(glob.glob(os.path.join(RES, "*", f"*__patched-{ver}"))):
        base = pd_dir.replace(f"__patched-{ver}", "")
        if not os.path.exists(os.path.join(base, "matrices.npz")):
            continue
        org = os.path.basename(os.path.dirname(pd_dir)); label = os.path.basename(base)
        s0, f0, w0, g0, c0 = load(base); s1, f1, w1, g1, c1 = load(pd_dir)
        assert c0 == c1, "patched and unpatched runs must have the same conditions"
        common = [g for g in g0 if g in set(g1)]
        dropped = [g for g in g0 if g not in set(g1)] + [g for g in g1 if g not in set(g0)]
        i0 = [g0.index(g) for g in common]; i1 = [g1.index(g) for g in common]
        s0, f0, s1, f1 = s0[i0], f0[i0], s1[i1], f1[i1]; g0 = common
        grows = (w0 >= 1e-3) & (w1 >= 1e-3)
        st, ft = 1e-3, -2.0
        S0, S1 = s0[:, grows] < st, s1[:, grows] < st
        F = f0[:, grows]; imp = F < ft; okf = np.isfinite(F)
        newly_ess = S1 & ~S0 & okf; newly_disp = S0 & ~S1 & okf
        row = {"org": org, "model": label, "genes_compared": len(common), "genes_dropped_by_rich_medium_step": ";".join(dropped),
               "wt_grows_before": int((w0 >= st).sum()), "wt_grows_after": int((w1 >= st).sum()),
               "pairs_changed_to_essential": int(newly_ess.sum()), "of_which_experimentally_important": int((newly_ess & imp).sum()),
               "pairs_changed_to_dispensable": int(newly_disp.sum()), "of_which_experimentally_important": int((newly_disp & imp).sum())}
        for name, fn in [("aucpr_bernstein", M.aucpr_bernstein), ("mcc", M.mcc), ("balanced_accuracy", M.balanced_accuracy)]:
            row[f"{name}_before"] = round(fn(s0[:, grows], F), 4); row[f"{name}_after"] = round(fn(s1[:, grows], F), 4)
        c0m, c1m = M.confusion(s0[:, grows], F), M.confusion(s1[:, grows], F)
        row.update({"fp_before": c0m["fp"], "fp_after": c1m["fp"], "fn_before": c0m["fn"], "fn_after": c1m["fn"]})
        rows.append(row)
        fb = load_organism(org); desc = dict(zip(fb.genes["sysName"], fb.genes["desc"]))
        for i, g in enumerate(g0):
            ne, nd = int(newly_ess[i].sum()), int(newly_disp[i].sum())
            if ne or nd:
                changes.append({"org": org, "model": label, "gene": g, "desc": desc.get(g, "")[:60],
                                "conditions_now_essential": ne, "of_which_important": int((newly_ess[i] & imp[i]).sum()),
                                "conditions_now_dispensable": nd, "of_which_important_": int((newly_disp[i] & imp[i]).sum()),
                                "mean_fitness_growth_conditions": float(np.nanmean(F[i]))})
    df = pd.DataFrame(rows); ch = pd.DataFrame(changes)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    print(df.to_string(index=False)); print(); print(ch.to_string(index=False) if len(ch) else "no gene-level changes")
    df.to_csv(os.path.join(RES, f"patch_{ver}_evaluation.tsv"), sep="\t", index=False)
    ch.to_csv(os.path.join(RES, f"patch_{ver}_gene_changes.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
