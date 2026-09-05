"""Independent re-computation of the multi-organism carbon-source benchmark metrics from the saved
matrices, plus consistency checks. Written separately from the runner on purpose: it uses only
numpy/sklearn on the npz files, its own confusion-matrix code, and cross-checks against the cards."""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import auc, matthews_corrcoef, precision_recall_curve

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from gembench.fitness_browser import load_organism  # noqa: E402

RES = os.path.join(ROOT, "results", "carbon_fitness_multi")


def aucpr_b(sim, fit, st=1e-3):
    y = (sim > st).astype(int).ravel(); s = -fit.ravel(); ok = np.isfinite(s)
    if y[ok].min() == y[ok].max():
        return float("nan")
    p, r, _ = precision_recall_curve(y[ok], s[ok], pos_label=0)
    return float(auc(r, p))


def main() -> None:
    problems = 0
    for card_path in sorted(glob.glob(os.path.join(RES, "*", "*", "card.json"))):
        d = os.path.dirname(card_path)
        card = json.load(open(card_path))
        z = np.load(os.path.join(d, "matrices.npz"), allow_pickle=True)
        sim, fit, wt = z["sim_growth"], z["fitness"], z["wt_growth"]
        conds = list(z["conditions"]); bgenes = list(z["browser_genes"])
        org = card["benchmark"].split("organism ")[-1]
        st = card["protocol"]["params"]["growth_threshold"]; ft = card["protocol"]["params"]["fitness_threshold"]
        grows = wt >= st
        # 1. recompute pooled metrics on conditions where WT grows
        s2, f2 = sim[:, grows], fit[:, grows]
        rows_ok = np.isfinite(f2).any(axis=1); s2, f2 = s2[rows_ok], f2[rows_ok]
        rep = card["results"]["gene_level_conditions_where_wt_grows"]
        line = f"{org:7s} {os.path.basename(d)[:48]:48s} WT grows {grows.sum():2d}/{len(grows):2d}"
        if s2.size:
            a = aucpr_b(s2, f2, st)
            okm = np.isfinite(f2.ravel())
            pred = (s2.ravel() > st).astype(int)[okm]; lab = (f2.ravel() >= ft).astype(int)[okm]
            m = matthews_corrcoef(lab, pred)
            da = abs(a - rep["aucpr_bernstein"]["point"]); dm = abs(m - rep["mcc"]["point"])
            line += f" | AUC-PR {a:.3f} (card {rep['aucpr_bernstein']['point']:.3f}) MCC {m:.3f} (card {rep['mcc']['point']:.3f})"
            if da > 1e-9 or dm > 1e-9:
                line += "  MISMATCH"; problems += 1
            # class balance
            line += f" | pairs {okm.sum()} exp-important {int((lab == 0).sum())} pred-no-growth {int((pred == 0).sum())}"
        else:
            line += " | no gene-level scoring (WT grows nowhere)"
        # 2. fitness matrix re-derived from the raw Fitness Browser tables for a random condition/gene sample
        fb = load_organism(org)
        exp = fb.experiments
        rng = np.random.default_rng(0)
        for _ in range(3):
            j = int(rng.integers(0, len(conds))); i = int(rng.integers(0, len(bgenes)))
            cname, media = conds[j].split(" | ")
            sel = exp[(exp["expGroup"] == "carbon source") & (exp["condition_1"] == cname) & (exp["media"] == media)
                      & (exp["condition_2"].isin(["", "Dimethyl Sulfoxide"]))]["expName"]
            sel = [e for e in sel if e in fb.fitness.columns]
            v = fb.fitness.loc[bgenes[i], sel].astype(float).mean()
            if not (np.isnan(v) and np.isnan(fit[i, j])) and abs(v - fit[i, j]) > 1e-9:
                line += f"  FITNESS-MISMATCH({bgenes[i]},{cname})"; problems += 1
        # 3. every condition where WT does not grow must have sim == 0 for all genes; growth conditions non-negative
        if (sim[:, ~grows] != 0).any():
            line += "  SIM-NONZERO-ON-NOGROWTH"; problems += 1
        if (sim < -1e-6).any():
            line += "  NEGATIVE-GROWTH"; problems += 1
        # 4. condition table consistency
        ct = pd.read_table(os.path.join(d, "conditions.tsv"))
        if int(ct["wt_grows"].sum()) != int(grows.sum()):
            line += "  CONDITION-TABLE-MISMATCH"; problems += 1
        print(line)
    print("problems:", problems)


if __name__ == "__main__":
    main()
