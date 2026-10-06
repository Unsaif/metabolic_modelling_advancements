"""How similar are the corrected draft's (U') and the curated model's calls on shared genes and conditions?
Agreement, Cohen's kappa, overlap of 'important' calls, and the composition of the shared errors, for the three
organisms compared in Paper 1. Usage: python -I call_agreement.py   (writes call_agreement.json)
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from common import ROOT  # noqa: E402
from recompute_comparison import FT, ST, grows, load, table  # noqa: E402

PAIRS = {
    "P. putida (iJN1463)": ("results/transfer_v1/development/Putida/UNQ/matrices.npz", "results/transfer_v1/reference_models/Putida/iJN1463/matrices.npz"),
    "S. oneidensis (iSO783)": ("results/transfer_v1/development/MR1/UNQ/matrices.npz", "results/transfer_v1/reference_models/MR1/iSO783/matrices.npz"),
    "B. thetaiotaomicron (iAH991, exploratory)": ("results/transfer_v1/development/Btheta/UNQ/matrices.npz",
                                                  "results/transfer_v1/reference_models/Btheta/iAH991_exploratory_B12_H2S/matrices.npz"),
}


def main():
    out = {}
    for name, (ua, ca) in PAIRS.items():
        u, c = load(os.path.join(ROOT, ua)), load(os.path.join(ROOT, ca))
        genes = sorted(set(u["browser_genes"]) & set(c["browser_genes"]))
        conds = sorted(grows(u) & grows(c))
        (su, sc), fit = table([u, c], genes, conds, absent_as_wt=False)
        ok = np.isfinite(su) & np.isfinite(sc) & np.isfinite(fit)
        pu, pc, obs = (su < ST) & ok, (sc < ST) & ok, (fit <= FT) & ok
        n = int(ok.sum())
        agree = int(((pu == pc) & ok).sum())
        po = agree / n
        pe = (pu.sum() / n) * (pc.sum() / n) + ((n - pu.sum()) / n) * ((n - pc.sum()) / n)
        kappa = (po - pe) / (1 - pe)
        wu, wc = ok & (pu != obs), ok & (pc != obs)
        out[name] = {"genes": len(genes), "conditions": len(conds), "cells": n, "agreement": po, "cohen_kappa": float(kappa),
                     "important_Uprime": int(pu.sum()), "important_curated": int(pc.sum()), "important_both": int((pu & pc).sum()),
                     "jaccard_important": float((pu & pc).sum() / (pu | pc).sum()),
                     "observed_important": int(obs.sum()),
                     "shared_errors": int((wu & wc).sum()), "shared_missed_important": int((wu & wc & obs).sum()),
                     "shared_false_important": int((wu & wc & ~obs).sum()),
                     "Uprime_errors": int(wu.sum()), "curated_errors": int(wc.sum()),
                     "expected_shared_if_independent": float(wu.sum() * wc.sum() / n),
                     "share_of_Uprime_errors_shared": float((wu & wc).sum() / wu.sum())}
        print(name, json.dumps(out[name], default=float))
    json.dump(out, open(os.path.join(HERE, "call_agreement.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
