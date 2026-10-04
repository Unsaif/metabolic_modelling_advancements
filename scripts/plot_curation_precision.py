"""Benchmark paper Figure 4: how often new "important gene" calls created by blind curation are confirmed.

For every evaluable organism and each curation arm (M_R1R2: removals and joins; M_R6: gene assignments; M: all),
compare with U′ (UNQ) on the union of both arms' genes and the conditions where both wild types grow (as in
scripts/transfer_compare.union_pairs). A new call is a (gene, condition) cell predicted to grow in U′ and not in
the curated arm; it is confirmed when the measured fitness is below the threshold (-2). For context the script
also reports, on the same cells, the precision of U′'s own "important" calls and the base rate of measured
importance. Curation removed no "important" call in any organism (reported as a check).
No simulation is run. Writes results/transfer_v1/curation_precision.json and docs/paper/fig4_curation_precision.png.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from gembench.comparison import load_run  # noqa: E402
from transfer_compare import union_pairs  # noqa: E402

PHASES = [("development", "Development", ["Btheta", "Putida", "MR1", "Smeli"]),
          ("evaluation_panel_A", "Panel A", ["Cola", "Dino", "Dyella79", "MycoTube", "RPal_CGA009"]),
          ("evaluation_panel_B", "Panel B", ["Caulo", "Cup4G11", "Marino", "PV4", "SB2B"])]
ARMS = [("M_R1R2", "Removals and joins (R1/R2)"), ("M_R6", "Gene assignments (R6)"), ("M", "All curation")]


def counts(role, org, arm):
    base = os.path.join(ROOT, "results", "transfer_v1", role, org)
    ra, rb = load_run(os.path.join(base, "UNQ")), load_run(os.path.join(base, arm))
    sa, sb, fit, genes, grows = union_pairs(ra, rb)
    st, ft = ra["params"]["growth_threshold"], ra["params"]["fitness_threshold"]
    ok = np.isfinite(fit) & np.isfinite(sa) & np.isfinite(sb)
    imp = ok & (fit < ft)
    new = ok & (sa >= st) & (sb < st)
    removed = ok & (sa < st) & (sb >= st)
    pred_u = ok & (sa < st)
    return {"new_calls": int(new.sum()), "new_confirmed": int((new & imp).sum()),
            "new_call_genes": int(new.any(axis=1).sum()), "removed_calls": int(removed.sum()),
            "u_calls": int(pred_u.sum()), "u_confirmed": int((pred_u & imp).sum()),
            "cells": int(ok.sum()), "important_cells": int(imp.sum())}


def main() -> None:
    table = []
    for role, label, orgs in PHASES:
        for org in orgs:
            for arm, name in ARMS:
                c = counts(role, org, arm)
                table.append({"phase": label, "org": org, "arm": arm, "kind": name, **c})
    out_json = os.path.join(ROOT, "results", "transfer_v1", "curation_precision.json")
    with open(out_json, "w") as fh:
        json.dump(table, fh, indent=1)
        fh.write("\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9, 5.2), sharey=True)
    orgs = [(label, org) for _, label, olist in PHASES for org in olist]
    ypos = {o: len(orgs) - 1 - i for i, o in enumerate(orgs)}
    colours = {"Development": "#7f7f7f", "Panel A": "#1f77b4", "Panel B": "#d62728"}
    for ax, (arm, name) in zip(axes, ARMS[:2]):
        for row in (t for t in table if t["arm"] == arm):
            y = ypos[(row["phase"], row["org"])]
            n, k = row["new_calls"], row["new_confirmed"]
            ax.barh(y, k, color=colours[row["phase"]], alpha=0.85)
            ax.barh(y, n - k, left=k, color="white", edgecolor=colours[row["phase"]], hatch="///", lw=0.8)
            if n:
                ax.text(n + 1, y, f"{k}/{n}", va="center", fontsize=7)
        tot_n = sum(t["new_calls"] for t in table if t["arm"] == arm and t["phase"] != "Development")
        tot_k = sum(t["new_confirmed"] for t in table if t["arm"] == arm and t["phase"] != "Development")
        no_tb = sum(t["new_calls"] for t in table if t["arm"] == arm and t["phase"] != "Development" and t["org"] != "MycoTube")
        ax.set_title(f"{name}\nnew panels: {tot_k} of {tot_n} new calls confirmed ({tot_k} of {no_tb} without M. tuberculosis)",
                     fontsize=8.5)
        ax.set_xlabel("New 'important gene' calls (gene x condition cells)\nfilled: confirmed by fitness < -2; hatched: not confirmed", fontsize=8)
        ax.tick_params(axis="x", labelsize=8)
    axes[0].set_yticks(list(ypos.values()))
    axes[0].set_yticklabels([f"{org} ({label})" for label, org in ypos], fontsize=8)
    xmax = max(t["new_calls"] for t in table if t["arm"] in ("M_R1R2", "M_R6")) * 1.25
    for ax in axes:
        ax.set_xlim(0, xmax)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "docs", "paper", "fig4_curation_precision.png"), dpi=200)
    for row in table:
        base = row["important_cells"] / row["cells"] if row["cells"] else float("nan")
        up = row["u_confirmed"] / row["u_calls"] if row["u_calls"] else float("nan")
        print(f"{row['phase']:12s} {row['org']:12s} {row['arm']:7s} new {row['new_confirmed']:3d}/{row['new_calls']:3d} "
              f"removed {row['removed_calls']} | U' precision {up:.2f} base rate {base:.3f}")


if __name__ == "__main__":
    main()
