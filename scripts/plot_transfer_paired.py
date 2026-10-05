"""Figure 2 of the benchmark paper: per-organism paired MCC differences for H1 (U' vs B0) and H2 (M vs U').

Reads only existing results: results/transfer_v1/<role>/paired.json (per-organism differences and gene-bootstrap
intervals) and aggregate_union.json (organism means and organism-bootstrap intervals), plus the pooled file for
panels A and B. Colours match Figures 1, 3 and 4 (development grey, panel A blue, panel B red).

Usage: python scripts/plot_transfer_paired.py   -> docs/paper/fig2_paired_differences.png
"""
from __future__ import annotations

import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "transfer_v1")
PHASES = [("development", "Development", ["MR1", "Btheta", "Smeli", "Putida"]),
          ("evaluation_panel_A", "Panel A", ["Cola", "Dino", "Dyella79", "MycoTube", "RPal_CGA009"]),
          ("evaluation_panel_B", "Panel B", ["Caulo", "Cup4G11", "Marino", "PV4", "SB2B"])]
NAMES = {"MR1": "S. oneidensis MR-1", "Btheta": "B. thetaiotaomicron", "Smeli": "S. meliloti", "Putida": "P. putida",
         "Cola": "E. vietnamensis", "Dino": "D. shibae", "Dyella79": "D. japonica", "MycoTube": "M. tuberculosis",
         "RPal_CGA009": "R. palustris", "Caulo": "C. crescentus", "Cup4G11": "C. basilensis", "Marino": "M. adhaerens",
         "PV4": "S. loihica", "SB2B": "S. amazonensis"}
COLOURS = {"Development": "#7f7f7f", "Panel A": "#1f77b4", "Panel B": "#d62728"}
PAIRS = [(("B0", "UNQ"), "H1: automatic rules (U′ vs B0)"), (("UNQ", "M"), "H2: blind AI curation (M vs U′)")]


def load(path):
    with open(path) as fh:
        return json.load(fh)


def main() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = []      # (label, phase, kind, {pair: (delta, lo, hi)})
    for role, phase, orgs in PHASES:
        paired = load(os.path.join(R, role, "paired.json"))
        agg = load(os.path.join(R, role, "aggregate_union.json"))
        for org in orgs:
            vals = {}
            for (a, b), _ in PAIRS:
                hit = [p for p in paired if p["org"] == org and p["A"] == a and p["B"] == b]
                if hit:
                    vals[(a, b)] = (hit[0]["union_delta"], *hit[0]["union_delta_ci95"])
            rows.append((NAMES[org], phase, "organism", vals))
        means = {}
        for (a, b), _ in PAIRS:
            p = [q for q in agg["pairs"] if q["A"] == a and q["B"] == b][0]
            means[(a, b)] = (p["mean_delta"], *p["organism_bootstrap_ci95"])
        rows.append((f"{phase}: mean of {len(orgs)}", phase, "mean", means))
    pooled = load(os.path.join(R, "pooled_panels_A_B_aggregate_union.json"))
    pm = {}
    for (a, b), _ in PAIRS:
        p = [q for q in pooled["pairs"] if q["A"] == a and q["B"] == b][0]
        pm[(a, b)] = (p["mean_delta"], *p["organism_bootstrap_ci95"])
    rows.append(("Panels A + B: mean of 10", "Pooled", "pooled", pm))

    fig, axes = plt.subplots(1, 2, figsize=(10, 7.2), sharey=True)
    n = len(rows)
    ypos = []
    y = 0.0
    for i, (_, phase, kind, _) in enumerate(rows):
        if i and rows[i - 1][1] != phase:
            y += 0.6
        ypos.append(y)
        y += 1
    for ax, (pair, title) in zip(axes, PAIRS):
        for (label, phase, kind, vals), yy in zip(rows, ypos):
            if pair not in vals:
                continue
            d, lo, hi = vals[pair]
            colour = "black" if kind == "pooled" else COLOURS[phase]
            marker = "o" if kind == "organism" else "D"
            size = 5 if kind == "organism" else 7
            ax.plot([lo, hi], [yy, yy], color=colour, lw=1.4 if kind == "organism" else 2.2,
                    alpha=0.85 if kind == "organism" else 1)
            ax.plot([d], [yy], marker=marker, color=colour, ms=size)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Paired MCC difference, union of genes\n(organisms: 95% gene bootstrap; means: 95% organism bootstrap)",
                      fontsize=8.5)
        ax.tick_params(axis="x", labelsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_yticks(ypos)
    labels = []
    for label, phase, kind, _ in rows:
        labels.append(label if kind != "organism" else label)
    axes[0].set_yticklabels(labels, fontsize=8.5)
    for tick, (label, phase, kind, _) in zip(axes[0].get_yticklabels(), rows):
        if kind == "organism":
            tick.set_fontstyle("italic")
        else:
            tick.set_fontweight("bold")
        tick.set_color("black" if kind == "pooled" else COLOURS[phase])
    axes[0].invert_yaxis()
    fig.text(0.01, 0.01, "Not evaluable under the frozen rules: D. suillum (panel A, gap-fill), "
             "D. vulgaris Miyazaki F (panel B, reference condition). Development organisms are retrospective.",
             fontsize=7.5, color="#555555")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out = os.path.join(ROOT, "docs", "paper", "fig2_paired_differences.png")
    fig.savefig(out, dpi=200)
    print(out)


if __name__ == "__main__":
    main()
