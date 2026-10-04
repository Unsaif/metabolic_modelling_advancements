"""Attribution of the transfer-study gain to its components (benchmark paper, Figure 3).

For each component comparison and each phase (development organisms, panel A, panel B) the script takes the
per-organism union-gene paired MCC differences from results/transfer_v1/<phase>/paired.json and summarises them
with scripts/transfer_aggregate.aggregate (unweighted mean over evaluable organisms, 10,000-resample organism
bootstrap, seed 0), exactly as in the results documents. No new simulation is run.

Writes results/transfer_v1/attribution_by_component.json and docs/paper/fig3_attribution.png.
"""
from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from transfer_aggregate import aggregate  # noqa: E402

PHASES = [("development", "Development (4 organisms)"), ("evaluation_panel_A", "Panel A (new)"),
          ("evaluation_panel_B", "Panel B (replication)")]
COMPONENTS = [
    ("B0", "Uonly_rxn", "Universal reaction patches"),
    ("B0", "Uonly_atp", "ATP synthase if annotated"),
    ("U", "UN", "Rule normalisation"),
    ("U", "UQ", "Menaquinone rule"),
    ("B0", "UNQ", "All automatic rules (H1)"),
    ("UNQ", "M_R1R2", "Curation: removals and joins (R1/R2)"),
    ("UNQ", "M_R6", "Curation: gene assignments (R6)"),
    ("UNQ", "M", "All blind curation (H2)"),
]


def main() -> None:
    table = []
    for phase, label in PHASES:
        rows = json.load(open(os.path.join(ROOT, "results", "transfer_v1", phase, "paired.json")))
        agg = {(r["A"], r["B"]): r for r in aggregate(rows)}
        for a, b, name in COMPONENTS:
            r = agg[(a, b)]
            table.append({"phase": phase, "component": name, "A": a, "B": b, "n_evaluable": r["n_evaluable"],
                          "mean": r["mean_delta"], "ci95": r["organism_bootstrap_ci95"],
                          "up_down_unchanged": [r["n_improved"], r["n_worsened"], r["n_unchanged"]],
                          "per_organism": {o["org"]: o["delta"] for o in r["per_organism"] if o["delta"] is not None}})
    out_json = os.path.join(ROOT, "results", "transfer_v1", "attribution_by_component.json")
    with open(out_json, "w") as fh:
        json.dump(table, fh, indent=1)
        fh.write("\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colours = {"development": "#7f7f7f", "evaluation_panel_A": "#1f77b4", "evaluation_panel_B": "#d62728"}
    offsets = {"development": 0.25, "evaluation_panel_A": 0.0, "evaluation_panel_B": -0.25}
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    n = len(COMPONENTS)
    for i, (_, _, name) in enumerate(COMPONENTS):
        y0 = n - 1 - i
        if name.startswith("All"):
            ax.axhspan(y0 - 0.45, y0 + 0.45, color="#f2f2f2", zorder=0)
        for phase, label in PHASES:
            row = next(t for t in table if t["phase"] == phase and t["component"] == name)
            y = y0 + offsets[phase]
            lo, hi = row["ci95"]
            ax.plot([lo, hi], [y, y], color=colours[phase], lw=2, zorder=2)
            ax.plot(row["mean"], y, "o", color=colours[phase], ms=6, zorder=3,
                    label=label if i == 0 else None)
            ax.plot(list(row["per_organism"].values()), [y] * len(row["per_organism"]), "|",
                    color=colours[phase], ms=7, alpha=0.5, zorder=2)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(range(n))
    ax.set_yticklabels([c[2] for c in COMPONENTS][::-1], fontsize=9)
    ax.set_xlabel("Paired MCC difference, union of genes\n(dot: mean over organisms; bar: 95% organism bootstrap; ticks: organisms)", fontsize=9)
    ax.set_title("Components of the corrections, by phase", fontsize=10)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False)
    ax.tick_params(axis="x", labelsize=8)
    fig.tight_layout()
    out_png = os.path.join(ROOT, "docs", "paper", "fig3_attribution.png")
    fig.savefig(out_png, dpi=200)
    for t in table:
        print(f"{t['phase']:20s} {t['component']:40s} {t['mean']:+.3f} [{t['ci95'][0]:+.3f}, {t['ci95'][1]:+.3f}] "
              f"{t['up_down_unchanged']}")


if __name__ == "__main__":
    main()
