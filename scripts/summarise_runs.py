"""Tabulate every carbon-source fitness run (all organisms, variants, patch arms and scoring conventions).

Reads results/carbon_fitness_multi/<org>/<run>/card.json and writes results/carbon_fitness_multi/summary_all_runs.tsv,
one row per run, with the gene-level metrics pooled over the conditions where the wild-type model grows.
The 'arm' column encodes what was applied: variant, universe patches, gene-rule patch versions, medium completion,
and whether rich-medium essentials were dropped (Bernstein convention) or every gene with fitness data was scored.
"""
from __future__ import annotations

import glob
import json
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results", "carbon_fitness_multi")


def main() -> None:
    rows = []
    for card_path in sorted(glob.glob(os.path.join(RES, "*", "*", "card.json"))):
        org = os.path.basename(os.path.dirname(os.path.dirname(card_path)))
        run = os.path.basename(os.path.dirname(card_path))
        c = json.load(open(card_path))
        pr = c["protocol"]
        pt = pr.get("patches") or {}
        params = pr.get("params", {})
        r = c["results"]
        g = r["gene_level_conditions_where_wt_grows"]
        cl = r["condition_level"]
        applied = pt.get("applied") or []
        rows.append({
            "org": org, "model": c["model"]["model_id"], "variant": pr.get("variant"),
            "universe_patches": os.path.basename(pt["file"]) if pt.get("file") else "",
            "gpr_patches": ",".join(os.path.basename(f) for f in pt["gpr_file"].split(",")) if pt.get("gpr_file") else "",
            "reaction_patches_applied": ";".join(a["reaction"] for a in applied if "rule" not in a),
            "gene_rule_patches_applied": sum(1 for a in applied if "rule" in a),
            "medium_completion": bool(params.get("complete_medium_transport", False)),
            "scoring": "drop_rich_essentials" if params.get("drop_rich_medium_essentials", True) else "all_genes",
            "genes_scored": g.get("n_genes", 0), "cond_mapped": cl["n_conditions_mapped"], "cond_wt_grows": cl["n_conditions_wt_grows"],
            "pairs": g.get("n_gene_condition_pairs", 0),
            "aucpr_bernstein": _pt(g, "aucpr_bernstein"), "mcc": _pt(g, "mcc"), "mcc_ci_low": _ci(g, "mcc", 0), "mcc_ci_high": _ci(g, "mcc", 1),
            "balanced_accuracy": _pt(g, "balanced_accuracy"), "auroc": _pt(g, "auroc_standard"),
            "tp": g.get("confusion", {}).get("tp"), "tn": g.get("confusion", {}).get("tn"),
            "fp": g.get("confusion", {}).get("fp"), "fn": g.get("confusion", {}).get("fn"),
            "created": c.get("created", ""), "run_dir": os.path.relpath(os.path.dirname(card_path), ROOT),
        })
    df = pd.DataFrame(rows)
    out = os.path.join(RES, "summary_all_runs.tsv")
    df.to_csv(out, sep="\t", index=False, float_format="%.4f")
    print(df[["org", "variant", "reaction_patches_applied", "gene_rule_patches_applied", "medium_completion", "scoring",
              "cond_wt_grows", "cond_mapped", "genes_scored", "mcc", "balanced_accuracy", "aucpr_bernstein"]].to_string(index=False))
    print(f"\nwritten {out} ({len(df)} runs)")


def _pt(g, k):
    v = g.get(k)
    return v.get("point") if isinstance(v, dict) else v


def _ci(g, k, i):
    v = g.get(k)
    return v.get("ci95", [None, None])[i] if isinstance(v, dict) else None


if __name__ == "__main__":
    main()
