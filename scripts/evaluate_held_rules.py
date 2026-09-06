"""Re-test held gene-rule patches one at a time on top of the current best arm.

For every patch with status 'held' in the given gene-rule files, the model is built as in the best arm
(gap-filled draft + universe patches + accepted gene-rule patches + medium completion), the held rule is applied
in isolation, and only the genes named in the old or the new rule are knocked out across the carbon-source
conditions. Each gene-condition prediction that changes is compared with the fitness data (important =
fitness below the threshold), which is the same attribution the verifier used to hold the patch, now under
medium completion and the later universe patches.

Usage: python scripts/evaluate_held_rules.py [--orgs Btheta,Putida,MR1,Smeli] [--variant gapfilled]
           [--patch data/reference/universe_patches_v0.1.json]
           [--gpr-patch data/reference/gpr_patches_v0.2.json,data/reference/gpr_patches_v0.3.json]
Output: results/carbon_fitness_multi/held_rules_reevaluation.tsv (one row per held patch) and
        held_rules_reevaluation_genes.tsv (one row per gene whose predictions change).
"""
from __future__ import annotations

import argparse
import logging
import os
import re
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
logging.getLogger("cobra").setLevel(logging.ERROR)

from gembench.fitness_browser import carbon_source_conditions, load_organism  # noqa: E402
from gembench.patches import apply_gpr_patches, apply_model_patches, apply_universe_patches, load_patch_files, translate_rule  # noqa: E402
from gembench.protocols import carbon_fitness_generic as P  # noqa: E402
from run_carbon_fitness_generic import gene_map_for, load_model  # noqa: E402

OUT = os.path.join(ROOT, "results", "carbon_fitness_multi")
_TOKEN = re.compile(r"[A-Za-z][A-Za-z0-9_]*")


def main() -> None:
    import json
    ap = argparse.ArgumentParser()
    ap.add_argument("--orgs", default="Btheta,Putida,MR1,Smeli")
    ap.add_argument("--variant", default="gapfilled")
    ap.add_argument("--patch", default=os.path.join(ROOT, "data/reference/universe_patches_v0.1.json"))
    ap.add_argument("--gpr-patch", default=",".join(os.path.join(ROOT, f) for f in ["data/reference/gpr_patches_v0.2.json", "data/reference/gpr_patches_v0.3.json"]))
    ap.add_argument("--model-patch", default=os.path.join(ROOT, "data/reference/model_patches_v0.1.json"))
    ap.add_argument("--no-medium-completion", action="store_true")
    ap.add_argument("--processes", type=int, default=1)
    args = ap.parse_args()
    patches = json.load(open(args.patch))
    model_patches = json.load(open(args.model_patch)) if args.model_patch else None
    gpr = load_patch_files(args.gpr_patch.split(","))

    rows, gene_rows = [], []
    for org in args.orgs.split(","):
        fb = load_organism(org)
        model, _, _ = load_model(org, args.variant)
        apply_universe_patches(model, org, patches)
        gm = gene_map_for(org, model, set(fb.genes["sysName"]), args.variant)
        if model_patches and apply_model_patches(model, org, model_patches, gm, verbose=False):
            gm = gene_map_for(org, model, set(fb.genes["sysName"]), args.variant)
        acc = apply_gpr_patches(model, org, gm, gpr, statuses=("accepted",), verbose=False)
        if any(a.get("genes_added") for a in acc):
            gm = gene_map_for(org, model, set(fb.genes["sysName"]), args.variant)
        conds = carbon_source_conditions(fb)
        model_ids = {g.id for g in model.genes}
        fit_index = set(fb.fitness.index)
        held = [pt for pt in gpr["patches"] if pt["org"] == org and pt.get("status") == "held" and pt["reaction"] in model.reactions]
        print(f"== {org}: {len(held)} held patches to re-test on {model.id}", flush=True)
        for pt in held:
            r = model.reactions.get_by_id(pt["reaction"])
            old_rule = r.gene_reaction_rule
            new_rule = translate_rule(pt["new_rule"], gm, model_ids)
            toks = {t for t in _TOKEN.findall(old_rule + " " + new_rule) if t not in ("and", "or")}
            genes = sorted(t for t in toks if t in model_ids and gm.model_to_browser.get(t) in fit_index)
            if not genes:
                rows.append({"org": org, "reaction": pt["reaction"], "rule": pt.get("rule"), "new_rule": pt["new_rule"], "genes_tested": 0,
                             "note": "no gene of the old or new rule has fitness data"})
                continue
            params = P.GenericParams(drop_rich_medium_essentials=False, processes=args.processes, genes_subset=genes,
                                     complete_medium_transport=not args.no_medium_completion)
            base = P.run(model, fb, conds, gm, params, verbose=False)
            r.gene_reaction_rule = new_rule
            try:
                alt = P.run(model, fb, conds, gm, params, verbose=False)
            finally:
                r.gene_reaction_rule = old_rule
            thr, ft = params.growth_threshold, params.fitness_threshold
            grows = base.wt_growth >= thr
            eb, ea = base.sim_growth < thr, alt.sim_growth < thr
            to_ess = (~eb & ea) & grows[None, :]
            to_dis = (eb & ~ea) & grows[None, :]
            imp = base.fitness < ft
            row = {"org": org, "reaction": pt["reaction"], "rule": pt.get("rule"), "new_rule": pt["new_rule"], "old_rule_model_ids": old_rule,
                   "genes_tested": len(genes), "growth_conditions": int(grows.sum()),
                   "to_essential": int(to_ess.sum()), "to_essential_important": int((to_ess & imp).sum()),
                   "to_dispensable": int(to_dis.sum()), "to_dispensable_important": int((to_dis & imp).sum()),
                   "held_reason": pt.get("verifier", "")}
            # verdict on the same criterion as the verifier: changes must be supported by the fitness data
            n_ch = row["to_essential"] + row["to_dispensable"]
            n_ok = row["to_essential_important"] + (row["to_dispensable"] - row["to_dispensable_important"])
            row["changes_supported"] = n_ok
            row["verdict"] = "neutral" if n_ch == 0 else ("supported" if n_ok / n_ch >= 0.8 else ("contradicted" if n_ok / n_ch <= 0.2 else "mixed"))
            rows.append(row)
            for i, g in enumerate(base.model_genes):
                te, td = to_ess[i], to_dis[i]
                if te.any() or td.any():
                    gene_rows.append({"org": org, "reaction": pt["reaction"], "model_gene": g, "browser_gene": base.browser_genes[i],
                                      "to_essential": int(te.sum()), "to_essential_important": int((te & imp[i]).sum()),
                                      "to_dispensable": int(td.sum()), "to_dispensable_important": int((td & imp[i]).sum()),
                                      "mean_fitness_growth_conditions": float(np.nanmean(base.fitness[i][grows])) if grows.any() else float("nan"),
                                      "conditions_to_essential": ";".join(c.name for c, x in zip(base.conditions, te) if x)})
            print(f"   {pt['reaction']:10s} {pt.get('rule')}: genes={len(genes)} to_essential={row['to_essential']} "
                  f"(important {row['to_essential_important']}) to_dispensable={row['to_dispensable']} "
                  f"(important {row['to_dispensable_important']}) -> {row['verdict']}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "held_rules_reevaluation.tsv"), sep="\t", index=False)
    pd.DataFrame(gene_rows).to_csv(os.path.join(OUT, "held_rules_reevaluation_genes.tsv"), sep="\t", index=False)
    print(df[["org", "reaction", "rule", "genes_tested", "to_essential", "to_essential_important", "to_dispensable", "to_dispensable_important", "verdict"]].to_string(index=False))


if __name__ == "__main__":
    main()
