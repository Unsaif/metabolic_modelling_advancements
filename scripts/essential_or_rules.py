"""Which OR-rule reactions actually matter for the benchmark: reactions with an OR in their gene rule that
are essential (removal abolishes growth) in at least one carbon-source condition where the wild type
grows. Adjudicating gene rules is only worth doing where the rule can change a prediction.
Writes results/carbon_fitness_multi/or_rules_conditionally_essential_<org>.tsv with the genes' annotations.
"""
from __future__ import annotations
import glob, gzip, json, logging, os, re, sys
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
logging.getLogger("cobra").setLevel(logging.ERROR)
import cobra
from cobra.flux_analysis import single_reaction_deletion
from gembench.fitness_browser import base_medium, carbon_source_conditions, load_organism
from gembench.gene_mapping import accession_from_model_gene, load_genpept_table
from gembench.media import apply_medium

def apply_patches(m, org):
    up = json.load(open(os.path.join(ROOT, "data/reference/universe_patches_v0.1.json")))
    for pt in up["patches"]:
        if pt["reaction"] in m.reactions and "gpr" not in pt["change"] and "CBPS" in m.reactions:
            m.reactions.get_by_id(pt["reaction"]).bounds = (pt["change"]["lower_bound"], pt["change"]["upper_bound"])
    return m

def main():
    for org in (sys.argv[1].split(",") if len(sys.argv) > 1 else ["Btheta", "Putida", "MR1", "Smeli"]):
        with gzip.open(os.path.join(ROOT, "models", "gapfilled", f"{org}.xml.gz"), "rt") as fh:
            m = cobra.io.read_sbml_model(fh)
        m.solver = "glpk"; apply_patches(m, org)
        fb = load_organism(org); conds = [c for c in carbon_source_conditions(fb) if c.bigg_ids]
        gp = load_genpept_table(org).set_index("version")
        or_rxns = [r.id for r in m.reactions if " or " in r.gene_reaction_rule]
        ess = {}
        for c in conds:
            with m:
                apply_medium(m, base_medium(c.media))
                for ex in c.exchanges:
                    if ex in m.reactions: m.reactions.get_by_id(ex).lower_bound = -10.0
                wt = m.slim_optimize()
                if not wt or wt < 1e-3: continue
                res = single_reaction_deletion(m, or_rxns, processes=2)
                for ids, g in zip(res["ids"], res["growth"]):
                    rid = next(iter(ids))
                    if g is None or np.isnan(g) or g < 1e-3:
                        ess.setdefault(rid, []).append(c.name)
        rows = []
        for rid, cl in ess.items():
            r = m.reactions.get_by_id(rid)
            genes = sorted(set(re.findall(r"[A-Z]P_\d+_\d|BT\d+|PP_\d+|SO_?\d+|SM[a-c_]\w+", r.gene_reaction_rule)))
            ann = []
            for g in genes:
                acc = accession_from_model_gene(g)
                lt = gp.loc[acc, "locus_tags"].replace("_", "") if acc in gp.index else g
                d = gp.loc[acc, "definition"].split(" [")[0] if acc in gp.index else "?"
                ann.append(f"{lt}: {d}")
            rows.append({"reaction": rid, "name": r.name, "rule": r.gene_reaction_rule, "n_conditions_essential": len(cl),
                         "conditions": ";".join(cl[:5]), "genes": " || ".join(ann)})
        df = pd.DataFrame(rows).sort_values("n_conditions_essential", ascending=False)
        df.to_csv(os.path.join(ROOT, "results", "carbon_fitness_multi", f"or_rules_conditionally_essential_{org}.tsv"), sep="\t", index=False)
        print(f"{org}: {len(or_rxns)} OR-rule reactions, {len(df)} essential in >=1 growth condition", flush=True)

if __name__ == "__main__":
    main()
