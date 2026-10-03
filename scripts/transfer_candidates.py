"""Build blind adjudication packets for transfer study v1 (no fitness data are read).

For one organism and one arm (normally U), the model is prepared exactly as for scoring (arm transforms, medium
completion), the carbon-source conditions are taken from experiment metadata, and two kinds of candidate reactions are
listed: gene-less reactions and OR-rule reactions that are essential (deletion abolishes growth) in at least one
condition where the wild-type model grows. Each candidate carries its equation, its genes' annotations and their
genomic neighbourhood. The packet directory also receives the organism's full gene table and the frozen procedure.

Usage: python scripts/transfer_candidates.py --org X --arm U --out <packet-dir>
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
logging.getLogger("cobra").setLevel(logging.ERROR)

from cobra.flux_analysis import single_reaction_deletion  # noqa: E402

from gembench import fitness_browser as FB  # noqa: E402
from gembench import transfer as T  # noqa: E402
from gembench.gene_mapping import accession_from_model_gene, build_gene_map, load_genpept_table  # noqa: E402
from gembench.media import apply_medium  # noqa: E402
from gembench.protocols import carbon_fitness_generic as P  # noqa: E402
from run_transfer_study import MODELS, apply_arm, genes_table, load_config, read_sbml, reference_tables  # noqa: E402

PROCEDURE = os.path.join(ROOT, "docs", "studies", "transfer-adjudication-procedure-v1.md")
SKIP_PREFIX = ("EX_", "DM_", "SK_", "sink_", "MEDt_")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--org", required=True)
    ap.add_argument("--arm", default="U")
    ap.add_argument("--out", required=True)
    ap.add_argument("--processes", type=int, default=2)
    a = ap.parse_args()
    cfg = load_config(a.org)
    arms = json.load(open(os.path.join(ROOT, "data", "studies", "transfer_v1", "arms.json")))["arms"]
    genes = genes_table(cfg)
    sysnames = set(genes["sysName"])
    genpept = os.path.join(ROOT, cfg["genpept_map"])
    exps = pd.read_table(os.path.join(ROOT, cfg["fb_dir"], "experiments.tsv"), dtype=str, keep_default_na=False)
    model = read_sbml(os.path.join(MODELS, f"{a.org}_base.xml.gz"))
    model.solver = "glpk"

    def gm_builder(m):
        return build_gene_map(a.org, [g.id for g in m.genes], sysnames, genpept_path=genpept)

    with reference_tables(cfg) as (media_table, carbon_table):
        model, applied, gm = apply_arm(model, cfg, arms[a.arm], gm_builder)
        conds = [c for c in T.conditions_from_metadata(exps, FB.load_carbon_source_map(), FB.load_media_map()) if c["bigg_ids"]]
        params = P.GenericParams()
        P.complete_medium_transport(model, sorted({c["media"] for c in conds}), params.medium_completion_exclude)
        geneless = [r.id for r in model.reactions if not r.genes and not r.id.startswith(SKIP_PREFIX) and r.id != "Growth"]
        orrule = [r.id for r in model.reactions if " or " in r.gene_reaction_rule]
        ess = {}
        n_grow = 0
        for ex in model.exchanges:                      # as in the scoring protocol
            ex.lower_bound, ex.upper_bound = 0.0, 1000.0
        for c in conds:
            # as in the scoring protocol: open the carbon exchanges the model has; a condition is simulated even if
            # some (or all) of its carbon exchanges are absent
            exs = [f"EX_{b}_e" for b in c["bigg_ids"] if f"EX_{b}_e" in model.reactions]
            with model:
                apply_medium(model, FB.base_medium(c["media"]), close_all=True)
                for x in exs:
                    model.reactions.get_by_id(x).lower_bound = -10.0
                wt = model.slim_optimize()
                if wt is None or not np.isfinite(wt) or wt < 1e-3:
                    continue
                n_grow += 1
                res = single_reaction_deletion(model, geneless + orrule, processes=a.processes)
                for ids, g, st in zip(res["ids"], res["growth"], res["status"]):
                    rid = next(iter(ids))
                    if st == "infeasible" or (g is not None and np.isfinite(g) and g < 1e-3):
                        ess.setdefault(rid, []).append(c["key"])
    gp = load_genpept_table(a.org, genpept).set_index("version")
    order = genes.assign(_b=genes["begin"].astype(int)).sort_values(["scaffoldId", "_b"]).reset_index(drop=True)
    pos = {s: i for i, s in enumerate(order["sysName"])}

    def gene_info(gid):
        acc = accession_from_model_gene(gid)
        sysn = gm.model_to_browser.get(gid)
        d = {"model_gene": gid, "locus": sysn, "fitness_browser_desc": T.clean_desc(genes.set_index("sysName").at[sysn, "desc"]) if sysn else None,
             "refseq_definition": gp.at[acc, "definition"].split(" [")[0] if acc in gp.index else None}
        if sysn in pos:
            i = pos[sysn]
            nb = order.iloc[max(0, i - 3): i + 4]
            d["neighbourhood"] = [f"{r.sysName} ({r.strand}) {T.clean_desc(r.desc)}" for r in nb.itertuples() if r.scaffoldId == order.at[i, "scaffoldId"]]
        return d

    cands = []
    for rid, cl in sorted(ess.items()):
        r = model.reactions.get_by_id(rid)
        eq = r.build_reaction_string(use_metabolite_names=True)
        ann = {k: v for k, v in r.annotation.items() if k in ("ec-code", "kegg.reaction", "metanetx.reaction", "rhea", "bigg.reaction")}
        cands.append({"reaction": rid, "name": r.name, "equation_ids": r.reaction, "equation_names": eq,
                      "annotation": ann, "type": "A_geneless" if not r.genes else "B_or_rule",
                      "current_rule": r.gene_reaction_rule, "n_growth_conditions_essential": len(cl),
                      "n_growth_conditions": n_grow, "genes": [gene_info(g.id) for g in sorted(r.genes, key=lambda x: x.id)]})
    os.makedirs(a.out, exist_ok=True)
    json.dump({"org": a.org, "arm": a.arm, "n_growth_conditions_simulated": n_grow, "candidates": cands}, open(os.path.join(a.out, "candidates.json"), "w"), indent=1)
    genes[["locusId", "sysName", "scaffoldId", "begin", "end", "strand", "desc"]].to_csv(os.path.join(a.out, "genes.tsv"), sep="\t", index=False)
    shutil.copy(PROCEDURE, os.path.join(a.out, "PROCEDURE.md"))
    print(f"{a.org}: {n_grow} growth conditions; {sum(c['type']=='A_geneless' for c in cands)} gene-less and "
          f"{sum(c['type']=='B_or_rule' for c in cands)} OR-rule essential candidates -> {a.out}", flush=True)


if __name__ == "__main__":
    main()
