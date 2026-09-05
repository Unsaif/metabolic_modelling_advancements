"""Find the reactions that let a draft model grow without an experimentally important gene.

For every gene that is important (fitness < -2) in at least half of the conditions where the wild-type
model grows but is never predicted essential, knock it out on the reference condition, run parsimonious
FBA, and record which reactions now produce the products of the gene's reactions (the 'bypass').
Aggregated over genes and organisms, the most frequent bypass reactions are the universe-level
suspects (reversibility or promiscuity errors) for the fix-at-source loop.

Output: results/carbon_fitness_multi/bypasses_<org>.tsv and bypass_summary.tsv
"""
from __future__ import annotations

import collections
import glob
import gzip
import json
import logging
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
logging.getLogger("cobra").setLevel(logging.ERROR)

import cobra  # noqa: E402
from cobra.flux_analysis import pfba  # noqa: E402

from gembench.fitness_browser import base_medium, load_organism  # noqa: E402
from gembench.gene_mapping import build_gene_map  # noqa: E402
from gembench.media import apply_medium  # noqa: E402

REFERENCE = {"Btheta": ("Varel_Bryant_medium", "EX_glc__D_e"), "Putida": ("MOPS minimal media_noCarbon", "EX_glc__D_e"),
             "MR1": ("ShewMM_noCarbon", "EX_lac__L_e"), "Smeli": ("RCH2_defined_noCarbon", "EX_glc__D_e")}
CURRENCY = {"h_c", "h2o_c", "pi_c", "ppi_c", "adp_c", "atp_c", "amp_c", "nad_c", "nadh_c", "nadp_c", "nadph_c", "co2_c",
            "nh4_c", "h_p", "h2o_p", "pi_p", "coa_c", "accoa_c", "glu__L_c", "akg_c", "o2_c", "h_e", "h2o_e"}
RES = os.path.join(ROOT, "results", "carbon_fitness_multi")


def main() -> None:
    orgs = sys.argv[1].split(",") if len(sys.argv) > 1 else list(REFERENCE)
    summary = collections.Counter(); summary_genes = collections.defaultdict(set)
    for org in orgs:
        d = glob.glob(os.path.join(RES, org, "*gapfilled__gapfilled"))[0]
        z = np.load(os.path.join(d, "matrices.npz"), allow_pickle=True)
        sim, fit, wt = z["sim_growth"], z["fitness"], z["wt_growth"]
        genes = list(z["browser_genes"]); mgenes = list(z["model_genes"])
        grows = wt >= 1e-3
        s, f = sim[:, grows], fit[:, grows]
        n_imp = np.nansum(f < -2, axis=1); pred_ess = (s < 1e-3).sum(axis=1)
        targets = [(mgenes[i], genes[i], int(n_imp[i])) for i in range(len(genes))
                   if n_imp[i] >= max(2, 0.5 * grows.sum()) and pred_ess[i] == 0]
        with gzip.open(os.path.join(ROOT, "models", "gapfilled", f"{org}.xml.gz"), "rt") as fh:
            m = cobra.io.read_sbml_model(fh)
        m.solver = "glpk"
        medname, ref_ex = REFERENCE[org]
        apply_medium(m, base_medium(medname)); m.reactions.get_by_id(ref_ex).lower_bound = -10.0
        fb = load_organism(org); desc = dict(zip(fb.genes["sysName"], fb.genes["desc"]))
        rows = []
        for mg, bg, nimp in targets:
            g = m.genes.get_by_id(mg)
            rxn_ids = [r.id for r in g.reactions]
            with m:
                g.knock_out()
                try:
                    sol = pfba(m)
                except Exception as e:  # noqa: BLE001
                    rows.append({"gene": bg, "model_gene": mg, "n_conditions_important": nimp, "desc": desc.get(bg, ""),
                                 "gene_reactions": ";".join(rxn_ids), "ko_growth": float("nan"), "bypass": f"error {e}"})
                    continue
                growth = float(sol.fluxes["Growth"]) if "Growth" in sol.fluxes.index else float("nan")
                bypass = []
                for rid in rxn_ids:
                    r = m.reactions.get_by_id(rid)
                    # the gene rule may still be satisfied by another gene (isozyme): record that
                    if r.bounds != (0.0, 0.0) and abs(sol.fluxes[rid]) > 1e-6:
                        bypass.append(f"{rid}:isozyme[{r.gene_reaction_rule[:40]}]")
                        continue
                    for met, coef in r.metabolites.items():
                        if met.id in CURRENCY:
                            continue
                        for rx in met.reactions:
                            if rx.id == rid or rx.id in rxn_ids:
                                continue
                            v = sol.fluxes[rx.id]
                            # a bypass produces what r produced (coef > 0) or consumes what r consumed (coef < 0):
                            # sign(stoichiometry * flux) must equal sign(coef)
                            if abs(v) > 1e-6 and rx.metabolites[met] * v * coef > 0:
                                role = "makes" if coef > 0 else "uses"
                                bypass.append(f"{rx.id}[{role} {met.id}; v={v:.3g}; rev={rx.reversibility}; gpr={'yes' if rx.gene_reaction_rule else 'none'}]")
                bypass = list(dict.fromkeys(bypass))
                rows.append({"gene": bg, "model_gene": mg, "n_conditions_important": nimp, "desc": desc.get(bg, ""),
                             "gene_reactions": ";".join(rxn_ids), "ko_growth": growth,
                             "bypass": " | ".join(bypass)})
                for b in bypass:
                    rid = b.split("[")[0].split(":")[0]
                    summary[(org, rid)] += 1; summary_genes[(org, rid)].add(bg)
        pd.DataFrame(rows).to_csv(os.path.join(RES, f"bypasses_{org}.tsv"), sep="\t", index=False)
        print(f"{org}: {len(targets)} target genes; top bypass reactions:",
              collections.Counter({k[1]: v for k, v in summary.items() if k[0] == org}).most_common(12), flush=True)
    agg = collections.Counter()
    for (org, rid), n in summary.items():
        agg[rid] += 1
    rows = [{"reaction": rid, "n_organisms": agg[rid],
             "n_genes_total": sum(len(summary_genes[(o, rid)]) for o in orgs if (o, rid) in summary_genes),
             "organisms": ";".join(o for o in orgs if (o, rid) in summary_genes)} for rid in agg]
    pd.DataFrame(rows).sort_values(["n_organisms", "n_genes_total"], ascending=False).to_csv(
        os.path.join(RES, "bypass_summary.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
