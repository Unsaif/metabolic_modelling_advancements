"""Check 3: gene mapping of the curated models to Fitness Browser sysNames.

iSO783 maps by identity; iGD1575 by smc/sma/smb -> SMc/SMa/SM_b. This script re-implements both rules, reports the
rates (all genes, genes with fitness data), lists the unmapped genes with any near match in the Fitness Browser
table (underscore variants, case, locusId), checks that no two model genes map to one Fitness Browser gene, and
spot-checks ten mapped genes per model with fitness data: the Fitness Browser description next to the reactions the
model gives that gene.

Usage: python -I check_gene_mapping.py   (writes check_gene_mapping.json and gene_spot_checks.tsv)
"""
from __future__ import annotations

import collections
import os
import random
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import MODELS, OUT, dump, p, read_model  # noqa: E402


def sm_prefix(g):
    if g[:3].lower() == "smb":
        return "SM_b" + g[3:]
    if g[:2].lower() == "sm":
        return "SM" + g[2:]
    return g


RULES = {"iSO783": lambda g: g, "iGD1575": sm_prefix}


def fitness_genes(org):
    fit = pd.read_table(p(f"data/fitness_browser/{org}/fit_logratios.tsv"), dtype=str, keep_default_na=False,
                        usecols=["locusId", "sysName"])
    return set(fit["sysName"].where(fit["sysName"] != "", fit["locusId"]))


def universe_reaction_names():
    import gzip
    text = gzip.open(p("results/ppnp_repair_2026_09_06/evidence/source/universe_bacteria.xml.gz"), "rt").read()
    return {m.group(1): m.group(2) for m in re.finditer(r'<reaction[^>]*\bid="R_([^"]+)"[^>]*\bname="([^"]*)"', text)}


def main():
    report = {}
    rows = []
    uni_rxn = universe_reaction_names()
    for label, cfg in MODELS.items():
        org = cfg["org"]
        genes = pd.read_table(p(f"data/fitness_browser/{org}/genes.tsv"), dtype=str, keep_default_na=False)
        sysnames = set(genes["sysName"])
        desc = dict(zip(genes["sysName"], genes["desc"]))
        fit = fitness_genes(org)
        m = read_model(cfg["view"])
        rule = RULES[label]
        ids = [g.id for g in m.genes]
        mapped = {g: rule(g) for g in ids if rule(g) in sysnames}
        unm = [g for g in ids if g not in mapped]
        dup = {k: v for k, v in collections.Counter(mapped.values()).items() if v > 1}
        with_fit = {g: s for g, s in mapped.items() if s in fit}
        near = []
        low = {s.lower().replace("_", ""): s for s in sysnames}
        loc = dict(zip(genes["locusId"], genes["sysName"]))
        for g in unm:
            cands = set()
            key = rule(g).lower().replace("_", "")
            if key in low:
                cands.add(low[key])
            if g in loc:
                cands.add(loc[g])
            near.append({"gene": g, "rule_result": rule(g), "near_matches": sorted(cands),
                         "reactions": [r.id for r in m.genes.get_by_id(g).reactions]})
        rep = {"model_genes": len(ids), "mapped": len(mapped), "mapped_rate": len(mapped) / len(ids),
               "mapped_with_fitness_data": len(with_fit), "mapped_with_fitness_rate": len(with_fit) / len(ids),
               "fitness_browser_genes_with_fitness": len(fit),
               "duplicate_targets": dup, "unmapped": near,
               "unmapped_with_a_near_match": [x for x in near if x["near_matches"]]}
        # spot checks
        rng = random.Random(20261006)
        pick = sorted(with_fit)
        rng.shuffle(pick)
        spots = []
        for g in pick[:10]:
            rx = m.genes.get_by_id(g).reactions
            names = []
            for r in sorted(rx, key=lambda r: r.id):
                nm = r.name if r.name and r.name != r.id else uni_rxn.get(r.id, "")
                names.append(f"{r.id} ({nm})" if nm else r.id)
            spots.append({"model_gene": g, "sysName": with_fit[g], "fitness_browser_desc": desc.get(with_fit[g], ""),
                          "model_reactions": names})
            rows.append({"model": label, "model_gene": g, "sysName": with_fit[g], "fitness_browser_desc": desc.get(with_fit[g], ""),
                         "model_reactions": " | ".join(names)})
        rep["spot_checks"] = spots
        report[label] = rep
        print(label, {k: rep[k] for k in ("model_genes", "mapped", "mapped_with_fitness_data", "duplicate_targets")},
              "unmapped:", [x["gene"] for x in near], "near:", rep["unmapped_with_a_near_match"], flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "gene_spot_checks.tsv"), sep="\t", index=False)
    dump(report, "check_gene_mapping.json")


if __name__ == "__main__":
    main()
