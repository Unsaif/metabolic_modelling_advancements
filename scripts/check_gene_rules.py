"""Flag suspicious gene rules in a draft model from the genome annotation alone (rules R1/R2 of
data/reference/gpr_patches_v0.2.json, made mechanical).

For every reaction whose rule has OR-alternatives, each gene's RefSeq product name is compared with
the reaction's enzyme name (from the model or the CarveMe universe). A gene is flagged when its
product name (a) contains a non-enzyme or family/domain word (hypothetical, uncharacterized, regulator,
domain, family, scaffold, elongation factor, ...), or (b) shares no informative word with the reaction
name while another alternative in the rule does. Reactions whose alternatives are annotated as 'small
subunit'/'large subunit' (or alpha/beta, component I/II) of the same enzyme are flagged as OR-joined
complexes. Output: a ranked TSV for adjudication — the tool proposes, a curator (or Claude) decides,
the benchmark verifies.

Usage: python scripts/check_gene_rules.py <org> [model.xml.gz]   (org in Btheta, Putida, MR1, Smeli)
"""
from __future__ import annotations

import gzip
import logging
import os
import re
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
logging.getLogger("cobra").setLevel(logging.ERROR)

import cobra  # noqa: E402

from gembench.gene_mapping import accession_from_model_gene, load_genpept_table  # noqa: E402

NON_ENZYME = ["hypothetical", "uncharacterized", "uncharacterised", "regulator", "domain", "family", "scaffold",
              "elongation factor", "transcription", "duf", "putative", "-like", "protein of unknown", "chaperone", "cbs "]
STOP = {"protein", "subunit", "putative", "the", "of", "and", "or", "a", "an", "type", "i", "ii", "iii", "chain", "component",
        "large", "small", "alpha", "beta", "gamma", "delta", "enzyme", "activity", "dependent", "specific", "family", "like",
        "precursor", "nadph", "nadh", "nad", "atp", "gtp", "reductase", "dehydrogenase", "synthase", "synthetase", "kinase",
        "transferase", "ligase", "hydrolase", "isomerase", "lyase", "mutase", "phosphatase", "oxidase", "aminotransferase",
        "carboxylase", "decarboxylase", "epimerase", "dehydratase", "aldolase", "cyclase", "polymerase", "reductoisomerase"}
SUBUNIT = re.compile(r"(small|large) subunit|subunit (alpha|beta|1|2|a|b)|component (i|ii|1|2)\b|(alpha|beta) (chain|subunit)", re.I)


def words(s: str):
    s = s.lower().split(" [")[0]
    toks = re.findall(r"[a-z0-9][a-z0-9'\-]*", s)
    out = set()
    for t in toks:
        t = t.strip("'-")
        if t and t not in STOP and len(t) > 2:
            out.add(t)
    return out


def main() -> None:
    org = sys.argv[1]
    default = {"Btheta": "Bacteroides_thetaiotaomicron_VPI_5482", "Putida": "Pseudomonas_putida_KT2440",
               "MR1": "Shewanella_oneidensis_MR_1", "Smeli": "Sinorhizobium_meliloti_1021"}
    path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "models", "gapfilled", f"{org}.xml.gz")
    with gzip.open(path, "rt") as fh:
        m = cobra.io.read_sbml_model(fh)
    with gzip.open(os.path.join(ROOT, "external", "carveme", "universe_bacteria.xml.gz"), "rt") as fh:
        u = cobra.io.read_sbml_model(fh)
    gp = load_genpept_table(org).set_index("version")
    lt = dict(zip(gp.index, gp["locus_tags"].str.replace("_", "", regex=False)))
    dfn = dict(zip(gp.index, gp["definition"]))
    rows = []
    for r in m.reactions:
        rule = r.gene_reaction_rule
        if " or " not in rule:
            continue
        rname = r.name if r.name and r.name != r.id else (u.reactions.get_by_id(r.id).name if r.id in u.reactions else "")
        rwords = words(rname)
        alts = [a.strip().strip("()") for a in re.split(r"\s+or\s+(?![^()]*\))", rule)]
        genes = sorted(set(re.findall(r"[A-Z]P_\d+_\d", rule)))
        info = {}
        for g in genes:
            acc = accession_from_model_gene(g)
            d = dfn.get(acc, "")
            info[g] = {"locus": lt.get(acc, g), "def": d.split(" [")[0], "words": words(d), "nonenzyme": any(k in d.lower() for k in NON_ENZYME)}
        overlaps = {g: len(info[g]["words"] & rwords) for g in genes}
        best = max(overlaps.values()) if overlaps else 0
        flagged = []
        for g in genes:
            reasons = []
            if info[g]["nonenzyme"]:
                reasons.append("non-enzyme/family/domain annotation")
            if best > 0 and overlaps[g] == 0:
                reasons.append("no word in common with the reaction name while another alternative has")
            if reasons:
                flagged.append((g, info[g]["locus"], info[g]["def"], "; ".join(reasons)))
        subunit_alts = [g for g in genes if SUBUNIT.search(info[g]["def"]) and info[g]["def"]]
        or_joined_complex = len(subunit_alts) >= 2 and any(len(a.split(" and ")) == 1 and a in subunit_alts for a in alts)
        if flagged or or_joined_complex:
            rows.append({"reaction": r.id, "reaction_name": rname, "rule": rule[:120],
                         "n_alternatives": len(alts), "flagged_genes": " | ".join(f"{loc} ({d}) — {why}" for _, loc, d, why in flagged),
                         "or_joined_subunits": " | ".join(f"{info[g]['locus']} ({info[g]['def']})" for g in subunit_alts) if or_joined_complex else "",
                         "score": len(flagged) + (2 if or_joined_complex else 0)})
    df = pd.DataFrame(rows).sort_values("score", ascending=False)
    out = os.path.join(ROOT, "results", "carbon_fitness_multi", f"gene_rule_check_{org}.tsv")
    df.to_csv(out, sep="\t", index=False)
    print(f"{org}: {len(df)} reactions flagged of {sum(1 for r in m.reactions if ' or ' in r.gene_reaction_rule)} with OR-rules; written to {os.path.relpath(out, ROOT)}")
    print(df.head(15)[["reaction", "reaction_name", "flagged_genes", "or_joined_subunits"]].to_string(index=False, max_colwidth=90))


if __name__ == "__main__":
    main()
