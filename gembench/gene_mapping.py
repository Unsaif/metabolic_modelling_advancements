"""Map model gene identifiers to Fitness Browser gene identifiers.

EMBL GEMs (CarveMe) name genes after RefSeq protein accessions ('NP_808922_1' for NP_808922.1);
the Fitness Browser names genes by locus tag ('BT0009', 'PP_0001', 'SO0020', 'SMc02791').
NCBI GenPept records carry the /locus_tag of the coding gene, so accession -> locus tag is read
from a parsed GenPept table (data/genpept/<org>_genpept_map.tsv: version, locus_tags,
old_locus_tags, gene_names, coded_by, definition) and locus tags are then normalised to the
Fitness Browser convention (underscore removed where the Browser omits it).
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GENPEPT_DIR = os.path.join(ROOT, "data", "genpept")


@dataclass
class GeneMap:
    org_id: str
    model_to_browser: Dict[str, str]           # model gene id -> Fitness Browser sysName
    unmapped_model_genes: List[str]
    stats: Dict[str, int] = field(default_factory=dict)
    provenance: Dict[str, str] = field(default_factory=dict)


def accession_from_model_gene(gene_id: str) -> str:
    """'NP_808922_1' -> 'NP_808922.1'; 'G_NP_808922_1' -> same; ids already dotted are returned as is."""
    g = gene_id[2:] if gene_id.startswith("G_") else gene_id
    if "." in g:
        return g
    head, _, ver = g.rpartition("_")
    if head and ver.isdigit():
        return f"{head}.{ver}"
    return g


def load_genpept_table(org_id: str, path: Optional[str] = None) -> pd.DataFrame:
    path = path or os.path.join(GENPEPT_DIR, f"{org_id}_genpept_map.tsv")
    df = pd.read_table(path, header=None, names=["version", "locus_tags", "old_locus_tags", "gene_names",
                                                  "coded_by", "definition"], dtype=str, keep_default_na=False)
    df["accession"] = df["version"].str.split(".").str[0]
    return df


def build_gene_map(org_id: str, model_gene_ids: List[str], browser_sysnames: Set[str],
                   genpept_path: Optional[str] = None) -> GeneMap:
    gp = load_genpept_table(org_id, genpept_path)
    by_version = dict(zip(gp["version"], gp["locus_tags"]))
    by_accession = dict(zip(gp["accession"], gp["locus_tags"]))
    old_by_version = dict(zip(gp["version"], gp["old_locus_tags"]))

    # candidate normalisations of a RefSeq locus tag to the Browser's sysName convention
    def candidates(tags: str) -> List[str]:
        out: List[str] = []
        for t in [x for x in tags.split(";") if x]:
            out += [t, t.replace("_", ""), t.replace("_", "", 1)]
        return out

    mapping: Dict[str, str] = {}
    unmapped: List[str] = []
    n_version, n_accession_only, n_no_record, n_no_sysname = 0, 0, 0, 0
    for gid in model_gene_ids:
        acc = accession_from_model_gene(gid)
        tags = by_version.get(acc)
        if tags is None:
            tags = by_accession.get(acc.split(".")[0])
            if tags is not None:
                n_accession_only += 1
        else:
            n_version += 1
        if tags is None:
            n_no_record += 1
            unmapped.append(gid)
            continue
        old = old_by_version.get(acc, "")
        hit = next((c for c in candidates(tags) + candidates(old) if c in browser_sysnames), None)
        if hit is None:
            n_no_sysname += 1
            unmapped.append(gid)
        else:
            mapping[gid] = hit
    stats = {"model_genes": len(model_gene_ids), "mapped": len(mapping), "matched_by_version": n_version,
             "matched_by_accession_only": n_accession_only, "no_genpept_record": n_no_record,
             "locus_tag_not_in_browser": n_no_sysname,
             "browser_genes_hit": len(set(mapping.values()))}
    prov = {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), "
                      "5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed",
            "genpept_table": os.path.relpath(genpept_path or os.path.join(GENPEPT_DIR, f"{org_id}_genpept_map.tsv"), ROOT)}
    return GeneMap(org_id=org_id, model_to_browser=mapping, unmapped_model_genes=unmapped, stats=stats, provenance=prov)
