"""Describe local curated quinone chemistry and an exact, provisional transfer.

No optimization, gene deletion, fitness values, or phenotype scores are used.
The output is evidence for a candidate arm, not an applied model correction.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
from types import SimpleNamespace

import cobra
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gembench.gene_mapping import load_genpept_table
from scripts.run_quinone_biomass_sensitivity import prepare

SOURCE = "models/bigg/iJN1463.xml"
FEATURES = "data/ncbi_feature_tables/GCF_000007565.2_Putida_feature_table.txt.gz"
CANDIDATES = ["OHPHM", "OMPHHX", "OMBZLM", "OMMBLHX", "DMQMT"]
GROUPS = {
    "aromatic_precursor": ["DDPA", "DHQS", "DHQTi_copy1", "DHQTi_copy2", "SHK3Dr", "SHKK",
                            "PSCVT_copy1", "PSCVT_copy2", "CHORS", "CHRPL"],
    "MEP_isoprenoid_precursor": ["DXPS", "DXPRIi", "MEPCT", "MEPCT_1", "CDPMEK", "MECDPS_copy1",
                                 "MECDPS_copy2", "MECDPDH5", "IPDPS", "DMPPS", "IPDDI"],
    "prenyl_chain_and_attachment": ["DMATT", "GRTT", "OCTDPS", "HBZOPT"],
    "quinone_ring_modification": ["OPHBDC", "OPHHX", *CANDIDATES],
}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def met_record(m):
    return {"id": m.id, "name": m.name, "formula": m.formula, "charge": m.charge,
            "compartment": m.compartment, "annotation": m.annotation}


def reaction_record(r, gm=None):
    record = {"id": r.id, "name": r.name, "equation": r.reaction,
              "metabolites": {m.id: coefficient for m, coefficient in sorted(r.metabolites.items(), key=lambda p: p[0].id)},
              "lower_bound": r.lower_bound, "upper_bound": r.upper_bound,
              "gene_reaction_rule": r.gene_reaction_rule,
              "gene_ids": sorted(g.id for g in r.genes), "annotation": r.annotation,
              "mass_charge_imbalance": r.check_mass_balance(),
              "missing_formula_or_charge": sorted(m.id for m in r.metabolites if not m.formula or m.charge is None)}
    if gm is not None:
        record["browser_gene_ids"] = sorted(gm.model_to_browser.get(g.id, g.id) for g in r.genes)
    return record


def equivalent_direction(a, b):
    sa = {m.id: v for m, v in a.metabolites.items()}
    sb = {m.id: v for m, v in b.metabolites.items()}
    if sa.keys() != sb.keys():
        return None
    key = next(iter(sa))
    ratio = sa[key] / sb[key]
    if not all(abs(sa[k] - ratio * sb[k]) < 1e-10 for k in sa):
        return None
    return {"stoichiometric_scale_source_over_draft": ratio,
            "source_forward_allowed_by_draft": bool(b.upper_bound > 0 if ratio > 0 else b.lower_bound < 0),
            "source_reverse_allowed_by_draft": bool(b.lower_bound < 0 if ratio > 0 else b.upper_bound > 0)}


def locus_records(loci, draft, gm, features, genpept):
    inverse = defaultdict(list)
    for model_id, locus in gm.model_to_browser.items():
        inverse[locus].append(model_id)
    records = []
    for locus in sorted(loci):
        gene_rows = features[(features["# feature"] == "gene") & features["attributes"].map(
            lambda s: locus in re.findall(r"PP_\d+", s))]
        current_tags = sorted(set(gene_rows["locus_tag"]))
        cds = features[(features["# feature"] == "CDS") & features["locus_tag"].isin(current_tags)]
        gp = genpept[genpept["locus_tags"].map(lambda s: locus in s.split(";")) |
                     genpept["old_locus_tags"].map(lambda s: locus in s.split(";"))]
        existing = sorted(inverse[locus])
        product_ids = sorted(set(x for x in cds["product_accession"] if x))
        target = existing[0] if len(existing) == 1 else (locus if not existing else None)
        records.append({"source_locus": locus, "existing_prepared_model_gene_ids": existing,
                        "target_model_gene_id": target,
                        "target_choice": "reuse_existing_model_gene" if len(existing) == 1 else
                                         ("new_locus_tag_id_for_existing_gene_mapper" if not existing else "ambiguous_do_not_apply"),
                        "current_locus_tags": current_tags, "current_refseq_protein_accessions": product_ids,
                        "current_refseq_model_style_ids": [x.replace(".", "_") for x in product_ids],
                        "refseq_ids_already_in_model": [x.replace(".", "_") for x in product_ids if x.replace(".", "_") in draft.genes],
                        "identity_mapping_ambiguous": len(existing) > 1 or len(current_tags) != 1 or len(product_ids) != 1,
                        "historical_genpept_records": gp.to_dict("records"),
                        "current_gene_records": gene_rows.to_dict("records"), "current_CDS_records": cds.to_dict("records")})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "results/quinone_repair_2026_09_06/evidence")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    destinations = [args.out / name for name in ("curated_transfer_candidate.json", "quinone_pathway_map.json",
                                                "quinone_pathway_map.tsv", "quinone_gene_aliases.tsv")]
    if any(p.exists() for p in destinations):
        raise FileExistsError("Evidence output already exists; use a fresh output directory")
    curated = cobra.io.read_sbml_model(str(ROOT / SOURCE))
    # The preparation helper needs gene names only. Do not load any fitness table.
    metadata = SimpleNamespace(genes=pd.read_table(ROOT / "data/fitness_browser/Putida/genes.tsv", dtype=str))
    draft, gm, _, _ = prepare(metadata)
    features = pd.read_table(ROOT / FEATURES, dtype=str, keep_default_na=False)
    genpept = load_genpept_table("Putida")
    ids = [rid for group in GROUPS.values() for rid in group]
    source_genes = set().union(*(set(g.id for g in curated.reactions.get_by_id(rid).genes) for rid in ids))
    genes = locus_records(source_genes, draft, gm, features, genpept)
    gene_index = {g["source_locus"]: g for g in genes}
    table = []
    for group, rids in GROUPS.items():
        for rid in rids:
            source = curated.reactions.get_by_id(rid)
            matches = []
            for target in draft.reactions:
                equivalent = equivalent_direction(source, target)
                if equivalent is not None:
                    matches.append({**reaction_record(target, gm), **equivalent})
            table.append({"group": group, "source": reaction_record(source), "same_id_present": rid in draft.reactions,
                          "exact_stoichiometric_matches": matches,
                          "structural_status": "represented_by_equivalent_reaction" if matches else "no_exact_stoichiometric_match"})
    candidate_rxns = [curated.reactions.get_by_id(rid) for rid in CANDIDATES]
    assert all(r.id not in draft.reactions for r in candidate_rxns)
    assert all(not r.boundary and not r.check_mass_balance() for r in candidate_rxns)
    missing_mets = sorted({m.id for r in candidate_rxns for m in r.metabolites if m.id not in draft.metabolites})
    candidate_genes = sorted({g.id for r in candidate_rxns for g in r.genes})
    aliases = {g: gene_index[g]["target_model_gene_id"] for g in candidate_genes}
    assert all(aliases.values()), "Ambiguous existing model identity; do not construct aliases"
    net = defaultdict(float)
    for r in candidate_rxns:
        for m, coefficient in r.metabolites.items():
            net[m.id] += coefficient
    inputs = [SOURCE, "models/gapfilled/Putida.xml.gz", "models/gapfilled/Putida_gapfill.json",
              FEATURES, "data/genpept/Putida_genpept_map.tsv", "data/fitness_browser/Putida/genes.tsv",
              "data/reference/universe_patches_v0.1.json", "data/reference/model_patches_v0.4.json",
              *[f"data/reference/gpr_patches_v0.{v}.json" for v in (2, 3, 4)],
              "scripts/map_quinone_pathway.py", "scripts/run_quinone_biomass_sensitivity.py",
              "scripts/run_carbon_fitness_generic.py", "gembench/patches.py", "gembench/gene_mapping.py"]
    provenance = {"created_utc": datetime.now(timezone.utc).isoformat(),
                  "base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                  "input_sha256": {p: sha(ROOT / p) for p in inputs}, "phenotype_values_or_scores_accessed": False,
                  "optimization_or_gene_deletion_performed": False,
                  "draft_preparation": "scripts.run_quinone_biomass_sensitivity.prepare with gene metadata only; before experimental-medium application/completion"}
    uncertainties = [
        "Exact local iJN1463 transfer is a provisional chemistry arm, not an independently validated KT2440 pathway or GPR correction.",
        "OMPHHX and OMMBLHX retain the source half-O2 equations without explicit reductant/flavin costs. Atom balance does not establish enzyme-level energetic stoichiometry. Existing OPHHX has the same simplification.",
        "Existing OPHHX is assigned to PP_5013 in both models, but the local current CDS annotation calls this protein regulatory kinase UbiB; the catalytic GPR is unresolved here.",
        "Source OPHBDC uses PP_5213 or PP_0548; current annotations identify UbiD and flavin prenyltransferase UbiX, respectively. This OR is not independent evidence of interchangeable catalytic enzymes. The draft currently retains PP_5213 alone.",
        "Source OMMBLHX uses PP_0427 or PP_5197. Their substrate/hydroxylation-site equivalence is not established by this transfer; see the separate biological-source audit before revising GPRs.",
        "q8/q8h2 is retained solely as the existing source/draft representation. Native KT2440 quinone chain length and quantitative dilution requirement remain separate questions.",
        "Source reaction upper bounds are preserved exactly (999999), not interpreted as measured capacities. No phenotype fit or numerical gap-fill selected this reaction set.",
        "Current RefSeq WP accessions are identity evidence, not automatically interchangeable model IDs. Reuse an existing draft ID when unique; otherwise the source PP locus is a new ID supported by the existing mapper."
    ]
    candidate = {"candidate_id": "ijn1463_exact_five_step_quinone_transfer_v1", "status": "provisional_not_applied",
                 "source_model": {"path": SOURCE, "sha256": sha(ROOT / SOURCE), "model_id": curated.id},
                 "reaction_ids": CANDIDATES, "gene_aliases": aliases,
                 "reaction_records": [reaction_record(r) for r in candidate_rxns],
                 "new_metabolites": {mid: met_record(curated.metabolites.get_by_id(mid)) for mid in missing_mets},
                 "all_involved_metabolites": {m.id: met_record(m) for r in candidate_rxns for m in r.metabolites},
                 "gene_identity_records": [gene_index[g] for g in candidate_genes],
                 "summed_candidate_stoichiometry": {m: v for m, v in sorted(net.items()) if v != 0},
                 "selection_basis": "Five consecutive missing reactions from the existing 2ohph_c intermediate to q8h2_c in the local curated source. All are atom/charge balanced in the source model; no artificial source or boundary reaction is included.",
                 "uncertainties": uncertainties, "provenance": provenance}
    write_json(destinations[0], candidate)

    # Record all adjacent existing alternatives at the pathway's precursor interfaces.
    interfaces = ["4hbz_c", "octdp_c", "chor_c", "h2mb4p_c", "2mecdp_c", "amet_c", "ahcys_c"]
    alternatives = {mid: [reaction_record(r, gm) for r in sorted(draft.metabolites.get_by_id(mid).reactions, key=lambda r: r.id)]
                    for mid in interfaces if mid in draft.metabolites}
    path_met_ids = sorted({m.id for rid in ids for m in curated.reactions.get_by_id(rid).metabolites})
    metabolite_comparison = []
    for mid in path_met_ids:
        source_met = curated.metabolites.get_by_id(mid)
        target_met = draft.metabolites.get_by_id(mid) if mid in draft.metabolites else None
        metabolite_comparison.append({"source": met_record(source_met), "draft": met_record(target_met) if target_met else None,
                                      "same_formula": source_met.formula == target_met.formula if target_met else None,
                                      "same_charge": source_met.charge == target_met.charge if target_met else None})
    pool = {"q8_c", "q8h2_c"}
    contributors = []
    for r in draft.reactions:
        if any(m.id in pool for m in r.metabolites):
            contributors.append({**reaction_record(r, gm), "net_pool_coefficient": sum(v for m, v in r.metabolites.items() if m.id in pool)})
    evidence = {"provenance": provenance,
                "scope": "Curated quinone route from central-carbon precursor interfaces (PEP/E4P for chorismate; G3P/pyruvate for MEP), prenyl chain synthesis, ring modification, and existing precursor/co-substrate alternatives. Upstream central metabolism and universal cofactor synthesis are outside this extraction boundary.",
                "reaction_comparison": table, "metabolite_comparison": metabolite_comparison,
                "gene_identity_records": genes, "draft_precursor_and_cosubstrate_adjacency": alternatives,
                "draft_quinone_pool_reactions": contributors,
                "all_draft_pool_reactions_conserve_combined_pool": all(r["net_pool_coefficient"] == 0 for r in contributors),
                "minimal_contiguous_source_candidate": CANDIDATES,
                "missing_candidate_metabolites": missing_mets,
                "interpretation": "No LP was used to establish minimality or growth rescue. Minimal means the contiguous absent tail in this specified source route. Exact-ID absence among upstream reactions may reflect duplicate suffixes or redox-cofactor alternatives, not missing chemistry.",
                "uncertainties": uncertainties}
    write_json(destinations[1], evidence)
    pd.DataFrame([{"group": row["group"], "source_reaction": row["source"]["id"], "source_equation": row["source"]["equation"],
                   "source_gpr": row["source"]["gene_reaction_rule"], "source_lower_bound": row["source"]["lower_bound"],
                   "source_upper_bound": row["source"]["upper_bound"], "source_imbalance": json.dumps(row["source"]["mass_charge_imbalance"]),
                   "draft_equivalents": ";".join(r["id"] for r in row["exact_stoichiometric_matches"]),
                   "draft_gprs_browser": ";".join(",".join(r["browser_gene_ids"]) for r in row["exact_stoichiometric_matches"]),
                   "status": row["structural_status"]} for row in table]).to_csv(destinations[2], sep="\t", index=False)
    pd.DataFrame([{"source_locus": g["source_locus"], "existing_model_ids": ";".join(g["existing_prepared_model_gene_ids"]),
                   "target_model_id": g["target_model_gene_id"], "target_choice": g["target_choice"],
                   "current_locus_tags": ";".join(g["current_locus_tags"]),
                   "current_refseq_proteins": ";".join(g["current_refseq_protein_accessions"]),
                   "mapping_ambiguous": g["identity_mapping_ambiguous"]} for g in genes]).to_csv(destinations[3], sep="\t", index=False)
    print(json.dumps({"candidate": str(destinations[0]), "reaction_ids": CANDIDATES, "gene_aliases": aliases,
                      "missing_metabolites": missing_mets, "source_route_reactions_mapped": len(table)}, indent=2))


if __name__ == "__main__":
    main()
