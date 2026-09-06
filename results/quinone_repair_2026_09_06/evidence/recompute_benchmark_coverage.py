"""Metadata-only coverage audit of the locally curated Putida quinone pathway.

Reads only metadata columns from the fitness export and only gene labels from
the saved simulation archive. Does not load numeric fitness or knockout arrays.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import cobra
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.run_quinone_biomass_sensitivity import prepare  # noqa: E402

PATHWAY = ["CHRPL", "DMATT", "GRTT", "OCTDPS", "HBZOPT", "OPHBDC",
           "OPHHX", "OHPHM", "OMPHHX", "OMBZLM", "OMMBLHX", "DMQMT"]
MISSING_TERMINAL = ["OHPHM", "OMPHHX", "OMBZLM", "OMMBLHX", "DMQMT"]
FIT_METADATA = ["orgId", "locusId", "sysName", "geneName", "desc"]
BASELINE = "results/quinone_biomass_2026_09_06/no_quinone_demand/matrices.npz"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path(__file__).with_name("benchmark_coverage.json"))
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(f"Choose a new --out path: {args.out}")
    genes = pd.read_table(ROOT / "data/fitness_browser/Putida/genes.tsv", dtype=str, keep_default_na=False)
    fit_metadata = pd.read_table(ROOT / "data/fitness_browser/Putida/fit_logratios.tsv",
                                 usecols=FIT_METADATA, dtype=str, keep_default_na=False)
    assert genes["sysName"].is_unique and fit_metadata["sysName"].is_unique
    annotated = set(genes["sysName"])
    fitness_rows = set(fit_metadata["sysName"].where(fit_metadata["sysName"] != "", fit_metadata["locusId"]))
    with np.load(ROOT / BASELINE, allow_pickle=False) as saved:
        baseline_genes = set(saved["browser_genes"])
    model, mapping, _, _ = prepare(SimpleNamespace(genes=genes))
    curated = cobra.io.read_sbml_model(str(ROOT / "models/bigg/iJN1463.xml"))
    eligible = set(mapping.model_to_browser.values()) & fitness_rows
    assert eligible == baseline_genes
    reactions = []
    for rid in PATHWAY:
        reaction = curated.reactions.get_by_id(rid)
        current = model.reactions.get_by_id(rid) if rid in model.reactions else None
        reactions.append({"reaction": rid, "name": reaction.name,
                          "curated_equation": reaction.reaction,
                          "curated_gpr": reaction.gene_reaction_rule,
                          "curated_genes": sorted(g.id for g in reaction.genes),
                          "present_in_prepared_baseline": current is not None,
                          "baseline_gpr": current.gene_reaction_rule if current is not None else None,
                          "baseline_browser_genes": sorted(mapping.model_to_browser[g.id] for g in current.genes
                                                            if g.id in mapping.model_to_browser) if current is not None else []})
    assert [r["reaction"] for r in reactions if not r["present_in_prepared_baseline"]] == MISSING_TERMINAL
    pathway_genes = sorted({gid for r in reactions for gid in r["curated_genes"]})
    terminal_genes = {gid for r in reactions if r["reaction"] in MISSING_TERMINAL for gid in r["curated_genes"]}
    records = []
    for gid in pathway_genes:
        hits = genes.loc[genes["sysName"] == gid]
        represented = sorted(mg for mg, bg in mapping.model_to_browser.items() if bg == gid)
        records.append({"gene": gid,
                        "annotated_gene_present": gid in annotated,
                        "annotation": hits.to_dict(orient="records"),
                        "curated_reactions_in_scope": [r["reaction"] for r in reactions if gid in r["curated_genes"]],
                        "on_missing_terminal_reactions": gid in terminal_genes,
                        "exported_fitness_row_present": gid in fitness_rows,
                        "baseline_model_gene_ids": represented,
                        "eligible_by_current_model_mapping_and_fitness_row": gid in eligible,
                        "in_saved_baseline_gene_axis": gid in baseline_genes,
                        "per_condition_finite_fitness_coverage": "not examined",
                        "original_mutant_library_membership": "unknown; insertion and barcode records unavailable locally",
                        "experimental_essentiality": "not inferred",
                        "knockout_prediction": "not evaluated in this metadata audit"})
    paths = ["models/bigg/iJN1463.xml", "models/gapfilled/Putida.xml.gz", "models/gapfilled/Putida_gapfill.json",
             "data/fitness_browser/Putida/genes.tsv", "data/fitness_browser/Putida/fit_logratios.tsv",
             "data/genpept/Putida_genpept_map.tsv", "data/reference/universe_patches_v0.1.json",
             "data/reference/model_patches_v0.4.json", "scripts/run_quinone_biomass_sensitivity.py",
             "scripts/run_carbon_fitness_generic.py", "gembench/gene_mapping.py", "gembench/patches.py",
             "gembench/protocols/carbon_fitness_generic.py", BASELINE]
    paths += [f"data/reference/gpr_patches_v0.{i}.json" for i in (2, 3, 4)]
    report = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "evaluation_role": "development",
        "organism": "Putida",
        "scope": "iJN1463 ubiquinone synthesis from chorismate and prenyl precursors; immediate DMATT/GRTT/OCTDPS support. Not all upstream MEP/shikimate genes or all model cofactors.",
        "selection_basis": "Local curated reaction identities and GPRs; no fitness values or simulation phenotypes used to choose this list.",
        "fitness_export_columns_loaded": FIT_METADATA,
        "saved_matrix_arrays_loaded": ["browser_genes"],
        "fitness_score_columns_loaded": False,
        "repository_commit_at_audit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "script_sha256": digest(Path(__file__)),
        "input_sha256": {path: digest(ROOT / path) for path in paths},
        "counts": {"annotated_gene_rows": len(genes), "exported_fitness_rows": len(fit_metadata),
                   "annotated_genes_without_exported_fitness_row": len(annotated - fitness_rows),
                   "baseline_model_genes": len(model.genes), "baseline_mapped_genes": len(mapping.model_to_browser),
                   "baseline_gene_axis": len(baseline_genes), "pathway_genes": len(records),
                   "pathway_genes_with_fitness_row": sum(r["exported_fitness_row_present"] for r in records),
                   "missing_terminal_reaction_genes": len(terminal_genes),
                   "missing_terminal_genes_with_fitness_row": len(terminal_genes & fitness_rows)},
        "model_mapping_statistics": mapping.stats,
        "reactions": reactions,
        "genes": records,
        "interpretation": [
            "An absent fitness row is unscored. It is neither a zero-fitness observation nor evidence of experimental nonessentiality.",
            "The downloaded gene annotation table is a genome annotation, not a mutant library inventory. Actual absence from the original mutant pool cannot be established from these inputs.",
            "Fitness row presence is metadata coverage only; it does not guarantee finite or reliable measurements in every benchmark condition.",
            "Newly adding terminal pathway GPRs cannot create experimental observations for the five distinct genes whose exported fitness rows are absent.",
            "Any later score change among covered genes measures development-set agreement; it cannot directly validate the unscored terminal pathway genes or establish generalization.",
            "A pathway patch copied from iJN1463 is supported by that reconstruction; agreement with that same reconstruction is not independent biological validation.",
            "Biological validation would require pathway/gene-specific evidence such as documented insertion coverage plus an appropriate essentiality analysis, independently measured disruption/complementation, or cofactor measurements. None was inferred here."
        ],
        "primary_method_sources": [
            {"title": "Wetmore et al. 2015, Rapid Quantification of Mutant Fitness in Diverse Bacteria by Sequencing Randomly Bar-Coded Transposons",
             "url": "https://journals.asm.org/doi/10.1128/mbio.00306-15",
             "relevance": "Library insertion coverage and genes with fitness estimates are distinct quantities (Table 1). These are general method principles, not Putida-specific insertion evidence."},
            {"title": "Price et al. 2018 author data page",
             "url": "https://genomics.lbl.gov/supplemental/bigfit/",
             "relevance": "Documents separate gene metadata, fitness-gene vectors, strain usage, and likely-essential analysis fields."}
        ]
    }
    with args.out.open("x") as fh:
        json.dump(report, fh, indent=2, allow_nan=False)
        fh.write("\n")
    print(json.dumps(report["counts"], indent=2))


if __name__ == "__main__":
    main()
