"""Reproduce the two-locus identity audit offline; never loads phenotypes or solves.

Run from the repository root. Prints deterministic JSON unless --out is supplied;
--out must name a fresh file. The narrowly scoped GenBank reader accepts only the
two complete, forward-strand CDS intervals downloaded for this audit.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import re


BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
SOURCES = BASE / "sequence_identity.sources"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def single(rows, predicate, description):
    found = [r for r in rows if predicate(r)]
    require(len(found) == 1, f"Expected one {description}, found {len(found)}")
    return found[0]


def read_tsv(path):
    with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open()) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def cds_record(path):
    text = path.read_text()
    starts = re.findall(r"^     CDS +(.+)$", text, re.M)
    require(len(starts) == 1, "GenBank snapshot must contain exactly one CDS")
    dna = re.sub(r"[\s0-9]", "", text.split("ORIGIN", 1)[1].split("//", 1)[0]).upper()
    require(set(dna) <= set("ACGT"), "Ambiguous nucleotide in CDS")
    require(starts[0] == f"1..{len(dna)}", "Only full forward CDS intervals are accepted")
    require('/codon_start=1' in text and '/transl_table=11' in text, "Unexpected translation settings")
    require(dna.startswith("ATG") and len(dna) % 3 == 0, "Unexpected CDS start or length")
    # NCBI table 11 has standard internal codons; this audit requires ATG start.
    amino_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    code = dict(zip(("".join(c) for c in itertools.product("TCAG", repeat=3)), amino_acids))
    translation = "".join(code[dna[i:i+3]] for i in range(0, len(dna), 3))
    require(translation.endswith("*") and "*" not in translation[:-1], "Invalid terminal/internal stop")
    aa = translation[:-1]
    deposited = re.sub(r"\s+", "", re.search(r'/translation="([^"]+)"', text).group(1))
    require(aa == deposited, "Translated DNA differs from deposited translation")
    region = re.search(r"^ACCESSION +\S+ REGION: (\d+)\.\.(\d+)$", text, re.M)
    qualifier = lambda name: re.search(rf'/{name}="([^"]+)"', text).group(1)
    return {
        "reference_version": re.search(r"^VERSION +(\S+)", text, re.M).group(1),
        "start": int(region.group(1)), "end": int(region.group(2)), "strand": "+",
        "locus_tag": qualifier("locus_tag"), "old_locus_tag": qualifier("old_locus_tag"),
        "protein_accession": qualifier("protein_id"), "product": qualifier("product"),
        "dna": dna, "translation": aa, "nucleotide_length_including_stop": len(dna),
        "protein_length": len(aa), "terminal_stop": dna[-3:],
        "dna_sha256_uppercase_no_whitespace": sha(dna.encode()),
        "protein_sha256_uppercase_no_whitespace": sha(aa.encode()),
        "translation_table": 11, "translation_matches_genbank": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.out is not None:
        require(not args.out.exists(), f"Refusing to overwrite {args.out}")
    inputs = set()
    for log in ["retrievals.json", "retrievals_followup.json"]:
        inputs.add(SOURCES / log)
        for record in json.loads((SOURCES / log).read_text()):
            if "sha256" not in record:
                continue
            path = ROOT / record["file"]
            require(sha(path.read_bytes()) == record["sha256"], f"Source changed: {path}")
            require(path.stat().st_size == record["size"], f"Source size changed: {path}")
            inputs.add(path)
    genes_path = ROOT / "data/fitness_browser/Putida/genes.tsv"
    features_path = ROOT / "data/ncbi_feature_tables/GCF_000007565.2_Putida_feature_table.txt.gz"
    genpept_path = ROOT / "data/genpept/Putida_genpept_map.tsv"
    context_path = ROOT / "results/ribosyl_disposal_2026_09_06/evidence/model_context.json"
    inputs.update([genes_path, features_path, genpept_path, context_path, Path(__file__).resolve()])
    genes, features = read_tsv(genes_path), read_tsv(features_path)
    context = json.loads(context_path.read_text())
    for path, digest in context["input_sha256"].items():
        require(sha((ROOT / path).read_bytes()) == digest, f"Model context input changed: {path}")
        inputs.add(ROOT / path)
    specifications = [
        ("PP_4248", "PP4248_NC0029474.gb", "Q88F51.json", "WP0032542781_retry.fasta"),
        ("PP_1777", "PP1777_NC0029474.gb", "PP1777_uniprot.json", "WP0109528141.fasta"),
    ]
    identities = []
    for locus, gb, uniprot_name, fasta in specifications:
        record = cds_record(SOURCES / gb)
        require(record["old_locus_tag"] == locus, "CDS locus mismatch")
        gene = single(genes, lambda r: r["locusId"] == locus, f"metadata row for {locus}")
        gene_feature = single(features, lambda r: r["# feature"] == "gene" and r["attributes"] == f"old_locus_tag={locus}", f"gene feature for {locus}")
        cds = single(features, lambda r: r["# feature"] == "CDS" and r["locus_tag"] == gene_feature["locus_tag"], f"CDS feature for {locus}")
        for row, start_key in [(gene, "begin"), (gene_feature, "start"), (cds, "start")]:
            require(int(row[start_key]) == record["start"] and int(row["end"]) == record["end"] and row["strand"] == record["strand"], f"Coordinate mismatch: {locus}")
        require(cds["genomic_accession"] == record["reference_version"] and cds["locus_tag"] == record["locus_tag"], "Reference mismatch")
        require(cds["product_accession"] == record["protein_accession"], "Protein accession mismatch")
        require(int(cds["product_length"]) == record["protein_length"], "Protein length mismatch")
        uniprot = json.loads((SOURCES / uniprot_name).read_text())
        if "results" in uniprot:
            uniprot = single(uniprot["results"], lambda r: any(x.get("value") == locus for g in r["genes"] for x in g.get("orderedLocusNames", [])), "UniProt locus result")
        protein_fasta = (SOURCES / fasta).read_text().splitlines()
        require(record["protein_accession"] in protein_fasta[0], "FASTA accession mismatch")
        fasta_aa = "".join(protein_fasta[1:])
        require(record["translation"] == uniprot["sequence"]["value"] == fasta_aa, f"Amino-acid sequence mismatch: {locus}")
        record.update({
            "local_gene_metadata": gene, "local_gene_feature": gene_feature, "local_cds_feature": cds,
            "local_scaffold_note": "genes.tsv uses unversioned AE015451; the archived feature table and newly retrieved CDS use NC_002947.4. Coordinates and old_locus_tag match. No pre-existing local genomic/CDS FASTA was available for a separate original-genome comparison.",
            "uniprot_accession": uniprot["primaryAccession"], "uniprot_entry_type": uniprot["entryType"],
            "uniprot_audit": uniprot["entryAudit"], "uniprot_protein_existence": uniprot["proteinExistence"],
            "uniprot_catalytic_annotations": [c for c in uniprot["comments"] if c["commentType"] == "CATALYTIC ACTIVITY"],
            "uniprot_references": uniprot["references"],
            "all_three_protein_sequences_exactly_identical": True,
            "comparison_sources": [gb, uniprot_name, fasta],
            "model_gene_identities": context["model_gene_identities"][locus],
            "grade": "Exact reference-sequence identity and local metadata identity; enzymatic function is evaluated separately.",
        })
        identities.append(record)
    rule = (SOURCES / "MF01537_rule.txt").read_text()
    require("Template: P0C037;" in rule and "# Version: 25" in rule, "Unexpected HAMAP snapshot")
    report = {
        "schema_version": 1, "study_role": "post-outcome source-driven candidate development",
        "scope": "Sequence/annotation metadata only. No numeric phenotype data, optimization, or model changes.",
        "all_identity_checks_passed": True, "identities": identities,
        "hamap_basis": {
            "rule": "MF_01537", "rule_version": 25, "rule_last_updated": "2024-09-03",
            "profile_version": 6, "profile_data_update": "2017-05-10", "profile_info_update": "2025-02-05",
            "template": "P0C037", "template_organism": "Escherichia coli K-12", "template_gene": "ppnP/yaiE/b0391",
            "target_annotation_evidence": "ECO:0000255 / HAMAP-Rule:MF_01537",
            "template_catalytic_evidence": "ECO:0000269 / PubMed:27941785",
            "target_matching_basis": "Q88F51 explicitly records the MF_01537 profile match and propagates its catalytic annotations; no independent profile scan or inferred match score was performed.",
            "claim_limit": "A family-rule inference from an experimentally characterized E. coli template, not direct KT2440 catalysis. Original forward/reverse assay evidence is audited separately in substrate_evidence.json.",
        },
        "ppm_evidence_limit": {
            "existing_model_reaction": next(r for r in context["metabolite_adjacency"]["r1p_c"]["reactions"] if r["id"] == "PPM"),
            "kt2440_annotation": "Q88LZ9/cpsG/PP_1777 is an unreviewed phosphomannomutase EC 5.4.2.8; its ARBA catalytic annotation supports mannose-1-P to mannose-6-P, not ribose-1-P to ribose-5-P.",
            "related_experiment": {"doi": "10.1128/jb.176.16.4851-4857.1994", "pubmed": "8050998", "organism": "Pseudomonas aeruginosa", "gene": "algC", "evidence_read": "Original primary-paper abstract reports purified AlgC converts ribose-1-P at much lower activity than its preferred substrates. This is another protein/species; neither identity nor KT2440 physiological flux is established."},
            "decision": "Existing PPM substrate/GPR assignment remains provisional. A PPM-blocking control can establish dependence of a model rescue on this edge, but cannot validate PP_1777 enzymology.",
        },
        "input_files": [{"path": str(p.relative_to(ROOT)), "size": p.stat().st_size, "sha256": sha(p.read_bytes())} for p in sorted(inputs)],
    }
    output = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.out is None:
        print(output, end="")
    else:
        with args.out.open("x") as handle:
            handle.write(output)


if __name__ == "__main__":
    main()
