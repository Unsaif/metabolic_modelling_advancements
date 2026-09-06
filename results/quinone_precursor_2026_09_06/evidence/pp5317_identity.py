"""Offline PP_5317 identity audit. Reads sequence/annotation metadata only.

Run from any directory; prints deterministic JSON. --out refuses existing files.
This deliberately narrow reader selects the full PP_5317 CDS, excluding the
downstream partial ubiA feature returned by NCBI for this interval.
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
SOURCES = BASE / "pp5317_identity.sources"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def single(rows, predicate, label):
    matches = [row for row in rows if predicate(row)]
    require(len(matches) == 1, f"Expected one {label}, found {len(matches)}")
    return matches[0]


def tsv(path):
    with gzip.open(path, "rt") if path.suffix == ".gz" else path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def extract():
    inputs = set(SOURCES.iterdir()) | {Path(__file__).resolve()}
    for log in SOURCES.glob("retrievals*.json"):
        for row in json.loads(log.read_text()):
            path = Path(row.get("path", "/not-a-snapshot"))
            if path.is_absolute() or "sha256" not in row:
                continue  # Temporary full articles are intentionally not distributed.
            path = ROOT / path
            require(sha(path.read_bytes()) == row["sha256"], f"Changed snapshot: {path}")
            require(path.stat().st_size == row["size"], f"Changed size: {path}")

    text = (SOURCES / "PP5317_NC0029474.gb").read_text()
    dna = re.sub(r"[\s0-9]", "", text.split("ORIGIN", 1)[1].split("//", 1)[0]).upper()
    require(len(dna) == 558 and set(dna) <= set("ACGT"), "Unexpected DNA length/alphabet")
    blocks = re.findall(r"^     CDS +([^\n]+)\n(.*?)(?=^     \S|^ORIGIN)", text, re.M | re.S)
    _, block = single(blocks, lambda x: x[0] == "1..558" and '/old_locus_tag="PP_5317"' in x[1], "complete PP_5317 CDS")
    require('/codon_start=1' in block and '/transl_table=11' in block, "Unexpected translation settings")
    require(dna[:3] == "GTG" and dna[-3:] == "TGA", "Unexpected terminal codons")
    code = dict(zip(("".join(c) for c in itertools.product("TCAG", repeat=3)),
                    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
    # NCBI genetic code 11: GTG initiates methionine; internal GTG remains valine.
    translated = "M" + "".join(code[dna[i:i+3]] for i in range(3, len(dna), 3))
    require(translated.endswith("*") and "*" not in translated[:-1], "Invalid stop")
    protein = translated[:-1]
    qualifier = lambda key: re.search(rf'/{key}="([^\"]+)"', block).group(1)
    deposited = re.sub(r"\s+", "", qualifier("translation"))
    uniprot = single(json.loads((SOURCES / "PP5317_uniprot.json").read_text())["results"],
                     lambda x: x["primaryAccession"] == "Q88C66", "Q88C66")
    sequences = {"translated_NC_002947.4_CDS": protein, "NCBI_CDS_translation": deposited,
                 "UniProt_Q88C66_sequence_version_1": uniprot["sequence"]["value"]}
    for name, accession in [("NP7474181.fasta", "NP_747418.1"), ("WP0109558101.fasta", "WP_010955810.1")]:
        lines = (SOURCES / name).read_text().splitlines()
        require(accession in lines[0], "Wrong FASTA accession")
        sequences[accession] = "".join(lines[1:])
    require(len(set(sequences.values())) == 1 and len(protein) == 185, "Protein identity mismatch")

    paths = {"genes": ROOT / "data/fitness_browser/Putida/genes.tsv",
             "features": ROOT / "data/ncbi_feature_tables/GCF_000007565.2_Putida_feature_table.txt.gz",
             "genpept": ROOT / "data/genpept/Putida_genpept_map.tsv",
             "mapping": ROOT / "results/quinone_repair_2026_09_06/evidence/quinone_pathway_map.json"}
    inputs.update(paths.values())
    genes, features = tsv(paths["genes"]), tsv(paths["features"])
    local_gene = single(genes, lambda r: r["locusId"] == "PP_5317", "gene metadata")
    current_gene = single(features, lambda r: r["# feature"] == "gene" and r["attributes"] == "old_locus_tag=PP_5317", "gene feature")
    current_cds = single(features, lambda r: r["# feature"] == "CDS" and r["locus_tag"] == current_gene["locus_tag"], "CDS feature")
    for row, key in [(local_gene, "begin"), (current_gene, "start"), (current_cds, "start")]:
        require((int(row[key]), int(row["end"]), row["strand"]) == (6063112, 6063669, "+"), "Coordinate mismatch")
    require(current_cds["product_accession"] == qualifier("protein_id") == "WP_010955810.1", "Accession mismatch")
    historical = single(list(csv.reader(paths["genpept"].open(), delimiter="\t")), lambda r: r[0] == "NP_747418.1", "historical GenPept map")
    require(historical[1] == "PP_5317" and historical[4] == "NC_002947.4:6063112..6063669", "Historical mapping mismatch")
    mapping = single(json.loads(paths["mapping"].read_text())["gene_identity_records"],
                     lambda r: r["source_locus"] == "PP_5317", "prepared model mapping")
    require(mapping["existing_prepared_model_gene_ids"] == ["NP_747418_1"] and not mapping["identity_mapping_ambiguous"], "Ambiguous model identity")

    forward, reverse = "TCGTACGAATCCCCG", "TCAGCGGTTTTCCTCCTTG"
    reverse_target = reverse.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    require(dna[3:3+len(forward)] == forward and dna.endswith(reverse_target), "Kitade primer end mismatch")
    rule = (SOURCES / "MF01632_rule.txt").read_text()
    require("Template: P26602;" in rule and "# Version: 27" in rule, "Unexpected HAMAP rule")
    alternative = single(json.loads((SOURCES / "KT2440_ubiC_chorismatase_search.json").read_text())["results"],
                         lambda r: r["primaryAccession"] == "Q88GE0", "alternative domain record")
    alt_gene = single(genes, lambda r: r["locusId"] == "PP_3784", "PP_3784 metadata")
    alt_feature = single(features, lambda r: r["# feature"] == "gene" and r["attributes"] == "old_locus_tag=PP_3784", "PP_3784 feature")

    return {
        "audit_date": "2026-09-06", "scope": "Sequence and source evidence only; no phenotype values or optimization.",
        "identity": {"locus": "PP_5317", "gene": "ubiC", "current_locus": qualifier("locus_tag"),
            "reference": "NC_002947.4", "coordinates": [6063112, 6063669], "strand": "+",
            "protein_accessions": list(sequences), "all_five_sequences_exactly_identical": True,
            "protein_length": len(protein), "protein_sequence": protein, "protein_sha256": sha(protein.encode()),
            "dna_length_including_stop": len(dna), "dna_sequence": dna, "dna_sha256": sha(dna.encode()),
            "start_codon": "GTG", "translation_table": 11, "initiator_amino_acid": "M",
            "local_gene_metadata": local_gene, "current_gene_feature": current_gene, "current_cds_feature": current_cds,
            "historical_genpept": historical, "model_mapping": mapping,
            "limitation": "Local annotation coordinates/accessions match newly retrieved reference sequences. No original local genomic/CDS FASTA was available for an independent original-genome comparison."},
        "uniprot": {"accession": "Q88C66", "entry_type": uniprot["entryType"], "audit": uniprot["entryAudit"],
            "protein_existence": uniprot["proteinExistence"], "description": uniprot["proteinDescription"],
            "catalytic_activity": [r for r in uniprot["comments"] if r["commentType"] == "CATALYTIC ACTIVITY"],
            "references": uniprot["references"], "grade": "Homology annotation, not a KT2440 catalytic experiment."},
        "hamap": {"rule": "MF_01632", "version": 27, "template": "P26602 / E. coli UbiC",
            "scope_note": "The rule assigns Probable outside Enterobacterales; template has enzyme experiments. Reviewed does not mean the KT2440 protein was assayed."},
        "kitade_2018": {"doi": "10.1128/AEM.02587-17", "pmid": "29305513", "pmcid": "PMC5835730",
            "evidence": "P. putida genomic ubiC PCR product expressed in wild-type C. glutamicum; crude-extract activity measured. Donor strain/accession not specified for the screen. S12 is named only for tolerance assays.",
            "table_1_activity": {"mean": 142, "sd": 12, "n": 5, "units": "nmol mg^-1 min^-1"},
            "table_5_primer_67": "CTCTCATATG" + forward, "table_5_primer_68": "CTCTCATATG" + reverse,
            "primer_checks": {"forward_genomic_positions_1_based": [4, 3+len(forward)],
                "reverse_genomic_positions_1_based": [len(dna)-len(reverse)+1, len(dna)],
                "both_binding_ends_exactly_match": True, "full_construct_sequence_identity_proven": False,
                "note": "NdeI introduces ATG in place of genomic GTG; expected initiator methionine is unchanged. End matching cannot establish donor strain or internal sequence."},
            "grade": "Direct heterologous activity for a P. putida donor; PP_5317-compatible primers, unresolved exact construct."},
        "alternative_candidate": {"locus": "PP_3784", "uniprot": "Q88GE0", "current_locus": alt_feature["locus_tag"],
            "local_metadata": alt_gene, "current_gene_feature": alt_feature,
            "description": alternative["proteinDescription"], "audit": alternative["entryAudit"],
            "protein_existence": alternative["proteinExistence"], "protein_length": alternative["sequence"]["length"],
            "sequence_sha256": sha(alternative["sequence"]["value"].encode()),
            "catalytic_annotations": [r for r in alternative.get("comments", []) if r["commentType"] == "CATALYTIC ACTIVITY"],
            "grade": "Distinct uncharacterized KT2440 locus with Pfam FkbO/Hyg5-like domain assignment. Neither UbiC alias nor demonstrated 4-hydroxybenzoate synthase; no replacement GPR proposed."},
        "source_conclusion": "PP_5317 identity is resolved; UbiC assignment is well-supported homology with compatible species-level enzyme evidence, but exact KT2440 catalysis remains unverified. No native alternative is established by this bounded search.",
        "other_primary_evidence": [
            {"doi": "10.3389/fbioe.2016.00090", "pmid": "27965953", "authors_year": "Yu et al. 2016",
             "organism_distinction": "KT2440 production host; heterologous E. coli K-12 W3110 UbiC, accession CAA40681.1. Does not assay native PP_5317."},
            {"doi": "10.1021/acssynbio.8b00465", "pmid": "30861344", "authors_year": "Jha et al. 2019",
             "organism_distinction": "KT2440 sensor/production host; engineered E. coli UbiC. Discussion identifies PP_5317 as an annotated homolog and speculates about expression; no native PP_5317 catalytic demonstration."},
            {"doi": "10.1111/mmi.12084", "pmid": "23113660", "authors_year": "Zhou et al. 2013",
             "organism_distinction": "Xanthomonas campestris pv. campestris XanB2/Xcc4014 produces 3- and 4-hydroxybenzoate and supplies ubiquinone precursor. Direct evidence in another organism, not a KT2440 locus assignment."},
            {"doi": "10.1021/jacs.5b05559", "pmid": "26247872", "authors_year": "Hubrich et al. 2015",
             "inference_limit": "Related chorismatase folds can yield different products; experimentally tested residues affect specificity. FkbO/Hyg5-like domain membership alone cannot assign the 4-hydroxybenzoate reaction."},
        ],
        "search_limit": "An annotation query is not an exhaustive sequence search and cannot exclude unannotated paralogs or bypasses.",
        "input_sha256": {str(p.relative_to(ROOT)): sha(p.read_bytes()) for p in sorted(inputs)},
        "article_snapshot_policy": "Official database records are saved. Full copyrighted articles were read temporarily, with retrieval hashes retained; they are not needed for offline sequence recomputation and are not distributed.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    payload = json.dumps(extract(), indent=2, sort_keys=True) + "\n"
    if args.out:
        with args.out.open("x") as handle:
            handle.write(payload)
    else:
        print(payload, end="")
