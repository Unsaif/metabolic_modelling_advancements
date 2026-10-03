"""Map EMBL GEMs model genes (RefSeq protein accessions) to Fitness Browser loci by protein sequence.

EMBL drafts of newer RefSeq releases name genes by non-redundant WP_ accessions, which carry no organism-specific
locus tag, and two panel organisms are served by the Fitness Browser on a newer genome assembly than the one the
draft was built from. Mapping therefore compares sequences: model protein sequences (NCBI efetch, FASTA) against the
Fitness Browser protein sequences of the organism (orgSeqs.cgi).

Tiers: exact identity; one sequence contained in the other covering >= 90% of the longer (start-codon differences);
equal length with >= 97% identical positions. A model protein matching several loci equally is left unmapped.
Writes data/genpept_panel/<org>_genpept_map.tsv in the six-column layout read by gembench.gene_mapping
(version, locus_tags, old_locus_tags, gene_names, coded_by, definition) and a JSON summary.

Usage: python scripts/build_panel_gene_maps.py --orgs Cola,Dino,...
"""
from __future__ import annotations

import argparse
import json
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def read_fasta(path):
    recs, head, seq = [], None, []
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if head is not None:
                recs.append((head, "".join(seq).replace("*", "").upper()))
            head, seq = line[1:], []
        elif line.strip():
            seq.append(line.strip())
    if head is not None:
        recs.append((head, "".join(seq).replace("*", "").upper()))
    return recs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--orgs", required=True)
    a = ap.parse_args()
    summary = {}
    for org in a.orgs.split(","):
        model = [(h.split()[0], h.split(" ", 1)[1] if " " in h else "", s)
                 for h, s in read_fasta(os.path.join(ROOT, "data", "genpept_panel", f"{org}_model_proteins.faa"))]
        fb = []
        for h, s in read_fasta(os.path.join(ROOT, "data", "fitness_browser_panel", org, "proteins.faa")):
            parts = h.split()
            locus = parts[0].split(":", 1)[1] if ":" in parts[0] else parts[0]
            sysname = parts[1] if len(parts) > 1 and not parts[1].startswith("(") else locus
            fb.append((locus, sysname, s))
        by_seq = defaultdict(list)
        for locus, sysname, s in fb:
            by_seq[s].append(sysname or locus)
        rows, tiers = [], defaultdict(int)
        for acc, title, s in model:
            hit, tier = None, None
            if s in by_seq:
                hits = by_seq[s]
                hit, tier = (hits[0], "exact") if len(hits) == 1 else (None, "ambiguous_exact")
            else:
                cands = []
                for locus, sysname, t in fb:
                    shorter, longer = (s, t) if len(s) <= len(t) else (t, s)
                    if len(shorter) >= 0.9 * len(longer) and shorter in longer:
                        cands.append((sysname or locus, "contained"))
                    elif len(s) == len(t) and len(s) > 0:
                        same = sum(1 for x, y in zip(s, t) if x == y)
                        if same >= 0.97 * len(s):
                            cands.append((sysname or locus, "near_identical"))
                if len(cands) == 1:
                    hit, tier = cands[0]
                elif len(cands) > 1:
                    tier = "ambiguous_fuzzy"
                else:
                    tier = "unmatched"
            tiers[tier] += 1
            definition = title.replace("MULTISPECIES: ", "")
            rows.append("\t".join([acc, hit or "", "", "", "", definition]))
        out = os.path.join(ROOT, "data", "genpept_panel", f"{org}_genpept_map.tsv")
        open(out, "w").write("\n".join(rows) + "\n")
        summary[org] = {"model_proteins": len(model), "fitness_browser_proteins": len(fb), "tiers": dict(tiers)}
        print(org, summary[org], flush=True)
    json.dump(summary, open(os.path.join(ROOT, "data", "genpept_panel", "mapping_summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
