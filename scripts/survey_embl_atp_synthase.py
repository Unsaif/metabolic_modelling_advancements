"""Survey ATP synthase presence across the EMBL GEMs collection (github.com/cdanielmachado/embl_gems).

Draws a seeded random sample from the collection's model_list.tsv, downloads each model through raw.githubusercontent.com
and greps the SBML for ATP synthase reaction ids (R_ATPS*), without parsing — cheap enough for hundreds of models.

Usage: python scripts/survey_embl_atp_synthase.py [--n 500] [--seed 20260906]
Output: results/embl_atp_synthase_survey_sample<n>.tsv (assembly, taxid, organism, file, bytes, atps_ids, n_reactions)
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import random
import re
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = "https://raw.githubusercontent.com/cdanielmachado/embl_gems/master/"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--seed", type=int, default=20260906)
    args = ap.parse_args()
    rows = list(csv.DictReader(io.StringIO(urllib.request.urlopen(RAW + "model_list.tsv", timeout=120).read().decode()), delimiter="\t"))
    random.seed(args.seed)
    sample = random.sample(rows, args.n)
    out = os.path.join(ROOT, "results", f"embl_atp_synthase_survey_sample{args.n}.tsv")
    with open(out, "w") as fh:
        fh.write("assembly\ttaxid\torganism\tfile\tbytes\tatps_ids\tn_reactions\n")
        for r in sample:
            try:
                data = urllib.request.urlopen(RAW + r["file_path"], timeout=120).read()
                text = gzip.decompress(data).decode(errors="replace")
                ids = sorted(set(re.findall(r'id="(R_ATPS[A-Za-z0-9_]*)"', text)))
                n = text.count("<reaction ")
                fh.write(f"{r['assembly_accession']}\t{r['taxid']}\t{r['organism_name']}\t{os.path.basename(r['file_path'])}\t{len(data)}\t{';'.join(ids)}\t{n}\n")
            except Exception as e:  # noqa: BLE001
                fh.write(f"{r['assembly_accession']}\t{r['taxid']}\t{r['organism_name']}\t{os.path.basename(r['file_path'])}\t0\tDOWNLOAD_FAILED: {e}\t0\n")
            fh.flush()
    print("written", out)


if __name__ == "__main__":
    main()
