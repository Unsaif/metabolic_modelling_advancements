"""Screen ATP synthase representations in a pinned, historical EMBL GEMs snapshot.

No network is used for --summarize PATH, which also accepts the archived surveys.
New surveys record reaction-ID, name and equation candidates separately, compressed
file hashes, and sampling provenance in a companion .metadata.json file. None of
these screens establishes biological absence, usable ATP synthesis, or its cause.
Existing results are preserved unless --overwrite is explicitly supplied.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import random
import re
import sys
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gembench.survey import screen_sbml, summarize_rows  # noqa: E402

# Last upstream commit at audit, dated 2019-03-05; RefSeq release 84 collection.
DEFAULT_REF = "260d0f133802adbc151de4d7b14fb34722e1d9f4"
REPOSITORY = "https://github.com/cdanielmachado/embl_gems"
FIELDS = ["assembly", "taxid", "organism", "file", "bytes", "atps_ids", "n_reactions",
          "has_atps", "genus", "atp_synthase_name_candidates", "ion_coupled_atp_candidates",
          "status", "error", "source_path", "sha256"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--seed", type=int, default=20260906)
    ap.add_argument("--ref", default=DEFAULT_REF, help="Full immutable upstream commit SHA")
    ap.add_argument("--output", type=Path)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--summarize", type=Path, metavar="TSV", help="Audit an existing TSV offline")
    args = ap.parse_args()
    if args.summarize:
        with args.summarize.open(newline="") as stream:
            print(json.dumps(summarize_rows(list(csv.DictReader(stream, delimiter="\t"))), indent=2))
        return
    if not re.fullmatch(r"[0-9a-fA-F]{40}", args.ref):
        ap.error("--ref must be a full 40-character commit SHA, not a mutable branch or tag")
    if args.n < 1:
        ap.error("--n must be positive")
    out = args.output or ROOT / "results" / f"embl_atp_synthase_survey_sample{args.n}.tsv"
    metadata_path = out.with_suffix(".metadata.json")
    if not args.overwrite and (out.exists() or metadata_path.exists()):
        ap.error("Output already exists; use --output for a new file or explicitly --overwrite")
    raw = f"https://raw.githubusercontent.com/cdanielmachado/embl_gems/{args.ref}/"
    model_list = urllib.request.urlopen(raw + "model_list.tsv", timeout=120).read()
    population = list(csv.DictReader(io.StringIO(model_list.decode()), delimiter="\t"))
    if args.n > len(population):
        ap.error(f"--n exceeds population size {len(population)}")
    sample = random.Random(args.seed).sample(population, args.n)
    metadata = {
        "repository": REPOSITORY, "upstream_commit": args.ref,
        "model_list_sha256": hashlib.sha256(model_list).hexdigest(),
        "population_size": len(population), "seed": args.seed, "sample_size": args.n,
        "sampling": "Python random.Random(seed).sample in upstream model_list.tsv order",
        "python_version": sys.version, "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "complete": False,
        "detectors": "SBML reaction ID prefix ATPS; ATP synthase/synthetase name; BiGG-style ion-coupled ATP equation candidates. No flux or gene-annotation test.",
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    results = []
    with out.open("w" if args.overwrite else "x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        for source in sample:
            path = source["file_path"]
            row = dict(assembly=source["assembly_accession"], taxid=source["taxid"],
                       organism=source["organism_name"], file=PurePosixPath(path).name,
                       bytes=0, atps_ids="", n_reactions=0, has_atps="",
                       genus=source["organism_name"].split()[0], source_path=path,
                       sha256="", status="failed", error="")
            try:
                if not path.startswith("models/") or ".." in PurePosixPath(path).parts:
                    raise ValueError("Invalid model path in upstream manifest")
                data = urllib.request.urlopen(raw + path, timeout=120).read()
                row.update(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
                row.update(screen_sbml(gzip.decompress(data)))
                row["status"] = "ok"
            except Exception as exc:  # Failed models stay in the manifest, outside the denominator.
                row["error"] = f"{type(exc).__name__}: {exc}"
            writer.writerow(row)
            stream.flush()
            results.append(row)
    metadata.update(complete=True, summary=summarize_rows(results),
                    output_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),
                    completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata["summary"], indent=2))
    print(f"written {out} and {metadata_path}")


if __name__ == "__main__":
    main()
