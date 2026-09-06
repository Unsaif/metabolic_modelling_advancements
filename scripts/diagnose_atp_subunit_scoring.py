"""Run a controlled missing-subunit experiment on two pinned CarveMe sources.

Requires an existing public CarveMe Git clone containing the two source blobs.
Reads source through Git and executes that inspected, pinned scoring module; it
never reads phenotype data, aligns genomes, reconstructs models, or downloads.
The output is evidence of algorithm behavior, not of the cause in any organism.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import warnings

import pandas as pd

BEFORE = "86e1ef8a39f3f1af32309162894e4f8d999bf955"
AFTER = "5038ece6e445862eb5e73ddf1bb9258e67a780ff"
SOURCE_PATH = "carveme/reconstruction/scoring.py"
SCENARIOS = ("complete", "one_missing", "one_weak", "all_missing", "alternative_complete")


def fixture(scenario: str):
    """Two eight-subunit reference complexes plus an unrelated scoring control."""
    if scenario not in SCENARIOS:
        raise ValueError(f"Unknown scenario: {scenario}")
    gprs, annotation = [], []
    for reference in ("refA", "refB"):
        for index in range(8):
            gene = f"subunit{index}"
            gprs.append(dict(gene="G_" + gene, protein="P_complex", reaction="R_ATPS4rpp", model=reference))
            missing = scenario == "all_missing" or (index == 7 and scenario == "one_missing")
            missing |= scenario == "alternative_complete" and reference == "refA" and index == 7
            if not missing:
                score = 1.0 if scenario == "one_weak" and index == 7 else 80.0
                annotation.append(dict(BiGG_gene=f"{reference}.{gene}", query_gene=f"target_{gene}", score=score))
    gprs.append(dict(gene="G_control", protein="P_control", reaction="R_control", model="refA"))
    annotation.append(dict(BiGG_gene="refA.control", query_gene="target_control", score=80.0))
    return pd.DataFrame(annotation), pd.DataFrame(gprs)


def read_blob(repository: Path, commit: str, path: str) -> bytes:
    env = dict(os.environ, GIT_NO_LAZY_FETCH="1")
    return subprocess.run(["git", "-C", str(repository), "show", f"{commit}:{path}"],
                          check=True, capture_output=True, env=env).stdout


def run_source(source: bytes) -> list[dict]:
    namespace = {"__name__": "pinned_carveme_scoring"}
    exec(compile(source, "<pinned CarveMe scoring.py>", "exec"), namespace)
    output = []
    for scenario in SCENARIOS:
        annotation, gprs = fixture(scenario)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            scored = namespace["reaction_scoring"](annotation, gprs)
        match = scored.loc[scored.reaction == "R_ATPS4rpp"]
        output.append(dict(scenario=scenario, reaction_scored=not match.empty,
                           raw_score=None if match.empty else float(match.iloc[0].score),
                           gpr=None if match.empty else str(match.iloc[0].GPR)))
    return output


def validate(result: dict) -> None:
    """Check pre-specified consequences against unmodified upstream functions."""
    expected = {"before": [80.0, None, 1.0, None, 80.0],
                "after": [80.0, 70.0, 70.125, None, 80.0]}
    for version, scores in expected.items():
        if tuple(row["scenario"] for row in result[version]["scenarios"]) != SCENARIOS:
            raise AssertionError(f"Incomplete or reordered {version} experiment")
        for row, score in zip(result[version]["scenarios"], scores):
            actual = row["raw_score"]
            if (score is None and actual is not None) or (score is not None and (actual is None or not math.isclose(actual, score))):
                raise AssertionError(f"Unexpected {version}/{row['scenario']} score: {actual}, expected {score}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--carveme-repo", type=Path, required=True)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    result = {"experiment": "Synthetic alignment/GPR perturbation; not an organism reconstruction",
              "pandas_version": pd.__version__, "scoring_checks_passed": False}
    for label, commit in (("before", BEFORE), ("after", AFTER)):
        source = read_blob(args.carveme_repo, commit, SOURCE_PATH)
        result[label] = dict(commit=commit, source_sha256=hashlib.sha256(source).hexdigest(),
                             url=f"https://github.com/cdanielmachado/carveme/blob/{commit}/{SOURCE_PATH}",
                             scenarios=run_source(source))
    gpr_blob = read_blob(args.carveme_repo, BEFORE, "carveme/data/generated/bigg_gprs.csv.gz")
    metadata_blob = read_blob(args.carveme_repo, BEFORE, "carveme/data/input/bigg_models.csv")
    templates = [r for r in csv.DictReader(io.StringIO(gzip.decompress(gpr_blob).decode())) if r["reaction"] == "R_ATPS4rpp"]
    models = {r["bigg_id"]: r for r in csv.DictReader(io.StringIO(metadata_blob.decode()))}
    model_ids = {r["model"] for r in templates}
    result["historical_atps_templates"] = dict(
        commit=BEFORE, gpr_file_sha256=hashlib.sha256(gpr_blob).hexdigest(),
        model_metadata_sha256=hashlib.sha256(metadata_blob).hexdigest(),
        gene_rows=len(templates), models=len(model_ids),
        complexes=len({(r["model"], r["protein"]) for r in templates}),
        models_per_genus=dict(sorted(Counter(models[m]["organism"].split()[0] for m in model_ids).items())),
        templates_include_development_organism="Pseudomonas putida KT2440" in {models[m]["organism"] for m in model_ids})
    validate(result)
    result["scoring_checks_passed"] = True
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as stream:
            stream.write(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
