"""Check 6: did the choice of models follow the plan, and was the plan committed before any result existed?

Read-only git and file-system queries: the curated model named per organism in data/fitness_browser/PROVENANCE.md
at 6cabc1e; commit times of the plan and of the scripts; file modification times of the downloads, the translation
outputs, the plan file, the scoring outputs and comparison.json; the card's own creation time; the sha256 of every
file the plan lists.

Usage: python -I check_process.py   (writes check_process.json)
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT, dump  # noqa: E402


def git(*args):
    return subprocess.run(["git", "-C", ROOT, *args], capture_output=True, text=True, check=True).stdout


def mtime(rel):
    path = os.path.join(ROOT, rel)
    if not os.path.exists(path):
        return None
    return dt.datetime.fromtimestamp(os.stat(path).st_mtime, dt.timezone.utc).isoformat(timespec="milliseconds")


def sha(rel):
    return hashlib.sha256(open(os.path.join(ROOT, rel), "rb").read()).hexdigest()


def main():
    rep = {}
    prov = git("show", "6cabc1e:data/fitness_browser/PROVENANCE.md")
    rows = {}
    for line in prov.splitlines():
        m = re.match(r"^\|\s*(\w+)\s*\|.*\|\s*(.*)\|\s*$", line)
        if m and m.group(1) in ("MR1", "Smeli", "Btheta", "Putida", "Keio", "DvH", "SynE"):
            rows[m.group(1)] = m.group(2).strip()
    rep["provenance_at_6cabc1e"] = {"commit_time": git("show", "-s", "--format=%cI", "6cabc1e").strip(), "rows": rows,
                                    "changed_since": git("log", "--format=%h", "6cabc1e..HEAD", "--", "data/fitness_browser/PROVENANCE.md").split()}
    rep["development_organisms"] = sorted(json.load(open(os.path.join(ROOT, "data/studies/transfer_v1/organisms.json")))["organisms"])
    rep["commits"] = [dict(zip(("hash", "time", "subject"), l.split("\t", 2))) for l in
                      git("log", "-6", "--format=%h\t%cI\t%s").splitlines()]
    rep["plan_commit_files"] = git("show", "--stat", "--format=", "260b92f").strip().splitlines()[-1]
    rep["last_commit_touching"] = {f: git("log", "-1", "--format=%h %cI", "--", f).strip() for f in
                                   ("docs/studies/transfer-v1-curated-references-plan.md", "scripts/translate_curated_model.py",
                                    "scripts/score_reference_model.py", "scripts/compare_reference_models.py",
                                    "gembench/protocols/carbon_fitness_generic.py", "gembench/checks.py",
                                    "data/reference/fitness_browser_media_bigg.tsv", "data/reference/fitness_browser_carbon_sources_bigg.tsv")}
    files = ["models/curated/iSO783/MODEL1507180036_url.xml", "models/curated/iGD1575/41467_2016_BFncomms12219_MOESM594_ESM.zip",
             "models/curated/iGD1575/ncomms12219-s7.xml", "scripts/translate_curated_model.py",
             "data/reference/namespace/curated_model_id_overrides.tsv", "models/curated/iGD1575/iGD1575_bigg_view.xml.gz",
             "models/curated/iSO783/iSO783_bigg_view.xml.gz", "scripts/score_reference_model.py",
             "docs/studies/transfer-v1-curated-references-plan.md", "results/transfer_v1/reference_models/Smeli/iGD1575",
             "results/transfer_v1/reference_models/MR1/iSO783/card.json", "results/transfer_v1/reference_models/MR1/iSO783/matrices.npz",
             "results/transfer_v1/reference_models/comparison.json"]
    rep["mtimes_utc"] = sorted(([mtime(f), f] for f in files), key=lambda x: x[0] or "")
    rep["card_created"] = json.load(open(os.path.join(ROOT, "results/transfer_v1/reference_models/MR1/iSO783/card.json")))["created"]
    rep["iGD1575_results_dir_contents"] = os.listdir(os.path.join(ROOT, "results/transfer_v1/reference_models/Smeli/iGD1575"))
    rep["iGD1575_logs"] = [f for f in os.listdir(os.path.join(ROOT, "logs")) if "GD1575" in f or "Smeli" in f]
    rep["uncommitted"] = git("status", "--porcelain").splitlines()
    plan = open(os.path.join(ROOT, "docs/studies/transfer-v1-curated-references-plan.md")).read()
    listed = dict(re.findall(r"`([^`]+\.(?:xml|zip))`[^|]*\|\s*(?:zip )?([0-9a-f]{64})", plan))
    rep["plan_sha256"] = {
        "iSO783 xml": [sha("models/curated/iSO783/MODEL1507180036_url.xml"), "2d4d60b9522bfe9d329ca614949e3af016ecf3a635d6df667ae8f9a1075d56ee"],
        "iGD1575 zip": [sha("models/curated/iGD1575/41467_2016_BFncomms12219_MOESM594_ESM.zip"), "1f01e9c7b0d5f03d690b53e1b3e2bb028656cd0e11daf94a6fe78381660fcb1b"],
        "iGD1575 xml": [sha("models/curated/iGD1575/ncomms12219-s7.xml"), "49a0b6e1d7a675f3fe6bc245772557d5fc4d2d58c579455832a2a9926d021894"],
        "iSO783 bigg view (card)": [sha("models/curated/iSO783/iSO783_bigg_view.xml.gz"),
                                    json.load(open(os.path.join(ROOT, "results/transfer_v1/reference_models/MR1/iSO783/card.json")))["model"]["sha256"]]}
    rep["plan_sha256_all_match"] = all(a == b for a, b in rep["plan_sha256"].values())
    rep["plan_text_regex_hits"] = listed
    for k, v in rep.items():
        print(k, json.dumps(v)[:400])
    dump(rep, "check_process.json")


if __name__ == "__main__":
    main()
