"""Freeze and complete the eight previously specified cycle 6/7 development runs.

Run configurations and input hashes are written before simulation. Runs execute
sequentially with two deletion workers. Only validated complete attempts become
reusable results; partial/failed attempts are preserved rather than overwritten.
These data were inspected during curation, so this is retrospective development.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gembench.comparison import aligned_pairs, compare_runs, load_run
from gembench.protocols.carbon_fitness_generic import GenericParams

ORGS = ("Btheta", "Putida", "MR1", "Smeli")
GPR = ",".join(f"data/reference/gpr_patches_v0.{v}.json" for v in (2, 3, 4))
CORE = ["gembench/__init__.py", "gembench/cards.py", "gembench/checks.py", "gembench/comparison.py",
        "gembench/fitness_browser.py", "gembench/gene_mapping.py", "gembench/media.py", "gembench/metrics.py",
        "gembench/patches.py", "gembench/protocols/__init__.py", "gembench/protocols/carbon_fitness_generic.py",
        "gembench/protocols/ecoli_carbon_fitness.py", "scripts/run_carbon_fitness_generic.py",
        "scripts/run_development_sprint.py"]
REQUIRED = ("card.json", "card.md", "matrices.npz", "conditions.tsv", "per_condition_metrics.tsv",
            "gene_map.tsv", "gene_map_unmapped.txt", "dropped_rich_medium_essentials.txt")


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe(value):
    if isinstance(value, dict):
        return {str(k): safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe(v) for v in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(safe(value), indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def environment():
    return {"python": sys.version, "packages": {name: version(name) for name in
            ("cobra", "numpy", "scipy", "pandas", "scikit-learn", "optlang", "swiglpk", "python-libsbml")}}


def input_paths():
    paths = {ROOT / p for p in CORE}
    paths.update((ROOT / "data/reference").glob("*.tsv"))
    paths.update(ROOT / p for p in GPR.split(","))
    paths.update(ROOT / f"data/reference/model_patches_v0.{v}.json" for v in (3, 4))
    paths.add(ROOT / "data/reference/universe_patches_v0.1.json")
    for org in ORGS:
        paths.update(p for p in (ROOT / f"data/fitness_browser/{org}").rglob("*") if p.is_file())
        paths.update(ROOT / p for p in (f"models/gapfilled/{org}.xml.gz", f"models/gapfilled/{org}_gapfill.json",
                                       f"data/genpept/{org}_genpept_map.tsv"))
    return sorted(paths)


def prepare(out):
    path = out / "prespecified_manifest.json"
    if path.exists():
        manifest = json.loads(path.read_text())
        verify_inputs(manifest)
        return manifest
    params = asdict(GenericParams(processes=2, solver="glpk", drop_rich_medium_essentials=False,
                                  complete_medium_transport=True))
    configurations = []
    for cycle, patch in ((6, "0.3"), (7, "0.4")):
        for org in ORGS:
            command = [sys.executable, "scripts/run_carbon_fitness_generic.py", "--processes", "2", "--solver", "glpk",
                       "--orgs", org, "--variant", "gapfilled", "--no-drop-rich", "--complete-medium-transport",
                       "--patch", "data/reference/universe_patches_v0.1.json", "--model-patch",
                       f"data/reference/model_patches_v{patch}.json", "--gpr-patch", GPR,
                       "--output-dir", "{attempt_directory}/artifacts"]
            configurations.append({"run_id": f"{org}_cycle{cycle}", "organism": org, "cycle": cycle,
                                   "model_patch_version": patch, "params": params, "command": command})
    manifest = {"schema_version": 1, "frozen_at_utc": now(), "base_commit": subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "evaluation_role": "retrospective_development", "environment": environment(),
                "configurations": configurations,
                "input_sha256": {str(p.relative_to(ROOT)): sha(p) for p in input_paths()},
                "analysis": {"comparison": "cycle7 minus cycle6 within each organism; shared genes, growing conditions and finite observations",
                             "primary_metric": "MCC", "bootstrap": {"unit": "gene", "n": 500, "seed": 0},
                             "coverage": "report every mapped condition, including wild-type failures",
                             "material_regression_screen": "any lost growing condition or negative matched MCC difference; descriptive attribution only",
                             "selection_policy": "Run all eight configurations. Do not amend, select, accept or reject patches from these outcomes.",
                             "independence": "These organisms/phenotypes were used during patch development. Freezing this completion does not make them held-out data."}}
    write_json(path, manifest)
    print(f"FROZEN: {len(configurations)} runs, {len(manifest['input_sha256'])} input hashes at {path}", flush=True)
    return manifest


def verify_inputs(manifest):
    changed = [p for p, digest in manifest["input_sha256"].items() if not (ROOT / p).is_file() or sha(ROOT / p) != digest]
    if changed:
        raise RuntimeError(f"Frozen inputs changed; refusing to run/resume: {changed}")
    if environment() != manifest["environment"]:
        raise RuntimeError("Frozen Python/dependency environment changed; refusing to run/resume")


def validate(directory, config, manifest):
    marker = directory / "completion.json"
    if marker.exists():
        saved = json.loads(marker.read_text())
        if saved["run_id"] != config["run_id"]:
            raise RuntimeError("Completed run identifier differs from configuration")
        for name, digest in saved["artifact_sha256"].items():
            if sha(directory / name) != digest:
                raise RuntimeError(f"Completed run artifact changed: {name}")
    cards = list((directory / "artifacts" / config["organism"]).glob("*/card.json"))
    if len(cards) != 1:
        raise RuntimeError(f"Expected exactly one completed card for {config['run_id']}")
    card_path = cards[0]
    for name in REQUIRED:
        if not card_path.with_name(name).is_file():
            raise RuntimeError(f"Missing run artifact: {name}")
    card = json.loads(card_path.read_text())
    if card["protocol"]["params"] != config["params"]:
        raise RuntimeError("Completed run parameters differ from prespecified configuration")
    if card["protocol"]["patches"]["model_file"] != f"data/reference/model_patches_v{config['model_patch_version']}.json":
        raise RuntimeError("Completed run used a different model patch")
    fingerprints = card["protocol"]["input_sha256"]
    for path, digest in manifest["input_sha256"].items():
        if path in fingerprints and fingerprints[path] != digest:
            raise RuntimeError(f"Completed run input fingerprint differs: {path}")
    run = load_run(card_path.parent)
    if not np.isfinite(run["sim_growth"]).all() or not np.isfinite(run["wt_growth"]).all():
        raise RuntimeError("Nonfinite simulation output")
    if not np.isfinite(run["fitness"]).any():
        raise RuntimeError("No finite experimental observations")
    if not np.isclose(card["results"]["condition_level"]["n_conditions_wt_grows"],
                      np.sum(run["wt_growth"] >= config["params"]["growth_threshold"])):
        raise RuntimeError("Card/array coverage mismatch")
    return card_path.parent


def run_all(out, manifest):
    for config in manifest["configurations"]:
        verify_inputs(manifest)
        final = out / "runs" / config["run_id"]
        marker = final / "completion.json"
        if final.exists():
            if not marker.exists():
                raise RuntimeError(f"Existing result lacks completion marker: {final}")
            saved = json.loads(marker.read_text())
            if saved["manifest_sha256"] != sha(out / "prespecified_manifest.json"):
                raise RuntimeError("Existing result belongs to another manifest")
            validate(final, config, manifest)
            for name, digest in saved["artifact_sha256"].items():
                if sha(final / name) != digest:
                    raise RuntimeError(f"Existing result artifact changed: {name}")
            print(f"RESUME: validated {config['run_id']}", flush=True)
            continue
        attempt = out / "attempts" / f"{config['run_id']}_{time.time_ns()}"
        attempt.mkdir(parents=True)
        command = [part.replace("{attempt_directory}", str(attempt)) for part in config["command"]]
        start = now(); t0 = time.monotonic()
        write_json(attempt / "invocation.json", {"command": command, "cwd": str(ROOT), "started_utc": start,
                                                "manifest_sha256": sha(out / "prespecified_manifest.json")})
        print(f"START: {config['run_id']} {start}", flush=True)
        with (attempt / "execution.log").open("x") as log:
            process = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        record = {"run_id": config["run_id"], "started_utc": start, "finished_utc": now(),
                  "seconds": time.monotonic() - t0, "exit_code": process.returncode,
                  "manifest_sha256": sha(out / "prespecified_manifest.json")}
        try:
            if process.returncode:
                raise RuntimeError(f"Runner exited with {process.returncode}; see {attempt / 'execution.log'}")
            verify_inputs(manifest)
            validate(attempt, config, manifest)
            record["artifact_sha256"] = {str(p.relative_to(attempt)): sha(p) for p in sorted((attempt / "artifacts").rglob("*")) if p.is_file()}
            write_json(attempt / "completion.json", record)
            final.parent.mkdir(parents=True, exist_ok=True)
            attempt.rename(final)
            print(f"COMPLETE: {config['run_id']} ({record['seconds']:.1f}s)", flush=True)
        except Exception as error:
            record["failure"] = str(error)
            write_json(attempt / "failure.json", record)
            print(f"FAILED: {config['run_id']}: {error}", flush=True)


def summarize(out, manifest):
    verify_inputs(manifest)
    rows, by_id = [], {}
    for config in manifest["configurations"]:
        directory = out / "runs" / config["run_id"]
        if not (directory / "completion.json").exists():
            rows.append({"run_id": config["run_id"], "status": "incomplete_or_failed"})
            continue
        run_dir = validate(directory, config, manifest)
        card = json.loads((run_dir / "card.json").read_text())
        completion = json.loads((directory / "completion.json").read_text())
        results = card["results"]
        rows.append({"run_id": config["run_id"], "organism": config["organism"], "cycle": config["cycle"],
                     "status": "complete", "seconds": completion["seconds"], "run_dir": str(run_dir.relative_to(ROOT)),
                     "coverage": results["condition_level"], "gene_metrics": results["gene_level_conditions_where_wt_grows"],
                     "counts": results["counts"]})
        by_id[config["run_id"]] = load_run(run_dir)
    paired = []
    for org in ORGS:
        if not all(f"{org}_cycle{c}" in by_id for c in (6, 7)):
            continue
        a, b = (by_id[f"{org}_cycle{c}"] for c in (6, 7))
        report = compare_runs(a, b, n_boot=500, seed=0)
        report["organism"] = org
        st, ft = a["params"]["growth_threshold"], a["params"]["fitness_threshold"]
        ag = {str(c) for c, v in zip(a["conditions"], a["wt_growth"]) if v >= st}
        bg = {str(c) for c, v in zip(b["conditions"], b["wt_growth"]) if v >= st}
        report["conditions_gained"] = sorted(bg - ag); report["conditions_lost"] = sorted(ag - bg)
        sa, sb, fit, genes, conditions, _ = aligned_pairs(a, b)
        changed = np.isfinite(fit) & ((sa >= st) != (sb >= st))
        changes = [{"browser_gene": genes[i], "condition": conditions[j], "cycle6_growth": sa[i,j],
                    "cycle7_growth": sb[i,j], "fitness": fit[i,j],
                    "change": "to_growth" if sb[i,j] >= st else "to_no_growth",
                    "cycle7_agrees_with_threshold": bool((sb[i,j] >= st) == (fit[i,j] >= ft))}
                   for i, j in zip(*np.where(changed))]
        report["n_prediction_changes_on_shared_pairs"] = len(changes)
        report["n_changes_agree_with_fitness_threshold"] = sum(c["cycle7_agrees_with_threshold"] for c in changes)
        write_json(out / "prediction_changes" / f"{org}_cycle7_minus_cycle6.json", changes)
        paired.append(report)
    summary = {"manifest_sha256": sha(out / "prespecified_manifest.json"), "generated_utc": now(),
               "evaluation_role": "retrospective_development", "runs": rows, "paired_comparisons": paired,
               "interpretation": "All prespecified results are retained. Confidence intervals condition on previously selected patches and fixed conditions; they do not establish independent predictive improvement."}
    write_json(out / "summary.json", summary)
    print(json.dumps(safe({"runs": [{"run_id": r["run_id"], "status": r["status"], "coverage": r.get("coverage"),
                                    "mcc": r.get("gene_metrics", {}).get("mcc")} for r in rows],
                           "paired": [{"organism": p["organism"], "delta": p["paired_mcc_difference_B_minus_A"],
                                       "gained": p["conditions_gained"], "lost": p["conditions_lost"]} for p in paired]}), indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/development_sprint_2026_09_06")
    parser.add_argument("--mode", choices=("prepare", "run", "summarize", "all"), default="all")
    args = parser.parse_args()
    out = args.output_dir.resolve()
    manifest = prepare(out)
    if args.mode in ("run", "all"):
        run_all(out, manifest)
    if args.mode in ("summarize", "all"):
        summarize(out, manifest)
    if args.mode in ("run", "all") and any(not (out / "runs" / c["run_id"] / "completion.json").exists()
                                            for c in manifest["configurations"]):
        raise SystemExit("At least one prespecified run is incomplete or failed; preserved attempts can be inspected and resumed")


if __name__ == "__main__":
    main()
