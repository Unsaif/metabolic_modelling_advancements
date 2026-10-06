"""Independent provenance, freeze and time-order checks for IEM study v0.4. Read-only git; no gembench imports.

- results files: 57 records, statuses, one run fingerprint, fingerprint == sha256(json.dumps(provenance, sort_keys=True))
  (the recipe the v0.3 verifier used), senses, lp_bounds_sha256, model/protocol/constraint-input hashes vs files on disk;
  summary files agree; wall times (summary session_wall_s and sum of per-IEM time_s); inferred start = mtime - wall.
- freeze: my own recomputation of the manifest's content fingerprint (recipe read, not imported, from gembench/study.py),
  every frozen file's size and sha256 now and at the freeze commit, files unchanged from the freeze commit to HEAD.
- git: commit times of plan, freeze, first Harvetta v0.4 IEM result, first biomarker minimum, MATLAB run; ancestry.
- MATLAB job records: UTC times of both attempts, code hash vs the frozen job script.
Writes provenance_v04.json here.
"""
import hashlib
import json
import os
import subprocess
import sys
from collections import Counter, OrderedDict
from datetime import datetime, timedelta, timezone

from vcommon import ROOT, dump, load_json, sha256_file

RESULTS = OrderedDict([
    ("A_harvetta_max", ("results/wbm_iem/Harvetta_1_03d_iem_results_v0.4.json", "results/wbm_iem/Harvetta_1_03d_iem_summary_v0.4.json")),
    ("B_harvey_min", ("results/wbm_iem/Harvey_1_03d_iem_results_v0.4_min.json", "results/wbm_iem/Harvey_1_03d_iem_summary_v0.4_min.json")),
    ("C_harvetta_min", ("results/wbm_iem/Harvetta_1_03d_iem_results_v0.4_min.json", "results/wbm_iem/Harvetta_1_03d_iem_summary_v0.4_min.json")),
    ("v0.3_harvey_max", ("results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json", "results/wbm_iem/Harvey_1_03d_iem_summary_v0.3.json")),
])
MODEL_FILES = {"Harvey": "external/COBRA.models/mat/Harvey_1_03d.mat", "Harvetta": "external/COBRA.models/mat/Harvetta_1_03d.mat"}
MANIFEST = "results/study_freezes/wbm_iem_v0.4_plan.json"


def git(*args):
    return subprocess.run(["git", "-C", ROOT, *args], capture_output=True, text=True, check=False)


def commits(path):
    out = git("log", "--format=%H|%cI|%aI|%s", "--", path).stdout.strip().splitlines()
    return [dict(zip(("hash", "committed", "authored", "subject"), l.split("|", 3))) for l in out]


def anc(a, b):
    return git("merge-base", "--is-ancestor", a, b).returncode == 0


def mtime_utc(rel):
    return datetime.fromtimestamp(os.stat(os.path.join(ROOT, rel)).st_mtime, timezone.utc)


def main():
    rep = OrderedDict()
    disk = {"model_" + k: sha256_file(v) for k, v in MODEL_FILES.items()}
    disk["protocol"] = sha256_file("data/iem/iem_protocol_v0.2.json")
    disk["constraint_inputs"] = sha256_file("data/iem/wbm_constraint_inputs_v0.3.json")
    rep["disk_sha256"] = disk
    v03_lp = None
    runs = OrderedDict()
    for label, (path, spath) in RESULTS.items():
        recs = load_json(path)
        summ = load_json(spath)
        fps = sorted({r["run_fingerprint"] for r in recs})
        provs = {json.dumps(r["provenance"], sort_keys=True) for r in recs}
        prov = recs[0]["provenance"]
        recomputed = hashlib.sha256(json.dumps(prov, sort_keys=True).encode()).hexdigest()
        ms = prov["model_setup"]
        model = "Harvetta" if "Harvetta" in path else "Harvey"
        sum_time = sum(r["time_s"] for r in recs)
        end = mtime_utc(spath)
        start = end - timedelta(seconds=summ["session_wall_s"])
        runs[label] = OrderedDict(
            n_records=len(recs), iems=len({r["iem"] for r in recs}), statuses=dict(Counter(r["status"] for r in recs)),
            notes=sorted({json.dumps(r["notes"]) for r in recs}),
            n_fingerprints=len(fps), fingerprint=fps[0] if len(fps) == 1 else fps, n_provenance_blocks=len(provs),
            fingerprint_equals_sha256_of_provenance=recomputed == fps[0],
            summary_fingerprint_equal=summ["run_fingerprint"] == fps[0],
            summary_provenance_equal=json.dumps(summ["provenance"], sort_keys=True) == json.dumps(prov, sort_keys=True),
            senses=prov.get("senses", "(field absent)"), provenance_keys=sorted(prov),
            model_sha256_matches_file=prov["model_sha256"] == disk["model_" + model],
            protocol_sha256_matches_file=prov["protocol_sha256"] == disk["protocol"],
            constraint_inputs_sha256_matches_file=ms.get("constraint_inputs_sha256") == disk["constraint_inputs"],
            solver=prov["solver"], method=prov["method"], feas_tol=prov["feas_tol"], opt_tol=prov["opt_tol"],
            time_limit=prov["time_limit_per_solve_s"], min_flux_healthy=prov["min_flux_healthy"],
            bile_duct=prov["bile_duct"], model_setup=ms["model_setup"], toolbox_commit=ms["toolbox_commit"],
            n_lb_changed=ms["n_lb_changed"], n_ub_changed=ms["n_ub_changed"], sex=ms["parameters"]["sex"],
            global_constraints=prov["global_constraints"], lp_bounds_sha256=prov["lp_bounds_sha256"],
            sum_time_s=sum_time, sum_time_h=sum_time / 3600, session_wall_s=summ.get("session_wall_s"),
            session_wall_h=(summ.get("session_wall_s") or 0) / 3600, session_solves=summ.get("session_solves"),
            total_solves=summ.get("total_solves_in_recorded_attempts"), sum_n_solves=sum(r["n_solves"] for r in recs),
            summary_mtime_utc=end.isoformat(), inferred_start_utc=start.isoformat(),
            summary_counts=dict(n_scored=summ.get("n_biomarkers_scored"), n_correct=summ.get("n_correct"),
                                n_optimal_minima=summ.get("n_biomarkers_with_optimal_minima")))
        if label == "v0.3_harvey_max":
            v03_lp = prov["lp_bounds_sha256"]
    runs["B_harvey_min"]["lp_bounds_equals_v0.3"] = runs["B_harvey_min"]["lp_bounds_sha256"] == v03_lp
    runs["C_harvetta_min"]["lp_bounds_equals_A"] = runs["C_harvetta_min"]["lp_bounds_sha256"] == runs["A_harvetta_max"]["lp_bounds_sha256"]
    # A vs v0.3: same provenance keys (plan: 'a max-only run has the same provenance fields as v0.3')
    runs["A_harvetta_max"]["same_provenance_keys_as_v0.3"] = runs["A_harvetta_max"]["provenance_keys"] == runs["v0.3_harvey_max"]["provenance_keys"]
    a_bm_keys = sorted({k for r in load_json(RESULTS["A_harvetta_max"][0]) for b in r["biomarkers"] for k in b})
    v3_bm_keys = sorted({k for r in load_json(RESULTS["v0.3_harvey_max"][0]) for b in r["biomarkers"] for k in b})
    runs["A_harvetta_max"]["biomarker_keys_extra_vs_v0.3"] = sorted(set(a_bm_keys) - set(v3_bm_keys))
    rep["runs"] = runs
    # ---- freeze: my own fingerprint
    man = load_json(MANIFEST)
    files = []
    for rec in man["files"]:
        p = os.path.join(ROOT, rec["path"])
        files.append(dict(path=rec["path"], size_ok=os.path.getsize(p) == rec["size_bytes"], sha_ok=sha256_file(rec["path"]) == rec["sha256"]))
    plan = json.loads(json.dumps(man["plan"]))
    plan["paths"] = sorted(plan["paths"])
    content = {"schema_version": man["schema_version"], "manifest_type": man["manifest_type"], "plan": plan,
               "files": sorted(man["files"], key=lambda x: x["path"])}
    fp = hashlib.sha256(json.dumps(content, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                                   allow_nan=False).encode("utf-8")).hexdigest()
    plan_json = load_json("data/studies/wbm_iem_v0.4_plan.json")
    plan_json_norm = json.loads(json.dumps(plan_json))
    plan_json_norm["paths"] = sorted(plan_json_norm["paths"])
    freeze_c = commits(MANIFEST)[-1]["hash"]
    rep["freeze"] = OrderedDict(
        created_at=man["created_at"], repository=man["repository"], n_files=len(files),
        n_size_ok=sum(f["size_ok"] for f in files), n_sha_ok=sum(f["sha_ok"] for f in files),
        failing=[f for f in files if not (f["size_ok"] and f["sha_ok"])],
        recorded=man["content_fingerprint"], recomputed=fp, match=fp == man["content_fingerprint"],
        plan_json_equals_manifest_plan=plan_json_norm == man["plan"],
        frozen_model_sha_vs_plan=dict(Harvetta=plan_json["models"]["Harvetta_1_03d"]["sha256"] == disk["model_Harvetta"],
                                      Harvey=plan_json["models"]["Harvey_1_03d"]["sha256"] == disk["model_Harvey"]),
        frozen_paths_changed_freeze_commit_to_HEAD=git("diff", "--name-only", freeze_c, "HEAD", "--", *man["plan"]["paths"]).stdout.split(),
        frozen_paths_changed_in_worktree=git("diff", "--name-only", "HEAD", "--", *man["plan"]["paths"]).stdout.split(),
        plan_json_changed_since_plan_commit=git("log", "--format=%h", "--", "data/studies/wbm_iem_v0.4_plan.json").stdout.split(),
        freeze_verify_tool=load_json("results/wbm_iem/independent_verification_v0.4/freeze_verify_output.json"))
    # ---- git time order
    paths = OrderedDict([
        ("plan_md", "docs/studies/wbm-iem-v0.4-plan.md"), ("plan_json", "data/studies/wbm_iem_v0.4_plan.json"),
        ("freeze_manifest", MANIFEST),
        ("A_results", RESULTS["A_harvetta_max"][0]), ("B_results", RESULTS["B_harvey_min"][0]), ("C_results", RESULTS["C_harvetta_min"][0]),
        ("M_matlab_results", "results/wbm_iem/matlab_reference/harvetta_step2/matlab_iem_results_Harvetta_1_03d.mat"),
        ("harvetta_setup_bounds", "results/wbm_iem/matlab_reference/harvetta_step1/matlab_setup_bounds_Harvetta_1_03d.mat"),
        ("harvetta_port_check", "results/wbm_iem/Harvetta_1_03d_constraint_port_check_v0.3.json"),
        ("error_anatomy_script", "scripts/iem_error_anatomy.py"), ("v04_report_script", "scripts/v04_report.py"),
        ("v04_report_json", "results/wbm_iem/v0.4_report.json"),
        ("results_note", "docs/studies/wbm-iem-v0.4-results.md"), ("paper2", "docs/paper/paper2-draft.md")])
    g = OrderedDict((k, commits(p)) for k, p in paths.items())
    plan_c, freeze_c = g["plan_md"][-1]["hash"], g["freeze_manifest"][-1]["hash"]
    first = {k: g[k][-1]["hash"] for k in ("A_results", "B_results", "C_results", "M_matlab_results")}
    # first appearance of any path containing Harvetta IEM results or minima, across all refs
    allpaths = git("log", "--all", "--format=@@%H|%cI", "--name-only").stdout.splitlines()
    hits, cur = [], None
    for line in allpaths:
        if line.startswith("@@"):
            cur = line[2:]
        elif line and (("Harvetta" in line and "iem_results" in line) or ("iem_results" in line and "_min" in line)
                       or "matlab_iem_results_Harvetta" in line):
            hits.append((cur, line))
    rep["git"] = OrderedDict(
        commits={k: v for k, v in g.items()},
        plan_commit=plan_c, freeze_commit=freeze_c,
        plan_ancestor_of_first_results={k: anc(plan_c, c) and c != plan_c for k, c in first.items()},
        freeze_ancestor_of_first_results={k: anc(freeze_c, c) and c != freeze_c for k, c in first.items()},
        first_appearance_any_ref=sorted(hits, key=lambda x: x[0].split("|")[1])[:6],
        n_refs_commits=dict(all=len(git("rev-list", "--all").stdout.split()), head=len(git("rev-list", "HEAD").stdout.split())),
        head=git("rev-parse", "HEAD").stdout.strip(), status_porcelain=git("status", "--porcelain").stdout.splitlines())
    # ---- MATLAB job records
    q1 = load_json("results/wbm_iem/matlab_reference/harvetta_step2/queue_result_step2_first_attempt_failed.json")
    q2 = load_json("results/wbm_iem/matlab_reference/harvetta_step2/queue_result_step2b.json")
    utc = lambda t: datetime.fromtimestamp(t, timezone.utc).isoformat()
    frozen_job = next(f for f in man["files"] if f["path"] == "tools/matlab/runner/iem_reference_harvetta_step2_runiem.m")
    rep["matlab_jobs"] = OrderedDict(
        first=dict(script=q1["script"], status=q1["status"], submitted=utc(q1["submitted_at"]), finished=utc(q1["finished_at"]),
                   hours=(q1["finished_at"] - q1["submitted_at"]) / 3600, error=q1.get("error"), code_sha256=q1["code_sha256"],
                   automatic_retry=q1.get("automatic_retry"), may_have_executed=q1.get("may_have_executed")),
        second=dict(script=q2["script"], status=q2["status"], submitted=utc(q2["submitted_at"]), finished=utc(q2["finished_at"]),
                    hours=(q2["finished_at"] - q2["submitted_at"]) / 3600, elapsed_h=q2["runner_result"]["elapsed_seconds"] / 3600,
                    code_sha256=q2["code_sha256"], errors=q2["runner_result"]["error"]),
        code_equals_frozen_job_script=q1["code_sha256"] == q2["code_sha256"] == frozen_job["sha256"],
        agent_job_m_sha256=sha256_file("results/wbm_iem/matlab_reference/harvetta_step2/agent_job.m"),
        run_info=open(os.path.join(ROOT, "results/wbm_iem/matlab_reference/harvetta_step2/run_info_harvetta_step2.txt")).read().splitlines())
    dump("provenance_v04.json", rep)
    print(json.dumps(rep, indent=1, default=str)[:20000])


if __name__ == "__main__":
    sys.exit(main())
