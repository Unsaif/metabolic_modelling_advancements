"""Independent provenance and time-order checks (claim 6). Read-only git commands; no gembench imports.

Writes provenance.json next to this script.
"""
import hashlib
import json
import os
import subprocess
from collections import OrderedDict

ROOT = "/home/claude/mma"
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = OrderedDict([
    ("v0.2", "results/wbm_iem/Harvey_1_03d_iem_results_v0.2.json"),
    ("v0.2b", "results/wbm_iem/Harvey_1_03d_iem_results_v0.2b.json"),
    ("v0.3", "results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json"),
])
MANIFEST = "results/study_freezes/wbm_iem_v0.3_plan.json"
RECOVERY = "results/study_freezes/wbm_iem_v0.3_plan_recovery.json"


def sha256_file(path):
    h = hashlib.sha256()
    with open(os.path.join(ROOT, path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git(*args):
    return subprocess.run(["git", "-C", ROOT, *args], capture_output=True, text=True, check=False)


def commit_info(rev):
    out = git("show", "-s", "--format=%H%x09%cI%x09%aI%x09%s", rev).stdout.strip().split("\t")
    return {"rev": rev, "hash": out[0], "committed": out[1], "authored": out[2], "subject": out[3]}


def is_ancestor(a, b):
    return git("merge-base", "--is-ancestor", a, b).returncode == 0


def main():
    rep = OrderedDict()
    # 1. fingerprints and provenance in the results files
    runs = OrderedDict()
    for label, path in RUNS.items():
        recs = json.load(open(os.path.join(ROOT, path)))
        fps = sorted({r["run_fingerprint"] for r in recs})
        provs = {json.dumps(r["provenance"], sort_keys=True) for r in recs}
        prov = recs[0]["provenance"]
        recomputed = hashlib.sha256(json.dumps(prov, sort_keys=True).encode()).hexdigest()
        ms = prov.get("model_setup", {})
        runs[label] = OrderedDict(
            n_records=len(recs), n_distinct_fingerprints=len(fps), fingerprint=fps[0] if len(fps) == 1 else fps,
            n_distinct_provenance_blocks=len(provs),
            fingerprint_equals_sha256_of_provenance_json=recomputed == fps[0],
            model_setup=ms.get("model_setup") if isinstance(ms, dict) else ms, bile_duct=prov.get("bile_duct"),
            solver=prov.get("solver"), method=prov.get("method"), feas_tol=prov.get("feas_tol"), opt_tol=prov.get("opt_tol"),
            model_sha256=prov.get("model_sha256"), protocol_sha256=prov.get("protocol_sha256"),
            engine_version=prov.get("engine_version"), global_constraints=prov.get("global_constraints"),
            lp_bounds_sha256=prov.get("lp_bounds_sha256"),
            constraint_inputs_sha256=ms.get("constraint_inputs_sha256") if isinstance(ms, dict) else None,
            toolbox_commit=ms.get("toolbox_commit") if isinstance(ms, dict) else None,
            n_lb_changed=ms.get("n_lb_changed") if isinstance(ms, dict) else None,
            n_ub_changed=ms.get("n_ub_changed") if isinstance(ms, dict) else None,
            physiology_and_diet_ported=prov.get("physiology_and_diet_ported"),
            statuses=sorted({r["status"] for r in recs}))
    rep["runs"] = runs
    # 2. hashes of the inputs on disk
    rep["file_sha256"] = {
        "model": sha256_file("external/COBRA.models/mat/Harvey_1_03d.mat"),
        "protocol": sha256_file("data/iem/iem_protocol_v0.2.json"),
        "constraint_inputs": sha256_file("data/iem/wbm_constraint_inputs_v0.3.json")}
    try:
        import highspy
        rep["installed_highs"] = highspy.Highs().version()
    except Exception as e:  # pragma: no cover
        rep["installed_highs"] = repr(e)
    # 3. freeze manifest, recomputed with my own code (recipe read from gembench/study.py: sha256 of
    #    json.dumps({schema_version, manifest_type, plan (paths sorted), files sorted by path}, sort_keys,
    #    separators (',', ':'), ensure_ascii False))
    man = json.load(open(os.path.join(ROOT, MANIFEST)))
    files = []
    for rec in man["files"]:
        p = os.path.join(ROOT, rec["path"])
        files.append({"path": rec["path"], "size_ok": os.path.getsize(p) == rec["size_bytes"],
                      "sha256_ok": sha256_file(rec["path"]) == rec["sha256"]})
    plan = json.loads(json.dumps(man["plan"]))
    plan["paths"] = sorted(plan["paths"])
    content = {"schema_version": man["schema_version"], "manifest_type": man["manifest_type"], "plan": plan,
               "files": sorted(man["files"], key=lambda x: x["path"])}
    fp = hashlib.sha256(json.dumps(content, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                                   allow_nan=False).encode("utf-8")).hexdigest()
    recovery = json.load(open(os.path.join(ROOT, RECOVERY)))
    rep["freeze"] = OrderedDict(
        created_at=man["created_at"], repository=man["repository"], n_files=len(files),
        n_size_ok=sum(f["size_ok"] for f in files), n_sha256_ok=sum(f["sha256_ok"] for f in files),
        failing=[f for f in files if not (f["size_ok"] and f["sha256_ok"])],
        content_fingerprint_recorded=man["content_fingerprint"], content_fingerprint_recomputed=fp,
        recovery_original_fingerprint=recovery["original_freeze"]["content_fingerprint"],
        matches=fp == man["content_fingerprint"] == recovery["original_freeze"]["content_fingerprint"],
        plan_runs=man["plan"].get("runs"), plan_model_sha256=man["plan"].get("model_sha256"),
        frozen_constraint_inputs_sha256=next(r["sha256"] for r in man["files"] if r["path"] == "data/iem/wbm_constraint_inputs_v0.3.json"),
        frozen_protocol_sha256=next(r["sha256"] for r in man["files"] if r["path"] == "data/iem/iem_protocol_v0.2.json"),
        frozen_v0_2_results_sha256=next(r["sha256"] for r in man["files"] if r["path"] == RUNS["v0.2"]),
        v0_2_results_sha256_now=sha256_file(RUNS["v0.2"]))
    # 4. git time order
    commits = {
        "plan": git("log", "--format=%h", "--", "docs/studies/wbm-iem-v0.3-plan.md").stdout.split(),
        "freeze_manifest": git("log", "--format=%h", "--", MANIFEST).stdout.split(),
        "v0.3_results": git("log", "--format=%h", "--", RUNS["v0.3"]).stdout.split(),
        "v0.2b_results": git("log", "--format=%h", "--", RUNS["v0.2b"]).stdout.split(),
        "v0.2_results": git("log", "--format=%h", "--", RUNS["v0.2"]).stdout.split(),
        "matlab_step2": git("log", "--format=%h", "--",
                            "results/wbm_iem/matlab_reference/step2/matlab_iem_results_Harvey_1_03d.mat").stdout.split(),
    }
    info = {k: [commit_info(c) for c in v] for k, v in commits.items()}
    plan_c, freeze_c = commits["plan"][-1], commits["freeze_manifest"][-1]
    first_res = [commits["v0.3_results"][-1], commits["v0.2b_results"][-1]]
    rep["git"] = OrderedDict(
        commits=info,
        plan_commit_is_ancestor_of_first_results=[is_ancestor(plan_c, c) and c != plan_c for c in first_res],
        freeze_commit_is_ancestor_of_first_results=[is_ancestor(freeze_c, c) and c != freeze_c for c in first_res],
        freeze_manifest_created_at=man["created_at"],
        frozen_paths_changed_between_freeze_commit_and_HEAD=git("diff", "--name-only", freeze_c, "HEAD", "--",
                                                                 *man["plan"]["paths"]).stdout.split(),
        frozen_paths_changed_in_working_tree=git("diff", "--name-only", "HEAD", "--", *man["plan"]["paths"]).stdout.split(),
        results_changed_in_working_tree=git("diff", "--name-only", "HEAD", "--", RUNS["v0.3"], RUNS["v0.2b"], RUNS["v0.2"]).stdout.split(),
        status_porcelain=git("status", "--porcelain").stdout.splitlines(),
        head=commit_info("HEAD"))
    # 5. partial checkpoint: were the records committed at the checkpoint kept unchanged in the final files?
    ckpt = {}
    for label in ("v0.3", "v0.2b"):
        c = commits[f"{label}_results"][-1]
        old = json.loads(git("show", f"{c}:{RUNS[label]}").stdout)
        new = {(r["iem"], r["call_index"]): r for r in json.load(open(os.path.join(ROOT, RUNS[label])))}
        same = sum(json.dumps(r, sort_keys=True) == json.dumps(new[(r["iem"], r["call_index"])], sort_keys=True) for r in old)
        ckpt[label] = {"commit": c, "n_records_at_checkpoint": len(old), "n_identical_in_final": same,
                       "fingerprints_at_checkpoint": sorted({r["run_fingerprint"] for r in old}),
                       "iems_at_checkpoint": [r["iem"] for r in old]}
    rep["checkpoint"] = ckpt
    with open(os.path.join(HERE, "provenance.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
        fh.write("\n")
    print(json.dumps(rep, indent=1)[:12000])


if __name__ == "__main__":
    main()
