#!/usr/bin/env python3
"""Independent check of Supplementary Table S1 and the Methods paragraph "Locking, time order and verification".

Uses only read-only git commands (git log, git show, git cat-file) and the JSON records; no gembench or scripts/
imports. External-record times are taken from the task (earlier verified documents), not from the repository.
Writes check_timeline_and_locks.json next to this file.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone

ROOT = "/home/claude/mma"
OUT = os.path.join(ROOT, "results", "transfer_v1", "reference_models", "independent_check")

MANIFESTS = ["transfer_v1_method", "transfer_v1_inputs_panel_A", "transfer_v1_replication_plan", "transfer_v1_inputs_panel_B"]
EXTERNAL = {"transfer_v1_method": "2026-10-03T09:41:01Z", "transfer_v1_inputs_panel_A": "2026-10-03T10:11:50Z",
            "transfer_v1_replication_plan": "2026-10-03T12:31:17Z", "transfer_v1_inputs_panel_B": "2026-10-03T13:06:45Z"}
TABLE_S1 = {"transfer_v1_method": "09:40:45", "transfer_v1_inputs_panel_A": "10:11:43",
            "transfer_v1_replication_plan": "12:31:09", "transfer_v1_inputs_panel_B": "13:06:21",
            "panel_A_download": ("10:12:18", "10:12:40"), "panel_B_download": ("13:07:00", "13:07:20")}
DOWNLOADS = {"panel_A_download": "data/fitness_browser_panel/outcomes_download_2026-10-03.json",
             "panel_B_download": "data/fitness_browser_panel/outcomes_download_panel_B_2026-10-03.json"}


def git(*args, binary=False):
    r = subprocess.run(["git", "-C", ROOT, *args], capture_output=True)
    if r.returncode != 0:
        return None
    return r.stdout if binary else r.stdout.decode().strip()


def ts(s):
    s = s.replace("Z", "+00:00")
    return datetime.fromisoformat(s).astimezone(timezone.utc)


def sha(b):
    return hashlib.sha256(b).hexdigest()


def first_commit(path):
    out = git("log", "--diff-filter=A", "--format=%H %cI", "--", path)
    return out.splitlines()[-1].split() if out else None


def verify_manifest(man, commit):
    """Hashes of every listed file at `commit` (git blob) and at HEAD/working tree."""
    res = {"n_files": len(man["files"]), "at_commit_mismatch": [], "at_commit_untracked": [], "at_head_mismatch": [],
           "worktree_mismatch": []}
    for f in man["files"]:
        p, h = f["path"], f["sha256"]
        blob = git("show", f"{commit}:{p}", binary=True)
        if blob is None:
            res["at_commit_untracked"].append(p)
            wt = os.path.join(ROOT, p)
            if os.path.exists(wt) and sha(open(wt, "rb").read()) != h:
                res["worktree_mismatch"].append(p)
        elif sha(blob) != h:
            res["at_commit_mismatch"].append(p)
        hb = git("show", f"HEAD:{p}", binary=True)
        if hb is None:
            wt = os.path.join(ROOT, p)
            if not os.path.exists(wt) or sha(open(wt, "rb").read()) != h:
                res["at_head_mismatch"].append(p + " (untracked/worktree)")
        elif sha(hb) != h:
            res["at_head_mismatch"].append(p)
    return res


def main():
    out = {"manifests": {}, "downloads": {}, "runs": {}}
    for m in MANIFESTS:
        path = f"results/study_freezes/{m}.json"
        man = json.load(open(os.path.join(ROOT, path)))
        commit, ctime = first_commit(path)
        created = ts(man["created_at"])
        ext = ts(EXTERNAL[m])
        rec = {"created_at": man["created_at"], "commit": commit, "commit_time": ctime,
               "table_s1": TABLE_S1[m], "table_s1_matches_created_at": man["created_at"][11:19] == TABLE_S1[m],
               "external": EXTERNAL[m],
               "external_minus_created_s": (ext - created).total_seconds(),
               "external_minus_commit_s": (ext - ts(ctime)).total_seconds(),
               "commit_minus_created_s": (ts(ctime) - created).total_seconds()}
        rec["verify"] = verify_manifest(man, commit)
        out["manifests"][m] = rec

    # filter_report.json: earlier entries unchanged between the replication-plan commit and HEAD?
    rp = out["manifests"]["transfer_v1_replication_plan"]
    changed = [p for p in rp["verify"]["at_head_mismatch"]]
    detail = {}
    for p in changed:
        old = git("show", f"{rp['commit']}:{p}")
        new = git("show", f"HEAD:{p}")
        try:
            jo, jn = json.loads(old), json.loads(new)
        except Exception:  # noqa: BLE001
            detail[p] = "not JSON"
            continue

        def flat(x, pre=""):
            if isinstance(x, dict):
                d = {}
                for k, v in x.items():
                    d.update(flat(v, f"{pre}/{k}"))
                return d
            return {pre: json.dumps(x, sort_keys=True)}
        fo, fn = flat(jo), flat(jn)
        detail[p] = {"keys_old": len(fo), "keys_new": len(fn),
                     "old_keys_changed": [k for k in fo if k in fn and fo[k] != fn[k]][:20],
                     "old_keys_removed": [k for k in fo if k not in fn][:20],
                     "new_top_level_keys": sorted(set(jn) - set(jo)) if isinstance(jn, dict) else None,
                     "commits_touching_after_plan": git("log", "--format=%h %cI %s", f"{rp['commit']}..HEAD", "--", p)}
    out["replication_plan_head_diff"] = detail

    for k, f in DOWNLOADS.items():
        d = json.load(open(os.path.join(ROOT, f)))
        starts = [ts(v["fetch_started"]) for v in d["files"].values()]
        ends = [ts(v["fetch_finished"]) for v in d["files"].values()]
        hash_ok = {p: sha(open(os.path.join(ROOT, p), "rb").read()) == v["sha256"] for p, v in d["files"].items()}
        commit, ctime = first_commit(f)
        out["downloads"][k] = {"n_files": len(d["files"]), "first_start": min(starts).isoformat(), "last_finish": max(ends).isoformat(),
                               "files_hash_ok": all(hash_ok.values()), "hash_failures": [p for p, ok in hash_ok.items() if not ok],
                               "record_commit_time": ctime, "table_s1": TABLE_S1[k],
                               "fit_files_first_commit": sorted({(first_commit(p) or ["?", "?"])[1] for p in d["files"]})}

    # runs: arm cards of each panel, created after the corresponding download, one card per arm directory
    for phase, dk in (("evaluation_panel_A", "panel_A_download"), ("evaluation_panel_B", "panel_B_download")):
        created = []
        for c in glob.glob(os.path.join(ROOT, "results", "transfer_v1", phase, "*", "*", "card.json")):
            created.append(ts(json.load(open(c))["created"]))
        first_dl = ts(out["downloads"][dk]["first_start"])
        out["runs"][phase] = {"n_arm_cards": len(created), "first": min(created).isoformat(), "last": max(created).isoformat(),
                              "all_after_first_download": all(t > first_dl for t in created),
                              "all_after_last_download": all(t > ts(out["downloads"][dk]["last_finish"]) for t in created)}
    dev = [ts(json.load(open(c))["created"]) for c in glob.glob(os.path.join(ROOT, "results/transfer_v1/development/*/*/card.json"))]
    out["runs"]["development"] = {"n_arm_cards": len(dev), "first": min(dev).isoformat(), "last": max(dev).isoformat()}

    # order checks
    m = out["manifests"]
    out["order"] = {
        "method_lock_before_panel_A_lock": ts(m["transfer_v1_method"]["created_at"]) < ts(m["transfer_v1_inputs_panel_A"]["created_at"]),
        "panel_A_lock_before_panel_A_download": ts(m["transfer_v1_inputs_panel_A"]["commit_time"]) < ts(out["downloads"]["panel_A_download"]["first_start"]),
        "plan_lock_before_panel_B_lock": ts(m["transfer_v1_replication_plan"]["created_at"]) < ts(m["transfer_v1_inputs_panel_B"]["created_at"]),
        "panel_B_lock_before_panel_B_download": ts(m["transfer_v1_inputs_panel_B"]["commit_time"]) < ts(out["downloads"]["panel_B_download"]["first_start"]),
        "external_within_30s_of_commit": all(0 <= v["external_minus_commit_s"] < 30 for v in m.values()),
    }
    with open(os.path.join(OUT, "check_timeline_and_locks.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
        fh.write("\n")
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
