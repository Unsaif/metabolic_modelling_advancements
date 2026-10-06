"""Run the disease-ranking prediction matrix in parallel shards on one machine, with progress, then join them.

Usage, from the repository root:
  caffeinate -i python3 scripts/run_ranking_parallel.py [--shards N] [--threads T]

Each shard is scripts/run_wbm_iem_cross.py with the plan's settings (docs/studies/wbm-iem-ranking-plan.md):
  Harvey_1_03d --backend gurobi --order readout --warm concurrent --context protocol
  --extra-readouts data/iem/iem_ranking_extra_readouts_v0.2.txt --recheck-every 100 --out-suffix _ranking_v1
  --check-against results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json --quiet --shard K/N --threads T
Shard K computes readouts K, K+N, K+2N, ... of the panel. Every readout's LPs are the same as in a single run.
Defaults: N = cores // 3 (1 to 6) and T = cores // N. Logs go to results/wbm_iem/logs/.

Stop with Ctrl+C. Running the same command again resumes every shard where it stopped (N and T are taken from the
existing shard files). When every shard is done, the shards are joined into
results/wbm_iem/Harvey_1_03d_iem_cross_ranking_v1.json by scripts/merge_cross_shards.py.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import signal
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNNER = os.path.join(ROOT, "scripts", "run_wbm_iem_cross.py")
MERGER = os.path.join(ROOT, "scripts", "merge_cross_shards.py")
OUT = os.path.join(ROOT, "results", "wbm_iem")
LOGS = os.path.join(OUT, "logs")


def shard_files(model, suffix):
    pat = os.path.join(OUT, f"{model}_iem_cross{suffix}_shard*of*.json")
    found = {}
    for path in glob.glob(pat):
        m = re.search(r"_shard(\d+)of(\d+)\.json$", path)
        if m:
            found[(int(m.group(1)), int(m.group(2)))] = path
    return found


def progress(path):
    """(readouts done, readouts in the shard) from a shard file; (0, None) if it does not exist yet."""
    try:
        with open(path) as fh:
            recs = json.load(fh)
    except (OSError, ValueError):
        return 0, None
    if not recs:
        return 0, None
    # In readout order every IEM that was set up holds the same readouts at each checkpoint.
    done = max((len(r.get("readouts") or []) for r in recs), default=0)
    return done, recs[0]["provenance"].get("n_readouts")


def fmt(seconds):
    seconds = int(seconds)
    h, m = divmod(seconds // 60, 60)
    return f"{h}h{m:02d}m" if h else f"{m}m"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shards", type=int)
    ap.add_argument("--threads", type=int)
    ap.add_argument("--model", default="Harvey_1_03d")
    ap.add_argument("--out-suffix", default="_ranking_v1")
    ap.add_argument("--backend", default="gurobi")
    ap.add_argument("--warm", default="concurrent")
    ap.add_argument("--extra-readouts", default="data/iem/iem_ranking_extra_readouts_v0.2.txt")
    ap.add_argument("--check-against", default="results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json")
    ap.add_argument("--limit", type=int, default=0, help="testing only: first N IEMs")
    ap.add_argument("--readout-limit", type=int, default=0, help="testing only: first N readouts")
    ap.add_argument("--interval", type=int, default=120, help="seconds between progress lines")
    args = ap.parse_args()

    cores = os.cpu_count() or 4
    existing = shard_files(args.model, args.out_suffix)
    if existing:
        ns = {n for _, n in existing}
        if len(ns) != 1:
            sys.exit(f"shard files with different N exist: {sorted(existing.values())}")
        n_existing = ns.pop()
        if args.shards and args.shards != n_existing:
            sys.exit(f"shard files for N={n_existing} exist; resume with --shards {n_existing} (or move them away)")
        args.shards = n_existing
        with open(next(iter(existing.values()))) as fh:
            t_existing = json.load(fh)[0]["provenance"]["threads"]
        if args.threads and args.threads != t_existing:
            sys.exit(f"the existing shards used --threads {t_existing}; resume with that")
        args.threads = t_existing
        print(f"resuming {args.shards} shards with {args.threads} threads each", flush=True)
    n = args.shards or max(1, min(6, cores // 3))
    t = args.threads or max(1, cores // n)
    os.makedirs(LOGS, exist_ok=True)
    print(f"{cores} cores: {n} shards x {t} solver threads. Logs: {os.path.relpath(LOGS, ROOT)}/", flush=True)

    procs = {}
    for k in range(1, n + 1):
        cmd = [sys.executable, RUNNER, args.model, "--backend", args.backend, "--order", "readout", "--warm", args.warm,
               "--context", "protocol", "--extra-readouts", args.extra_readouts, "--recheck-every", "100",
               "--out-suffix", args.out_suffix, "--check-against", args.check_against, "--quiet",
               "--shard", f"{k}/{n}", "--threads", str(t)]
        if args.limit:
            cmd += ["--limit", str(args.limit)]
        if args.readout_limit:
            cmd += ["--readout-limit", str(args.readout_limit)]
        log = open(os.path.join(LOGS, f"{args.model}{args.out_suffix}_shard{k}of{n}.log"), "a")
        log.write(f"\n=== start {time.strftime('%Y-%m-%d %H:%M:%S')}: {' '.join(cmd)}\n")
        log.flush()
        procs[k] = (subprocess.Popen(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT), log)

    t0 = time.time()
    first_seen = None
    paths = {k: os.path.join(OUT, f"{args.model}_iem_cross{args.out_suffix}_shard{k}of{n}.json") for k in procs}
    start_done = None
    failed = {}
    try:
        while True:
            alive = [k for k, (p, _) in procs.items() if p.poll() is None]
            for k, (p, log) in procs.items():
                if p.poll() not in (None, 0) and k not in failed:
                    failed[k] = p.returncode
                    log.flush()
                    with open(log.name) as fh:
                        tail = fh.readlines()[-15:]
                    print(f"\nshard {k} stopped with exit code {p.returncode}. Last lines of its log:\n" + "".join(tail), flush=True)
            if failed and time.time() - t0 < 600 and alive:
                print("A shard failed during set-up; stopping the others. If Gurobi refused extra processes, "
                      "rerun with fewer, e.g. --shards 2.", flush=True)
                for k in alive:
                    procs[k][0].send_signal(signal.SIGINT)
                break
            done_total, total, known = 0, 0, True
            parts = []
            for k in sorted(procs):
                d, tot = progress(paths[k])
                done_total += d
                total += tot or 0
                known &= bool(tot)
                parts.append(f"{d}/{tot if tot else '?'}")
            stamp = time.strftime("%H:%M")
            if done_total and first_seen is None:
                first_seen, start_done = time.time(), done_total
            if not done_total:
                line = f"[{stamp}] setting up the diseases in each shard ({fmt(time.time() - t0)} so far)"
            else:
                rate = (done_total - start_done) / max(1.0, time.time() - first_seen)
                remaining = (total - done_total) if known else None
                eta = f"; about {fmt(remaining / rate)} to go" if rate > 0 and remaining else ""
                line = f"[{stamp}] {done_total}/{total if known else '?'} biomarkers done (shards: {', '.join(parts)}){eta}"
            print(line, flush=True)
            if not alive:
                break
            time.sleep(args.interval)
    except KeyboardInterrupt:
        print("\nstopping the shards (each keeps every biomarker it finished)...", flush=True)
        for p, _ in procs.values():
            if p.poll() is None:
                p.send_signal(signal.SIGINT)
        for p, _ in procs.values():
            try:
                p.wait(timeout=60)
            except subprocess.TimeoutExpired:
                p.kill()
        print("stopped. Run the same command again to resume.", flush=True)
        return
    for p, log in procs.values():
        p.wait()
        log.close()
    if failed or any(p.returncode != 0 for p, _ in procs.values()):
        print("not all shards finished; run the same command again to resume.", flush=True)
        return
    print("all shards finished; joining them...", flush=True)
    files = [paths[k] for k in sorted(paths)]
    rc = subprocess.call([sys.executable, MERGER, *files, "--check-against", args.check_against], cwd=ROOT)
    print("done." if rc == 0 else "the join failed; the shard files are kept. Tell Claude.", flush=True)


if __name__ == "__main__":
    main()
