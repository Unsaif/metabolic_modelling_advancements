"""Healthy maxima at several pin levels (the Toolbox's minRxnsFluxHealthy), for a set of IEMs and every panel readout.

Usage, from the repository root (Tim's Mac, Gurobi):
  caffeinate -i python3 scripts/run_wbm_iem_pin_sweep.py

Defaults: the development IEMs of data/iem/iem_ranking_split_v1.json, the 187-readout panel of the ranking matrix,
alpha = 1, 0.5, 0.1, 0.01, 0.001 and 0, Gurobi with concurrent warm starts, and N parallel shards (N = cores // 3, 1 to 6;
cores // N threads each). Progress lines every two minutes; Ctrl+C stops, and the same command resumes. When every
shard is done they are joined into results/wbm_iem/pin_sweep/Harvey_1_03d_pin_sweep_development_v1.json.

What is computed (docs/studies/wbm-iem-healthy-reference-dev-plan.md):
  - The LPs are rebuilt exactly as for the ranking matrix: same model, setup, global constraints, panel sinks and
    each IEM's tweaks and own sinks, with the IEM set-ups taken from the matrix's records (no set-up solve).
  - For each readout and IEM, the healthy state is solved at each alpha in turn: the pin is
    sum(IEM fluxes) >= alpha * v_max, truncated to six decimals as in checkIEM_WBM; at alpha = 1 it is the matrix's pin,
    so those values check the rebuild against the matrix. At alpha = 0 the pin is sum >= 0, which the disease state
    meets, so no readout can be higher in the disease state than there: those values show what the block removes.
  - Between the alpha levels only the pin's lower bound changes, so those solves warm-start.
  - The first solve of each readout is barrier with crossover; the rest use --warm. A warm solve that does not end
    optimal is solved again by barrier, and both are recorded. About 1 in --recheck-every warm solves is re-solved
    from scratch by barrier as a check.
  - The disease state does not depend on alpha; its values are the matrix's.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
import signal
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from gembench import wbm as W  # noqa: E402
from gembench import wbm_iem as I  # noqa: E402
from gembench import wbm_constraints as WC  # noqa: E402

spec = importlib.util.spec_from_file_location("run_wbm_iem_cross", os.path.join(ROOT, "scripts", "run_wbm_iem_cross.py"))
X = importlib.util.module_from_spec(spec)
spec.loader.exec_module(X)

ENGINE_VERSION = "iem-pin-sweep-v0.1"
OUT = os.path.join(ROOT, "results", "wbm_iem", "pin_sweep")
DEFAULT_ALPHAS = "1,0.5,0.1,0.01,0.001,0"


def pin_for(vmax, alpha, stored_pin=None):
    """The healthy pin at fraction alpha of v_max, truncated to six decimals towards zero as checkIEM_WBM does.
    At alpha = 1 the matrix's own pin is used, so the alpha = 1 LPs are exactly the matrix's."""
    if alpha == 1 and stored_pin is not None:
        return stored_pin
    v = alpha * vmax
    return math.floor(v * 1e6) / 1e6 if v > 0 else math.ceil(v * 1e6) / 1e6


def sweep_readout(hw, mover, setups, rid, alphas, warm, first, recheck_every, log=False):
    """Healthy maxima of readout rid for every set-up IEM at every alpha. Returns {iem: [entry per alpha]}."""
    out = {}
    if rid not in hw.rxn_pos:
        for s in setups:
            out[s.iem] = [{"alpha": a, "status": "absent", "value": None} for a in alphas]
        return out
    rcol = hw.rxn_pos[rid]
    hw.set_objective({rcol: 1.0}, "max")
    is_first = True
    for s in setups:
        entries = []
        lo_matrix = s.lo
        try:
            for a in alphas:
                s.lo = pin_for(s.rec["vmax_healthy"], a, lo_matrix)
                target, row_bounds = X.target_bounds(s, "healthy", rcol, mover.base_lb, mover.base_ub)
                mover.move(target, s.row, row_bounds)
                method = "ipm" if (is_first and first == "ipm") or warm == "ipm" else warm
                st, f, dt, info, fb = X.solve_logged(hw, method, f"{rid} {s.iem} alpha={a}", log)
                is_first = False
                entry = {"alpha": a, "pin": s.lo, "status": st, "time_s": round(dt, 3), "method": method,
                         "fallback": fb, "simplex_iterations": info.get("simplex_iterations"),
                         "ipm_iterations": info.get("ipm_iterations")}
                ok = st == "Optimal" and math.isfinite(f)
                entry["value"] = (0.0 if abs(f) <= X.TOL else f) if ok else None
                if method != "ipm" and X.rechecked(s.iem, f"healthy@{a}", rid, recheck_every):
                    hw.fresh()   # discard the solution, or the barrier solve would return it unchanged
                    st2, f2, dt2, rinfo, _ = X.solve_logged(hw, "ipm", f"{rid} {s.iem} alpha={a} recheck", log)
                    entry["recheck"] = {"status": st2, "value": f2, "time_s": round(dt2, 3),
                                        "ipm_iterations": rinfo.get("ipm_iterations"),
                                        "abs_diff": abs(f2 - f) if math.isfinite(f) and math.isfinite(f2) else None}
                entries.append(entry)
        finally:
            s.lo = lo_matrix
        out[s.iem] = entries
    return out


def build(matrix_path, extra_readouts, protocol_path, backend, threads, iems):
    """Rebuild the matrix's LP and the chosen IEMs' set-ups; checks every input against the matrix's provenance."""
    with open(matrix_path) as fh:
        recs = json.load(fh)
    prov = recs[0]["provenance"]
    with open(protocol_path) as fh:
        protocol = json.load(fh)
    if X.sha256_file(protocol_path) != prov["protocol_sha256"]:
        sys.exit("protocol differs from the matrix's")
    if X.sha256_file(extra_readouts) != prov.get("extra_readouts_sha256"):
        sys.exit("extra readouts differ from the matrix's")
    panel = list(dict.fromkeys(rid for p in protocol for rid, _ in p["biomarkers"]))
    with open(extra_readouts) as fh:
        panel = list(dict.fromkeys(panel + [line.strip() for line in fh if line.strip()]))
    by_iem = {r["iem"]: r for r in recs}
    unknown = set(iems) - set(by_iem)
    if unknown:
        sys.exit(f"IEMs not in the matrix: {sorted(unknown)}")
    model_file = os.path.join(ROOT, "external", "COBRA.models", "mat", "Harvey_1_03d.mat")
    m = W.load_wbm(model_file)
    if X.sha256_file(model_file) != prov["model_sha256"]:
        sys.exit("model file differs from the matrix's")
    inputs = WC.load_inputs(WC.DEFAULT_INPUTS)
    lb, ub, _ = WC.runiem_model_setup(m.rxns, m.mets, m.S, m.lb, m.ub, sex=m.meta.get("sex") or "male", inputs=inputs)
    m.lb, m.ub = lb, ub
    hw = X.make_backend(backend, m, threads)
    I.apply_runiem_global_constraints(hw, bile_duct="toolbox")
    bounds = hashlib.sha256(np.ascontiguousarray(hw.lb).tobytes() + np.ascontiguousarray(hw.ub).tobytes()).hexdigest()
    if bounds != prov["lp_bounds_sha256_before_sinks"]:
        sys.exit("LP bounds differ from the matrix's")
    sink_mets = list(dict.fromkeys([rid[3:] for rid in panel if rid.startswith("DM_")] +
                                   [met for p in protocol for met in (p.get("demand_metabolites") or [])]))
    for met in sink_mets:
        if met in hw.met_pos:
            hw.add_demand(met, lb=0.0, ub=0.0)
    chosen = [p for p in protocol if p["iem"] in set(iems)]
    setups = [X.prepare(hw, p, prov["context"], log=False, stored=by_iem[p["iem"]]) for p in chosen]
    bad = [s.iem for s in setups if not s.ok]
    if bad:
        sys.exit(f"IEMs without a usable set-up in the matrix: {bad}")
    return hw, setups, panel, prov


def shard_path(model, set_name, version, shard=None):
    tail = f"_shard{shard[0]}of{shard[1]}" if shard else ""
    return os.path.join(OUT, f"{model}_pin_sweep_{set_name}_{version}{tail}.json")


def child(args, iems, alphas):
    k, n = (int(x) for x in args.shard.split("/"))
    t0 = time.time()
    hw, setups, panel, prov = build(args.matrix, args.extra_readouts, args.protocol, args.backend, args.threads, iems)
    readouts = panel[k - 1::n]
    if args.readout_limit:
        readouts = readouts[:args.readout_limit]
    provenance = {"engine_version": ENGINE_VERSION, "matrix": os.path.relpath(os.path.abspath(args.matrix), ROOT),
                  "matrix_sha256": X.sha256_file(args.matrix), "model_sha256": prov["model_sha256"],
                  "protocol_sha256": prov["protocol_sha256"], "extra_readouts_sha256": prov["extra_readouts_sha256"],
                  "lp_bounds_sha256_before_sinks": prov["lp_bounds_sha256_before_sinks"], "context": prov["context"],
                  "solver": X.backend_version(args.backend), "warm": args.warm, "first": "ipm",
                  "feas_tol": 1e-7, "opt_tol": 1e-7, "tol": X.TOL, "recheck_every": args.recheck_every,
                  "threads": args.threads, "alphas": alphas, "iems": [s.iem for s in setups],
                  "split": args.split, "set": args.set, "shard": f"{k}/{n}",
                  "n_panel": len(panel) if not args.readout_limit else f"test: {args.readout_limit} per shard",
                  "panel_sha256": hashlib.sha256(json.dumps(panel).encode()).hexdigest(),
                  "readouts_sha256": hashlib.sha256(json.dumps(readouts).encode()).hexdigest(), "n_readouts": len(readouts)}
    out = shard_path(args.model, args.set, args.version, (k, n))
    result = {"provenance": provenance,
              "setups": {s.iem: {"vmax_healthy": s.rec["vmax_healthy"], "matrix_pin": s.lo,
                                 "pins": {str(a): pin_for(s.rec["vmax_healthy"], a, s.lo) for a in alphas}} for s in setups},
              "readouts": {}}
    if os.path.exists(out):
        with open(out) as fh:
            old = json.load(fh)
        if old["provenance"] != provenance:
            sys.exit(f"{out} was written with other settings; move it away or use another --version")
        result["readouts"] = old["readouts"]
    print(f"== shard {k}/{n}: {len(readouts)} readouts x {len(setups)} IEMs x {len(alphas)} alphas; "
          f"{len(result['readouts'])} readouts already done; setup {time.time() - t0:.0f}s", flush=True)
    mover = X.Mover(hw)
    for j, rid in enumerate(readouts):
        if rid in result["readouts"]:
            continue
        t_r, n_r = time.time(), hw.n_solves
        result["readouts"][rid] = sweep_readout(hw, mover, setups, rid, alphas, args.warm, "ipm", args.recheck_every,
                                                log=not args.quiet)
        X.write_json(out, result)
        print(f"  readout {j + 1}/{len(readouts)} {rid}: {hw.n_solves - n_r} solves {time.time() - t_r:.0f}s", flush=True)
    mover.reset()
    print("done", flush=True)


def merge(paths, out):
    parts = []
    for path in paths:
        with open(path) as fh:
            parts.append(json.load(fh))
    base = {k: v for k, v in parts[0]["provenance"].items() if k not in ("shard", "readouts_sha256", "n_readouts")}
    for p in parts[1:]:
        if {k: v for k, v in p["provenance"].items() if k not in ("shard", "readouts_sha256", "n_readouts")} != base:
            raise ValueError("shards were run with different settings")
        if p["setups"] != parts[0]["setups"]:
            raise ValueError("shards disagree on the IEM set-ups")
    readouts = {}
    for p in parts:
        n_expected = p["provenance"]["n_readouts"]
        if len(p["readouts"]) != n_expected:
            raise ValueError(f"a shard is incomplete: {len(p['readouts'])} of {n_expected} readouts")
        readouts.update(p["readouts"])
    if isinstance(base["n_panel"], int) and len(readouts) != base["n_panel"]:
        raise ValueError(f"the shards hold {len(readouts)} readouts, not the panel's {base['n_panel']}")
    merged = {"provenance": dict(base, shards=len(parts), merged_from=[os.path.relpath(p, ROOT) for p in paths]),
              "setups": parts[0]["setups"], "readouts": readouts}
    X.write_json(out, merged)
    return merged


def fmt(seconds):
    seconds = int(seconds)
    h, m = divmod(seconds // 60, 60)
    return f"{h}h{m:02d}m" if h else f"{m}m"


def parent(args):
    cores = os.cpu_count() or 4
    existing = {}
    for name in os.listdir(OUT) if os.path.isdir(OUT) else []:
        m = re.match(rf"{re.escape(args.model)}_pin_sweep_{re.escape(args.set)}_{re.escape(args.version)}_shard(\d+)of(\d+)\.json$", name)
        if m:
            existing[(int(m.group(1)), int(m.group(2)))] = os.path.join(OUT, name)
    if existing:
        ns = {n for _, n in existing}
        if len(ns) != 1:
            sys.exit(f"shard files with different N exist: {sorted(existing.values())}")
        n = ns.pop()
        with open(next(iter(existing.values()))) as fh:
            t = json.load(fh)["provenance"]["threads"]
        print(f"resuming {n} shards with {t} threads each", flush=True)
    else:
        n = args.parallel or max(1, min(6, cores // 3))
        t = args.threads or max(1, cores // n)
    os.makedirs(os.path.join(OUT, "logs"), exist_ok=True)
    print(f"{cores} cores: {n} shards x {t} solver threads", flush=True)
    procs = {}
    for k in range(1, n + 1):
        cmd = [sys.executable, os.path.abspath(__file__), "--shard", f"{k}/{n}", "--threads", str(t), "--quiet"]
        for flag in ("matrix", "split", "set", "alphas", "backend", "warm", "extra_readouts", "protocol", "version",
                     "recheck_every", "model", "iems", "readout_limit"):
            if getattr(args, flag):
                cmd += [f"--{flag.replace('_', '-')}", str(getattr(args, flag))]
        log = open(os.path.join(OUT, "logs", f"{args.model}_pin_sweep_{args.set}_{args.version}_shard{k}of{n}.log"), "a")
        log.write(f"\n=== start {time.strftime('%Y-%m-%d %H:%M:%S')}: {' '.join(cmd)}\n")
        log.flush()
        procs[k] = (subprocess.Popen(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT), log)
    paths = {k: shard_path(args.model, args.set, args.version, (k, n)) for k in procs}
    t0, first_seen, start_done = time.time(), None, 0
    try:
        while True:
            alive = [k for k, (p, _) in procs.items() if p.poll() is None]
            failed = [k for k, (p, _) in procs.items() if p.poll() not in (None, 0)]
            if failed:
                for k in failed:
                    procs[k][1].flush()
                    with open(procs[k][1].name) as fh:
                        tail = fh.readlines()[-15:]
                    print(f"\nshard {k} stopped with exit code {procs[k][0].returncode}. Last lines of its log:\n" + "".join(tail),
                          flush=True)
                for k in alive:
                    procs[k][0].send_signal(signal.SIGINT)
                print("stopping the other shards; run the same command again to resume, or tell Claude.", flush=True)
                break
            done, total = 0, 0
            for k in procs:
                try:
                    with open(paths[k]) as fh:
                        d = json.load(fh)
                    done += len(d["readouts"]); total += d["provenance"]["n_readouts"]
                except (OSError, ValueError, KeyError):
                    pass
            if done and first_seen is None:
                first_seen, start_done = time.time(), done
            if not done:
                line = f"setting up ({fmt(time.time() - t0)} so far)"
            else:
                rate = (done - start_done) / max(1.0, time.time() - first_seen)
                eta = f"; about {fmt((total - done) / rate)} to go" if rate > 0 and total > done else ""
                line = f"{done}/{total or '?'} readouts done{eta}"
            print(f"[{time.strftime('%H:%M')}] {line}", flush=True)
            if not alive:
                break
            time.sleep(args.interval)
    except KeyboardInterrupt:
        print("\nstopping the shards (each keeps every readout it finished)...", flush=True)
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
    if any(p.returncode != 0 for p, _ in procs.values()):
        return
    out = shard_path(args.model, args.set, args.version)
    merge([paths[k] for k in sorted(paths)], out)
    print(f"all shards finished and joined: {os.path.relpath(out, ROOT)}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", default=os.path.join("results", "wbm_iem", "Harvey_1_03d_iem_cross_ranking_v1.json"))
    ap.add_argument("--model", default="Harvey_1_03d")
    ap.add_argument("--split", default=os.path.join("data", "iem", "iem_ranking_split_v1.json"))
    ap.add_argument("--set", default="development", choices=["development", "held_out"])
    ap.add_argument("--alphas", default=DEFAULT_ALPHAS)
    ap.add_argument("--backend", choices=["highs", "gurobi"], default="gurobi")
    ap.add_argument("--warm", choices=["primal", "dual", "concurrent", "ipm"], default="concurrent")
    ap.add_argument("--extra-readouts", default=os.path.join("data", "iem", "iem_ranking_extra_readouts_v0.2.txt"))
    ap.add_argument("--protocol", default=os.path.join("data", "iem", "iem_protocol_v0.2.json"))
    ap.add_argument("--version", default="v1")
    ap.add_argument("--recheck-every", type=int, default=100)
    ap.add_argument("--parallel", type=int, help="number of shards (default cores // 3, 1 to 6)")
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--shard", help="K/N: run one shard (used by the parallel launcher)")
    ap.add_argument("--interval", type=int, default=120)
    ap.add_argument("--iems", help="testing only: comma-separated IEMs instead of the split's set")
    ap.add_argument("--readout-limit", type=int, default=0, help="testing only: first N readouts of each shard")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    os.chdir(ROOT)
    alphas = [float(a) for a in args.alphas.split(",")]
    if any(not 0 <= a <= 1 for a in alphas):
        ap.error("every alpha must be in [0, 1]")
    with open(args.split) as fh:
        iems = json.load(fh)[args.set]
    if args.iems:
        iems = args.iems.split(",")
    if args.iems or args.readout_limit:
        # Test runs never write to the real result files.
        args.version = args.version if args.version.startswith("test") else f"test_{args.version}"
    os.makedirs(OUT, exist_ok=True)
    if args.shard:
        child(args, iems, alphas)
    else:
        parent(args)


if __name__ == "__main__":
    main()
