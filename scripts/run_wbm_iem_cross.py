"""Cross-disease biomarker readouts on a whole-body model (disease-ranking study).

Usage:
  python scripts/run_wbm_iem_cross.py Harvey_1_03d [--backend highs|gurobi] [--warm primal|dual|ipm]
        [--context protocol|minimal] [--iems HIS PKU ...] [--limit N] [--out-suffix _cross_test]
        [--check-against results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json]

Every IEM is set up as gembench.wbm_iem._run_iem does (the runIEM_HH / checkIEM_WBM protocol): the IEM's
bound tweaks; the healthy reference v_max (maximum summed IEM flux, auxiliary bound +/-1e5); the healthy
pin sum(IEM) >= v_max truncated to 6 decimals; the disease state with lb = ub = 0 on the IEM reactions; and
the whole-body objective check in the disease state. The protocol then maximises only the IEM's own
biomarkers. Here every reaction in the readout set (default: the union of all biomarkers in the protocol,
in protocol order) is maximised, first in the healthy state and then in the disease state, with its upper
bound raised to 1e5 for its own solve, as the protocol does for a biomarker.

--context protocol  the IEM's own demand sinks (its demand_metabolites and DM_ biomarkers) are open at
                    ub 1000 throughout, as in the protocol. The own-biomarker readouts are then the same LPs
                    as the protocol's, so they reproduce a protocol run.
--context minimal   no demand sink is open except the one being maximised, for every IEM alike.
In both contexts a demand sink that is not open is closed (ub 0) except while it is being maximised.

--warm primal       the first readout in each state is solved by interior point with crossover. The others
                    start from the previous basis with primal simplex: only the objective, and the upper
                    bound of the readout being maximised, change between them. A warm solve that does not
                    end optimal is solved again by interior point, and this is recorded.
--warm ipm          interior point with crossover for every solve, as the protocol runs.

--order iem         (default) one IEM at a time: all readouts in its healthy state, then in its disease state.
--order readout     one readout at a time: every IEM's disease state, then every IEM's healthy state. The
                    objective stays the same between consecutive solves and only bounds change, so the warm
                    start (--warm dual) begins from a dual-feasible basis. The per-IEM protocol solves (v_max,
                    disease check, whole-body check) run first for all IEMs. Each LP is the same as in --order iem.
                    The first solve of each readout uses interior point (--first ipm) unless --first warm.

--extra-readouts F  adds the readouts listed in F (one reaction per line) to the panel, after the protocol's.
--recheck-every N   re-solves about one in N warm-started readouts from scratch by interior point (chosen by a
                    hash of IEM, state and reaction) and records the difference.

The decision rule is the protocol's: values with |f| <= 1e-6 are 0, and disease - healthy > 1e-6 is
Increased, < -1e-6 Decreased, otherwise Unchanged. Failed solves are NA, never zero.
Results are fingerprinted and written to results/wbm_iem/<model>_iem_cross<suffix>.json; --order iem resumes per IEM,
--order readout per readout.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gembench import wbm as W  # noqa: E402
from gembench import wbm_iem as I  # noqa: E402
from gembench import wbm_constraints as WC  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "wbm_iem")
ENGINE_VERSION = "iem-cross-v0.1"
TOL = 1e-6
TERMINAL_STATUSES = {"complete", "inactive", "no_reactions", "disease_infeasible", "wb_infeasible", "missing_wb_objective"}


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.floating,)):
        return json_safe(float(value))
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def write_json(path, value):
    temporary = path + ".tmp"
    with open(temporary, "w") as fh:
        json.dump(json_safe(value), fh, indent=1, allow_nan=False)
    os.replace(temporary, path)


def call(fh, fd):
    if not (math.isfinite(fh) and math.isfinite(fd)):
        return "NA"
    diff = fd - fh
    return "Increased" if diff > TOL else "Decreased" if diff < -TOL else "Unchanged"


def make_backend(name, m, threads):
    if name == "highs":
        return I.HighsWBM(m, threads=threads)
    if name == "gurobi":
        from gembench import wbm_iem_gurobi as G
        return G.GurobiWBM(m, threads=threads)
    raise ValueError(name)


def backend_version(name):
    if name == "highs":
        import highspy
        return f"HiGHS {highspy.Highs().version()}"
    import gurobipy
    return "Gurobi " + ".".join(str(v) for v in gurobipy.gurobi.version())


def solve_logged(hw, method, label, log, show_value=True):
    hw.method = method
    st, f, _, dt = hw.solve()
    info = dict(getattr(hw, "last_info", {}) or {})
    fallback = False
    if method != "ipm" and st != "Optimal":
        # A warm start that does not end optimal is solved again from scratch; both are recorded.
        hw.method = "ipm"
        st2, f2, _, dt2 = hw.solve()
        info = {"first_attempt": {"method": method, "status": st, "time_s": dt, **info},
                **dict(getattr(hw, "last_info", {}) or {})}
        st, f, dt, fallback = st2, f2, dt + dt2, True
    if log:
        it = info.get("simplex_iterations", "")
        shown = f"f={f:14.8g}" if show_value else " " * 16
        print(f"      {label:44s} {st:10s} {shown} {dt:7.2f}s it={it}{' (fallback ipm)' if fallback else ''}", flush=True)
    return st, f, dt, info, fallback


def rechecked(iem, state, rid, every):
    if not every:
        return False
    return int(hashlib.md5(f"{iem}|{state}|{rid}".encode()).hexdigest(), 16) % every == 0


def run_one(hw, p, readouts, context, warm, log=True, recheck_every=0):
    t0 = time.time(); n0 = hw.n_solves
    rx = hw.wbm.rxns
    rec = {"iem": p["iem"], "call_index": p["call_index"], "context": context, "status": "not_run", "notes": [],
           "vmax_healthy": float("nan"), "vmax_disease": float("nan"), "wb_objective_disease_feasible": False,
           "iem_reactions": [], "readouts": []}
    own = {rid: I._expected(label) for rid, label in p["biomarkers"]}
    with hw.temporary_state():
        for tw in p["bound_tweaks"]:
            idx = I.match_reactions(rx, [tw["pattern"]])
            hw.set_bounds(idx, **{tw["bound"]: [tw["value"]] * len(idx)})
        if context == "protocol":
            required = list(dict.fromkeys(list(p.get("demand_metabolites") or []) +
                                          [rid[3:] for rid, _ in p["biomarkers"] if rid.startswith("DM_")]))
            for met in required:
                rid = f"DM_{met}"
                if rid in hw.rxn_pos:
                    hw.set_bounds([hw.rxn_pos[rid]], lb=[0.0], ub=[1000.0])
                else:
                    rec["notes"].append(f"metabolite for {rid} not in model")
        iem_idx = I.match_reactions(rx, p["include_patterns"], p["exclude_patterns"])
        rec["iem_reactions"] = [rx[i] for i in iem_idx]
        if not iem_idx:
            rec["status"] = "no_reactions"
            return finish(rec, hw, n0, t0)
        saved_lb = hw.lb[iem_idx].copy(); saved_ub = hw.ub[iem_idx].copy()
        row = hw.add_row({i: 1.0 for i in iem_idx}, -1e5, 1e5)
        hw.set_objective({i: 1.0 for i in iem_idx}, "max")
        st, vmax, _, _, _ = solve_logged(hw, "ipm", "v_max (healthy)", log)
        rec["vmax_healthy"] = vmax
        if st != "Optimal" or not np.isfinite(vmax) or abs(vmax) <= TOL:
            rec["status"] = "inactive" if st == "Optimal" and np.isfinite(vmax) and abs(vmax) <= TOL else "healthy_solve_failed"
            return finish(rec, hw, n0, t0)
        lo = math.floor(vmax * 1e6) / 1e6 if vmax > 0 else math.ceil(vmax * 1e6) / 1e6
        rec["healthy_pin"] = lo

        def healthy_state():
            hw.set_bounds(iem_idx, lb=saved_lb, ub=saved_ub)
            hw.set_row_bounds(row, lo, 1e5)

        def disease_state():
            hw.set_bounds(iem_idx, lb=[0.0] * len(iem_idx), ub=[0.0] * len(iem_idx))
            hw.set_row_bounds(row, 0.0, 1e5)

        disease_state()
        st, vd, _, _, _ = solve_logged(hw, "ipm", "summed IEM flux (disease)", log)
        rec["vmax_disease"] = vd if st == "Optimal" and np.isfinite(vd) else float("nan")
        if st != "Optimal" or not np.isfinite(vd):
            rec["status"] = "disease_infeasible" if st == "Infeasible" else "disease_solve_failed"
            return finish(rec, hw, n0, t0)
        wb = hw.rxn_pos.get("Whole_body_objective_rxn")
        if wb is None:
            rec["status"] = "missing_wb_objective"
            return finish(rec, hw, n0, t0)
        hw.set_objective({wb: 1.0}, "max")
        st, fw, _, _, _ = solve_logged(hw, "ipm", "whole-body objective (disease)", log)
        rec["wb_objective_disease_feasible"] = bool(st == "Optimal" and np.isfinite(fw))
        if not rec["wb_objective_disease_feasible"]:
            rec["status"] = "wb_infeasible" if st == "Infeasible" else "wb_solve_failed"
            healthy_state(); hw.set_row_bounds(row, -np.inf, np.inf)
            return finish(rec, hw, n0, t0)

        values = {}
        for state, setter in (("healthy", healthy_state), ("disease", disease_state)):
            setter()
            first = True
            for k, rid in enumerate(readouts):
                if rid not in hw.rxn_pos:
                    values[(state, rid)] = ("absent", float("nan"), 0.0, {}, False)
                    continue
                j = hw.rxn_pos[rid]
                old_ub = hw.ub[j]
                hw.set_bounds([j], ub=[1e5])
                hw.set_objective({j: 1.0}, "max")
                method = "ipm" if first or warm == "ipm" else warm
                # Values are printed for the IEM's own biomarkers only: cross-disease values are not looked at
                # while the ranking plan is still open.
                st, f, dt, info, fb = solve_logged(hw, method, f"{p['iem']} {state} {k + 1}/{len(readouts)} {rid}", log,
                                                   show_value=rid in own)
                first = False
                if method != "ipm" and rechecked(p["iem"], state, rid, recheck_every):
                    st2, f2, dt2, _, _ = solve_logged(hw, "ipm", f"{p['iem']} {state} recheck {rid}", log,
                                                      show_value=rid in own)
                    info = dict(info, recheck={"status": st2, "value": f2, "time_s": round(dt2, 3),
                                               "abs_diff": abs(f2 - f) if np.isfinite(f) and np.isfinite(f2) else None})
                hw.set_bounds([j], ub=[old_ub])
                ok = st == "Optimal" and np.isfinite(f)
                f = (0.0 if abs(f) <= TOL else f) if ok else float("nan")
                values[(state, rid)] = (st, f, dt, info, fb)
        healthy_state(); hw.set_row_bounds(row, -np.inf, np.inf)

    unavailable = 0
    for rid in readouts:
        sth, fh, th, ih, fbh = values[("healthy", rid)]
        std, fd, td, idd, fbd = values[("disease", rid)]
        pred = call(fh, fd) if sth == std == "Optimal" else "NA"
        if sth != "absent" and pred == "NA":
            unavailable += 1
        rec["readouts"].append({"reaction": rid, "own": rid in own, "expected": own.get(rid),
                                "healthy": fh, "disease": fd, "predicted": pred,
                                "status_healthy": sth, "status_disease": std,
                                "time_healthy_s": round(th, 3), "time_disease_s": round(td, 3),
                                "solve_healthy": ih, "solve_disease": idd,
                                "fallback_healthy": fbh, "fallback_disease": fbd})
    rec["status"] = "partial" if unavailable else "complete"
    return finish(rec, hw, n0, t0)


def finish(rec, hw, n0, t0):
    rec["n_solves"] = hw.n_solves - n0
    rec["time_s"] = round(time.time() - t0, 1)
    return rec


class Setup:
    """One IEM's protocol set-up, kept so that its healthy and disease LPs can be rebuilt as bound overrides."""

    def __init__(self, p, context):
        self.p, self.iem, self.call_index, self.context = p, p["iem"], p["call_index"], context
        self.own = {rid: I._expected(label) for rid, label in p["biomarkers"]}
        self.ok = False
        self.rec = {"iem": p["iem"], "call_index": p["call_index"], "context": context, "status": "not_run", "notes": [],
                    "vmax_healthy": float("nan"), "vmax_disease": float("nan"), "wb_objective_disease_feasible": False,
                    "iem_reactions": [], "readouts": []}


def prepare(hw, p, context, log=True):
    """Protocol steps 1-3 for one IEM (as in run_one); records the overrides that define its two states."""
    t0, n0 = time.time(), hw.n_solves
    s = Setup(p, context)
    rx = hw.wbm.rxns
    with hw.temporary_state():
        tweaks = {}
        for tw in p["bound_tweaks"]:
            idx = I.match_reactions(rx, [tw["pattern"]])
            hw.set_bounds(idx, **{tw["bound"]: [tw["value"]] * len(idx)})
            for col in idx:
                tweaks[col] = (float(hw.lb[col]), float(hw.ub[col]))
        sinks = []
        if context == "protocol":
            required = list(dict.fromkeys(list(p.get("demand_metabolites") or []) +
                                          [rid[3:] for rid, _ in p["biomarkers"] if rid.startswith("DM_")]))
            for met in required:
                rid = f"DM_{met}"
                if rid in hw.rxn_pos:
                    hw.set_bounds([hw.rxn_pos[rid]], lb=[0.0], ub=[1000.0])
                    sinks.append(hw.rxn_pos[rid])
                else:
                    s.rec["notes"].append(f"metabolite for {rid} not in model")
        s.tweaks, s.sinks = tweaks, sinks
        s.iem_idx = I.match_reactions(rx, p["include_patterns"], p["exclude_patterns"])
        s.rec["iem_reactions"] = [rx[i] for i in s.iem_idx]
        if not s.iem_idx:
            s.rec["status"] = "no_reactions"
            return finish_setup(s, hw, n0, t0)
        s.saved = [(float(hw.lb[i]), float(hw.ub[i])) for i in s.iem_idx]
        s.row = hw.add_row({i: 1.0 for i in s.iem_idx}, -1e5, 1e5)
        hw.set_objective({i: 1.0 for i in s.iem_idx}, "max")
        st, vmax, _, _, _ = solve_logged(hw, "ipm", f"{s.iem} v_max (healthy)", log)
        s.rec["vmax_healthy"] = vmax
        if st != "Optimal" or not np.isfinite(vmax) or abs(vmax) <= TOL:
            s.rec["status"] = "inactive" if st == "Optimal" and np.isfinite(vmax) and abs(vmax) <= TOL else "healthy_solve_failed"
            return finish_setup(s, hw, n0, t0)
        s.lo = math.floor(vmax * 1e6) / 1e6 if vmax > 0 else math.ceil(vmax * 1e6) / 1e6
        s.rec["healthy_pin"] = s.lo
        hw.set_bounds(s.iem_idx, lb=[0.0] * len(s.iem_idx), ub=[0.0] * len(s.iem_idx))
        hw.set_row_bounds(s.row, 0.0, 1e5)
        st, vd, _, _, _ = solve_logged(hw, "ipm", f"{s.iem} summed IEM flux (disease)", log)
        s.rec["vmax_disease"] = vd if st == "Optimal" and np.isfinite(vd) else float("nan")
        if st != "Optimal" or not np.isfinite(vd):
            s.rec["status"] = "disease_infeasible" if st == "Infeasible" else "disease_solve_failed"
            return finish_setup(s, hw, n0, t0)
        wb = hw.rxn_pos.get("Whole_body_objective_rxn")
        if wb is None:
            s.rec["status"] = "missing_wb_objective"
            return finish_setup(s, hw, n0, t0)
        hw.set_objective({wb: 1.0}, "max")
        st, fw, _, _, _ = solve_logged(hw, "ipm", f"{s.iem} whole-body objective (disease)", log)
        s.rec["wb_objective_disease_feasible"] = bool(st == "Optimal" and np.isfinite(fw))
        if not s.rec["wb_objective_disease_feasible"]:
            s.rec["status"] = "wb_infeasible" if st == "Infeasible" else "wb_solve_failed"
            return finish_setup(s, hw, n0, t0)
        s.ok = True
        s.rec["status"] = "set_up"
    return finish_setup(s, hw, n0, t0)


def finish_setup(s, hw, n0, t0):
    s.rec["n_setup_solves"] = hw.n_solves - n0
    s.rec["setup_time_s"] = round(time.time() - t0, 1)
    return s


def target_bounds(s, state, rcol, base_lb, base_ub):
    """Column overrides of the base LP for IEM s, state and readout column rcol, in the protocol's order:
    bound tweaks, open own sinks, readout upper bound 1e5, then the state's bounds on the IEM reactions."""
    t = dict(s.tweaks)
    for col in s.sinks:
        t[col] = (0.0, 1000.0)
    lb_r = t.get(rcol, (float(base_lb[rcol]), float(base_ub[rcol])))[0]
    t[rcol] = (lb_r, 1e5)
    for k, col in enumerate(s.iem_idx):
        t[col] = s.saved[k] if state == "healthy" else (0.0, 0.0)
    return t, (s.lo, 1e5) if state == "healthy" else (0.0, 1e5)


class Mover:
    """Applies a target set of overrides to the persistent LP, changing only what differs from the current one."""

    def __init__(self, hw):
        self.hw = hw
        self.base_lb, self.base_ub = hw.lb.copy(), hw.ub.copy()
        self.current, self.row = {}, None

    def base(self, col):
        return float(self.base_lb[col]), float(self.base_ub[col])

    def move(self, target, row=None, row_bounds=None):
        cols, lbs, ubs = [], [], []
        for col in set(self.current) | set(target):
            want = target.get(col, self.base(col))
            have = self.current.get(col, self.base(col))
            if want != have:
                cols.append(col); lbs.append(want[0]); ubs.append(want[1])
        if cols:
            self.hw.set_bounds(cols, lb=lbs, ub=ubs)
        self.current = {c: v for c, v in target.items() if v != self.base(c)}
        if self.row is not None and self.row != row:
            self.hw.set_row_bounds(self.row, -np.inf, np.inf)
        if row is not None:
            self.hw.set_row_bounds(row, *row_bounds)
        self.row = row

    def reset(self):
        self.move({}, None)


def readout_entry(rid, own, values):
    sth, fh, th, ih, fbh = values["healthy"]
    std, fd, td, idd, fbd = values["disease"]
    pred = call(fh, fd) if sth == std == "Optimal" else "NA"
    return {"reaction": rid, "own": rid in own, "expected": own.get(rid), "healthy": fh, "disease": fd, "predicted": pred,
            "status_healthy": sth, "status_disease": std, "time_healthy_s": round(th, 3), "time_disease_s": round(td, 3),
            "solve_healthy": ih, "solve_disease": idd, "fallback_healthy": fbh, "fallback_disease": fbd}


def run_readout_major(hw, protocol, readouts, context, warm, first, log, recheck_every, done, on_readout):
    """All IEMs' protocol set-ups, then one readout at a time across every IEM's disease and healthy states.

    done: {(iem, readout): entry} already computed (resume); on_readout(setups, entries) is called after each readout.
    """
    setups = [prepare(hw, p, context, log) for p in protocol]
    overlap = {readouts[k] for k in range(len(readouts)) if readouts[k] in hw.rxn_pos
               for s in setups if s.ok and hw.rxn_pos[readouts[k]] in set(s.iem_idx)}
    if overlap:
        raise ValueError(f"readouts that are also IEM reactions are not supported: {sorted(overlap)}")
    mover = Mover(hw)
    entries = dict(done)
    ok = [s for s in setups if s.ok]
    for k, rid in enumerate(readouts):
        if all((s.iem, rid) in entries for s in ok):
            continue
        t_r, n_r = time.time(), hw.n_solves
        if rid not in hw.rxn_pos:
            for s in ok:
                entries[(s.iem, rid)] = readout_entry(rid, s.own, {"healthy": ("absent", float("nan"), 0.0, {}, False),
                                                                   "disease": ("absent", float("nan"), 0.0, {}, False)})
            on_readout(setups, entries)
            continue
        rcol = hw.rxn_pos[rid]
        hw.set_objective({rcol: 1.0}, "max")
        values = {}
        is_first = True
        for state in ("disease", "healthy"):
            for s in ok:
                target, row_bounds = target_bounds(s, state, rcol, mover.base_lb, mover.base_ub)
                mover.move(target, s.row, row_bounds)
                method = "ipm" if (is_first and first == "ipm") or warm == "ipm" else warm
                st, f, dt, info, fb = solve_logged(hw, method, f"{rid} {k + 1}/{len(readouts)} {s.iem} {state}", log,
                                                   show_value=rid in s.own)
                is_first = False
                if method != "ipm" and rechecked(s.iem, state, rid, recheck_every):
                    st2, f2, dt2, _, _ = solve_logged(hw, "ipm", f"{rid} {s.iem} {state} recheck", log, show_value=rid in s.own)
                    info = dict(info, recheck={"status": st2, "value": f2, "time_s": round(dt2, 3),
                                               "abs_diff": abs(f2 - f) if np.isfinite(f) and np.isfinite(f2) else None})
                okv = st == "Optimal" and np.isfinite(f)
                f = (0.0 if abs(f) <= TOL else f) if okv else float("nan")
                values[(s.iem, state)] = (st, f, dt, info, fb)
        for s in ok:
            entries[(s.iem, rid)] = readout_entry(rid, s.own, {"healthy": values[(s.iem, "healthy")],
                                                               "disease": values[(s.iem, "disease")]})
        times = [v[2] for v in values.values()]
        warm_times = sorted(v[2] for v in values.values() if (v[3] or {}).get("method") not in (None, "ipm"))
        print(f"  readout {k + 1}/{len(readouts)} {rid}: {hw.n_solves - n_r} solves {time.time() - t_r:.0f}s; "
              f"warm median {np.median(warm_times) if warm_times else float('nan'):.2f}s max {max(times):.1f}s; "
              f"fallbacks {sum(v[4] for v in values.values())}", flush=True)
        on_readout(setups, entries)
    mover.reset()
    return setups, entries


def assemble(setups, entries, readouts, final):
    """Per-IEM records in the same shape as --order iem; status complete/partial once every readout is in."""
    out = []
    for s in setups:
        rec = dict(s.rec)
        rec["readouts"] = [entries[(s.iem, r)] for r in readouts if (s.iem, r) in entries]
        if s.ok:
            all_in = len(rec["readouts"]) == len(readouts)
            unavailable = sum(1 for e in rec["readouts"] if e["status_healthy"] != "absent" and e["predicted"] == "NA")
            rec["status"] = ("partial" if unavailable else "complete") if all_in else "readouts_incomplete"
        out.append(rec)
    return out


def compare_with_reference(rec, ref_by_iem):
    """Own-biomarker readouts against a protocol run of the same IEM (e.g. v0.3)."""
    ref = ref_by_iem.get((rec["iem"], rec["call_index"]))
    if ref is None:
        return None
    by_rxn = {r["reaction"]: r for r in rec["readouts"]}
    rows = []
    for b in ref["biomarkers"]:
        r = by_rxn.get(b["reaction"])
        if r is None:
            continue
        def d(a, c):
            if a is None or c is None or not (math.isfinite(a) and math.isfinite(c)):
                return None
            return abs(a - c)
        rows.append({"reaction": b["reaction"], "ref_healthy": b["healthy"], "ref_disease": b["disease"],
                     "healthy": r["healthy"], "disease": r["disease"],
                     "abs_diff_healthy": d(b["healthy"], r["healthy"]), "abs_diff_disease": d(b["disease"], r["disease"]),
                     "ref_call": b["predicted"], "call": r["predicted"], "same_call": b["predicted"] == r["predicted"]})
    vmax_diff = abs(rec["vmax_healthy"] - ref["vmax_healthy"]) if ref.get("vmax_healthy") is not None else None
    return {"n": len(rows), "n_same_call": sum(x["same_call"] for x in rows), "vmax_abs_diff": vmax_diff,
            "max_abs_diff": max([x["abs_diff_healthy"] or 0 for x in rows] + [x["abs_diff_disease"] or 0 for x in rows] + [0]),
            "rows": rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model")
    ap.add_argument("--model-file")
    ap.add_argument("--protocol", default=os.path.join(ROOT, "data", "iem", "iem_protocol_v0.2.json"))
    ap.add_argument("--backend", choices=["highs", "gurobi"], default="highs")
    ap.add_argument("--warm", choices=["primal", "dual", "ipm"], default="primal")
    ap.add_argument("--context", choices=["protocol", "minimal"], default="protocol")
    ap.add_argument("--iems", nargs="*")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--out-suffix", default="_cross_test")
    ap.add_argument("--check-against", help="protocol results file to compare own-biomarker readouts with")
    ap.add_argument("--extra-readouts", help="file with one additional readout reaction per line")
    ap.add_argument("--recheck-every", type=int, default=0, help="re-solve about 1 in N warm solves from scratch")
    ap.add_argument("--order", choices=["iem", "readout"], default="iem")
    ap.add_argument("--first", choices=["ipm", "warm"], default="ipm", help="--order readout: first solve of each readout")
    ap.add_argument("--readout-limit", type=int, default=0, help="only the first N readouts of the panel (timing tests)")
    ap.add_argument("--readout-stride", type=int, default=0, help="only every K-th readout of the panel (timing tests)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    with open(args.protocol) as fh:
        full_protocol = json.load(fh)
    readouts = list(dict.fromkeys(rid for p in full_protocol for rid, _ in p["biomarkers"]))
    if args.extra_readouts:
        with open(args.extra_readouts) as fh:
            readouts = list(dict.fromkeys(readouts + [line.strip() for line in fh if line.strip()]))
    if args.readout_stride:
        readouts = readouts[::args.readout_stride]
    if args.readout_limit:
        readouts = readouts[:args.readout_limit]
    protocol = full_protocol
    if args.iems:
        unknown = set(args.iems) - {p["iem"] for p in protocol}
        if unknown:
            ap.error(f"IEMs not in this protocol: {sorted(unknown)}")
        order = {name: k for k, name in enumerate(args.iems)}
        protocol = sorted([p for p in protocol if p["iem"] in order], key=lambda p: order[p["iem"]])
    if args.limit:
        protocol = protocol[:args.limit]

    t0 = time.time()
    model_file = args.model_file or os.path.join(ROOT, "external", "COBRA.models", "mat", f"{args.model}.mat")
    m = W.load_wbm(model_file)
    inputs = WC.load_inputs(WC.DEFAULT_INPUTS)
    lb, ub, report = WC.runiem_model_setup(m.rxns, m.mets, m.S, m.lb, m.ub, sex=m.meta.get("sex") or "male", inputs=inputs)
    m.lb, m.ub = lb, ub
    hw = make_backend(args.backend, m, args.threads)
    g = I.apply_runiem_global_constraints(hw, bile_duct="toolbox")
    bounds_protocol = hashlib.sha256(np.ascontiguousarray(hw.lb).tobytes() + np.ascontiguousarray(hw.ub).tobytes()).hexdigest()
    # Every demand sink any IEM or readout can need, added once and closed, so column positions do not
    # depend on which IEMs run in this process.
    sink_mets = list(dict.fromkeys([rid[3:] for rid in readouts if rid.startswith("DM_")] +
                                   [met for p in full_protocol for met in (p.get("demand_metabolites") or [])]))
    missing_mets = []
    for met in sink_mets:
        if met in hw.met_pos:
            hw.add_demand(met, lb=0.0, ub=0.0)
        else:
            missing_mets.append(met)
    provenance = {"engine_version": ENGINE_VERSION, "model_sha256": sha256_file(model_file),
                  "protocol_sha256": sha256_file(args.protocol), "constraint_inputs_sha256": sha256_file(WC.DEFAULT_INPUTS),
                  "toolbox_commit": inputs["provenance"]["toolbox_commit"], "model_setup": "toolbox", "bile_duct": "toolbox",
                  "global_constraints": g, "lp_bounds_sha256_before_sinks": bounds_protocol,
                  "solver": backend_version(args.backend), "warm": args.warm, "context": args.context,
                  "feas_tol": 1e-7, "opt_tol": 1e-7, "time_limit_per_solve_s": 1800, "threads": args.threads,
                  "readouts_sha256": hashlib.sha256(json.dumps(readouts).encode()).hexdigest(), "n_readouts": len(readouts),
                  "sink_metabolites_absent_from_model": missing_mets, "tol": TOL,
                  "extra_readouts_sha256": sha256_file(args.extra_readouts) if args.extra_readouts else None,
                  "recheck_every": args.recheck_every}
    if args.order != "iem":
        # Recorded only when it differs from the first runs' default, so their fingerprints stay valid.
        provenance.update(order=args.order, first=args.first)
    fingerprint = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
    out_json = os.path.join(OUT, f"{args.model}_iem_cross{args.out_suffix}.json")
    results = []
    if os.path.exists(out_json):
        with open(out_json) as fh:
            results = json.load(fh)
        if any(r.get("run_fingerprint") != fingerprint for r in results):
            sys.exit("Existing results use different settings; choose a new --out-suffix.")
    done = {(r["iem"], r["call_index"]) for r in results if r.get("status") in TERMINAL_STATUSES}
    todo = [p for p in protocol if (p["iem"], p["call_index"]) not in done]
    ref_by_iem = {}
    if args.check_against:
        with open(args.check_against) as fh:
            ref_by_iem = {(r["iem"], r["call_index"]): r for r in json.load(fh)}
    print(f"== {args.model} {provenance['solver']} warm={args.warm} context={args.context} order={args.order}: "
          f"{len(readouts)} readouts; bounds {bounds_protocol[:12]}; "
          f"{len(protocol) if args.order == 'readout' else len(todo)} IEMs to run; setup {time.time() - t0:.0f}s", flush=True)

    if args.order == "readout":
        stored = {r["iem"]: r for r in results}
        done = {(r["iem"], e["reaction"]): e for r in results for e in r.get("readouts", [])}

        def checkpoint(setups, entries):
            for s in setups:
                old = stored.get(s.iem)
                if old is not None and s.ok and old.get("vmax_healthy") is not None and \
                        abs(old["vmax_healthy"] - s.rec["vmax_healthy"]) > 1e-6 * max(1.0, abs(s.rec["vmax_healthy"])):
                    sys.exit(f"{s.iem}: v_max {s.rec['vmax_healthy']} differs from the stored {old['vmax_healthy']}; not resuming")
            recs = assemble(setups, entries, readouts, final=False)
            for rec in recs:
                rec.update(run_fingerprint=fingerprint, provenance=provenance)
                if rec["status"] in ("complete", "partial") and ref_by_iem:
                    cmp = compare_with_reference(rec, ref_by_iem)
                    if cmp is not None:
                        rec["check_against_reference"] = {"file": os.path.relpath(args.check_against, ROOT), **cmp}
            write_json(out_json, sorted(recs, key=lambda item: item["call_index"]))
            return recs

        setups, entries = run_readout_major(hw, protocol, readouts, args.context, args.warm, args.first,
                                            not args.quiet, args.recheck_every, done, checkpoint)
        recs = checkpoint(setups, entries)
        cmps = [r["check_against_reference"] for r in recs if "check_against_reference" in r]
        rechecks = [e[f"solve_{st}"]["recheck"] for r in recs for e in r["readouts"] for st in ("healthy", "disease")
                    if "recheck" in (e[f"solve_{st}"] or {})]
        diffs = [c["abs_diff"] for c in rechecks if c["abs_diff"] is not None]
        print(f"statuses {sorted({r['status'] for r in recs})}; own biomarkers vs reference: "
              f"{sum(c['n_same_call'] for c in cmps)}/{sum(c['n'] for c in cmps)} same call, max |diff| "
              f"{max([c['max_abs_diff'] for c in cmps] + [0]):.3g}; rechecks {len(rechecks)} "
              f"(max |diff| {max(diffs + [0]):.2g}, non-optimal {sum(c['status'] != 'Optimal' for c in rechecks)})", flush=True)
        print(f"done in {time.time() - t0:.0f}s", flush=True)
        return

    for p in todo:
        rec = run_one(hw, p, readouts, args.context, args.warm, log=not args.quiet, recheck_every=args.recheck_every)
        rec.update(run_fingerprint=fingerprint, provenance=provenance)
        cmp = compare_with_reference(rec, ref_by_iem) if ref_by_iem else None
        if cmp is not None:
            rec["check_against_reference"] = {"file": os.path.relpath(args.check_against, ROOT), **cmp}
        times = [r["time_healthy_s"] for r in rec["readouts"]] + [r["time_disease_s"] for r in rec["readouts"]]
        warm_times = [r[f"time_{s}_s"] for r in rec["readouts"] for s in ("healthy", "disease")
                      if (r[f"solve_{s}"] or {}).get("method") not in (None, "ipm")]
        rechecks = [r[f"solve_{s}"]["recheck"] for r in rec["readouts"] for s in ("healthy", "disease")
                    if "recheck" in (r[f"solve_{s}"] or {})]
        diffs = [c["abs_diff"] for c in rechecks if c["abs_diff"] is not None]
        recheck_text = (f"; rechecks {len(rechecks)} (max |diff| {max(diffs):.2g}, non-optimal "
                        f"{sum(c['status'] != 'Optimal' for c in rechecks)})" if rechecks else "")
        print(f"  {p['iem']:8s} status={rec['status']} solves={rec['n_solves']} {rec['time_s']:.0f}s{recheck_text}; "
              f"readout solve time total {sum(times):.0f}s, warm median {np.median(warm_times) if warm_times else float('nan'):.2f}s "
              f"(n={len(warm_times)})" + (f"; own biomarkers vs reference: {cmp['n_same_call']}/{cmp['n']} same call, "
                                         f"max |diff| {cmp['max_abs_diff']:.3g}, vmax |diff| {cmp['vmax_abs_diff']:.3g}" if cmp else ""),
              flush=True)
        results = [old for old in results if (old["iem"], old["call_index"]) != (p["iem"], p["call_index"])]
        results.append(rec)
        results.sort(key=lambda item: item["call_index"])
        write_json(out_json, results)
    print(f"done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
