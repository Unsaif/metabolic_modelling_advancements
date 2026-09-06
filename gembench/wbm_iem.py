"""Inborn-error-of-metabolism biomarker protocol on whole-body models, in Python.

Port of the COBRA Toolbox's `checkIEM_WBM` (Thiele 2018–2021) and the decision rule at the end of
`runIEM_HH.m`, using a persistent HiGHS model so that the dozens of LPs per IEM warm-start from the
previous model structure. The default is interior point with crossover for each solve.

Protocol (checkIEM_WBM defaults as used in runIEM_HH: minRxnsFluxHealthy = 1, complete knockout):
  1. Healthy model: maximise the summed flux through the IEM reactions, bounded to +/-1e5
     as in the Toolbox's auxiliary variable -> v_max.
  2. If |v_max| > 1e-6: healthy model gets the constraint  sum(IEM fluxes) >= minRxnsFluxHealthy * v_max
     (truncated to 6 decimals); disease model gets lb = ub = 0 on every IEM reaction.
  3. Check the disease model can still satisfy the whole-body objective (fixed at 1).
  4. For each biomarker reaction: raise its upper bound to 1e5, maximise it in the healthy and in the
     disease model; values with |f| <= 1e-6 are set to 0.
  5. Direction call: disease - healthy > 1e-6 -> increased; < -1e-6 -> decreased; else unchanged.
     Failed/nonfinite solves are unknown, never zero-flux predictions.
"""
from __future__ import annotations

import math
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .wbm import WBM, stacked_constraints


class HighsWBM:
    """A whole-body LP held in HiGHS with incremental bound/objective changes and warm starts."""

    def __init__(self, wbm: WBM, feas_tol: float = 1e-7, opt_tol: float = 1e-7, threads: int = 0, time_limit: float = 1800.0):
        import highspy
        self.hs = highspy
        self.wbm = wbm
        self.n = wbm.n_rxns
        A, rl, ru = stacked_constraints(wbm)
        self.m0 = A.shape[0]
        lp = highspy.HighsLp()
        lp.num_col_ = self.n; lp.num_row_ = self.m0
        lp.col_cost_ = np.zeros(self.n)
        lp.col_lower_ = wbm.lb.copy(); lp.col_upper_ = wbm.ub.copy()
        lp.row_lower_ = rl; lp.row_upper_ = ru
        lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
        lp.a_matrix_.start_ = A.indptr.astype(np.int32)
        lp.a_matrix_.index_ = A.indices.astype(np.int32)
        lp.a_matrix_.value_ = A.data.astype(float)
        lp.sense_ = highspy.ObjSense.kMaximize
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        h.setOptionValue("solver", "simplex")
        h.setOptionValue("primal_feasibility_tolerance", feas_tol)
        h.setOptionValue("dual_feasibility_tolerance", opt_tol)
        h.setOptionValue("threads", threads)
        h.setOptionValue("time_limit", time_limit)
        h.passModel(lp)
        self.h = h
        # HiGHS shares one process-wide scheduler. Zero cooperates with an
        # existing scheduler; an explicit incompatible count is a run error.
        self.threads = threads
        self.per_solve_time_limit = time_limit
        self.rxn_pos = {r: i for i, r in enumerate(wbm.rxns)}
        self.met_pos = {m: i for i, m in enumerate(wbm.mets)}
        self.lb = wbm.lb.copy(); self.ub = wbm.ub.copy()
        self.n_extra_rows = 0
        self.n_solves = 0
        self.solve_time = 0.0
        self.method = "ipm"

    # --- structure -----------------------------------------------------------------
    def add_demand(self, met_id: str, lb: float = 0.0, ub: float = 1000.0) -> str:
        """Add a demand reaction DM_<met> (met -> ∅) if absent; returns its id."""
        rid = f"DM_{met_id}"
        if rid in self.rxn_pos:
            return rid
        mi = self.met_pos[met_id]
        self.h.addCol(0.0, lb, ub, 1, np.array([mi], dtype=np.int32), np.array([-1.0]))
        self.rxn_pos[rid] = self.n
        self.lb = np.append(self.lb, lb); self.ub = np.append(self.ub, ub)
        self.n += 1
        return rid

    def add_row(self, coefs: Dict[int, float], lo: float, hi: float) -> int:
        idx = np.array(sorted(coefs), dtype=np.int32)
        val = np.array([coefs[i] for i in idx], dtype=float)
        self.h.addRow(lo, hi, len(idx), idx, val)
        row = self.m0 + self.n_extra_rows
        self.n_extra_rows += 1
        return row

    def set_row_bounds(self, row: int, lo: float, hi: float) -> None:
        self.h.changeRowBounds(row, lo, hi)

    def set_bounds(self, cols: Sequence[int], lb: Optional[Sequence[float]] = None, ub: Optional[Sequence[float]] = None) -> None:
        cols = np.array(list(cols), dtype=np.int32)
        if len(cols) == 0:
            return
        lbv = self.lb[cols] if lb is None else np.array(lb, dtype=float)
        ubv = self.ub[cols] if ub is None else np.array(ub, dtype=float)
        self.h.changeColsBounds(len(cols), cols, lbv, ubv)
        self.lb[cols] = lbv; self.ub[cols] = ubv

    def set_objective(self, coefs: Dict[int, float], sense: str = "max") -> None:
        # clear previous objective then set new
        if getattr(self, "_obj_cols", None):
            prev = np.array(self._obj_cols, dtype=np.int32)
            self.h.changeColsCost(len(prev), prev, np.zeros(len(prev)))
        cols = np.array(sorted(coefs), dtype=np.int32)
        self.h.changeColsCost(len(cols), cols, np.array([coefs[i] for i in cols], dtype=float))
        self._obj_cols = cols.tolist()
        self._objective = dict(coefs)
        self._objective_sense = sense
        self.h.changeObjectiveSense(self.hs.ObjSense.kMaximize if sense == "max" else self.hs.ObjSense.kMinimize)

    @contextmanager
    def temporary_state(self):
        """Restore an IEM's bounds/objective and disable its added rows and demand sinks.

        runIEM_HH starts each IEM from modelO. Keeping added columns is cheap, but
        their flux must be zero outside the IEM that requested them.
        """
        old_n = self.n
        old_lb, old_ub = self.lb.copy(), self.ub.copy()
        old_rows = self.n_extra_rows
        old_objective = getattr(self, "_objective", {}).copy()
        old_sense = getattr(self, "_objective_sense", "max")
        try:
            yield
        finally:
            self.set_bounds(range(old_n), lb=old_lb, ub=old_ub)
            added = list(range(old_n, self.n))
            self.set_bounds(added, lb=[0.0] * len(added), ub=[0.0] * len(added))
            for row in range(self.m0 + old_rows, self.m0 + self.n_extra_rows):
                self.set_row_bounds(row, -np.inf, np.inf)
            self.set_objective(old_objective, old_sense)

    # --- solve ---------------------------------------------------------------------
    def solve(self) -> Tuple[str, float, np.ndarray, float]:
        t = time.time()
        # Interior point (with crossover) for every solve: on Harvey a cold IPM solve takes ~15-20 s whereas
        # dual-simplex warm starts after bound/objective changes were observed to take minutes.
        self.h.setOptionValue("solver", self.method)
        # HiGHS applies time_limit to the cumulative run time of the Highs object, so a persistent model
        # would start failing every solve once the total exceeds the limit; re-base it before each run.
        self.h.setOptionValue("time_limit", float(self.h.getRunTime()) + self.per_solve_time_limit)
        run_status = self.h.run()
        dt = time.time() - t
        self.n_solves += 1; self.solve_time += dt
        if run_status == self.hs.HighsStatus.kError:
            raise RuntimeError(
                f"HiGHS failed to run (threads={self.threads}); no optimum is available. "
                "HiGHS instances share a process-wide scheduler: use threads=0 or a separate "
                "process if their thread settings differ. Enable solver output to diagnose other run errors."
            )
        st = self.h.modelStatusToString(self.h.getModelStatus())
        sol = self.h.getSolution()
        x = np.array(sol.col_value) if sol.value_valid else np.full(self.n, np.nan)
        obj = float(self.h.getInfo().objective_function_value) if sol.value_valid else float("nan")
        return st, obj, x, dt


@dataclass
class BiomarkerResult:
    reaction: str
    expected: str            # Increased / Decreased / Unchanged
    healthy: float
    disease: float
    predicted: str
    correct: Optional[bool]  # None for an unscored label or unavailable/non-optimal solve
    status_healthy: str
    status_disease: str


@dataclass
class IEMResult:
    iem: str
    iem_reactions: List[str]
    vmax_healthy: float
    vmax_disease: float
    wb_objective_disease_feasible: bool
    biomarkers: List[BiomarkerResult]
    n_solves: int
    time_s: float
    notes: List[str] = field(default_factory=list)
    status: str = "not_run"


def match_reactions(wbm_rxns: np.ndarray, include: Sequence[str], exclude: Sequence[str] = ()) -> List[int]:
    idx = [i for i, r in enumerate(wbm_rxns) if any(p in r for p in include)]
    if exclude:
        idx = [i for i in idx if not any(p in wbm_rxns[i] for p in exclude)]
    return [i for i in idx if "Micro_" not in wbm_rxns[i]]


def apply_runiem_global_constraints(hw: HighsWBM) -> Dict[str, int]:
    """The unified reaction constraints set at the top of runIEM_HH.m before any IEM."""
    rx = hw.wbm.rxns
    irreversible = ["_ARGSL", "_GACMTRc", "_FUM", "_FUMm", "_HMR_7698", "_UAG4E", "_UDPG4E", "_GALT", "_G6PDH2c", "_G6PDH2r",
                    "_G6PDH2rer", "_GLUTCOADHm", "_r0541", "_ACOAD8m", "_RE2410C", "_RE2410N"]
    excluded = ["_FUMt", "_FUMAC", "_FUMS", "BBB"]
    irr = match_reactions(rx, irreversible, excluded)
    hw.set_bounds(irr, lb=[0.0] * len(irr))
    closed = match_reactions(rx, ["_r0784", "_r0463"])
    hw.set_bounds(closed, lb=[0.0] * len(closed), ub=[0.0] * len(closed))
    bile = [i for i, r in enumerate(rx) if r.startswith("BileDuct_EX_") and r.endswith("[bd]_[luSI]")]
    hw.set_bounds(bile, ub=[100.0] * len(bile))
    return {"set_irreversible": len(irr), "closed": len(closed), "bile_duct_ub_100": len(bile)}


def run_iem(hw: HighsWBM, iem: str, include_patterns: Sequence[str], exclude_patterns: Sequence[str],
            biomarkers: Sequence[Tuple[str, str]], min_flux_healthy: float = 1.0, tol: float = 1e-6,
            verbose: bool = True, demand_metabolites: Optional[Sequence[str]] = None) -> IEMResult:
    """biomarkers: list of (reaction id, expected label text such as 'Increased (blood)')."""
    if not np.isfinite(min_flux_healthy) or not 0 <= min_flux_healthy <= 1:
        raise ValueError("min_flux_healthy must be a finite fraction between 0 and 1")
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol must be finite and positive")
    started, n0 = time.time(), hw.n_solves
    with hw.temporary_state():
        result = _run_iem(hw, iem, include_patterns, exclude_patterns, biomarkers,
                          min_flux_healthy, tol, verbose, demand_metabolites)
    result.n_solves = hw.n_solves - n0
    result.time_s = time.time() - started
    return result


def _run_iem(hw: HighsWBM, iem: str, include_patterns: Sequence[str], exclude_patterns: Sequence[str],
             biomarkers: Sequence[Tuple[str, str]], min_flux_healthy: float, tol: float,
             verbose: bool, demand_metabolites: Optional[Sequence[str]]) -> IEMResult:
    t0 = time.time()
    n0 = hw.n_solves
    rx = hw.wbm.rxns
    iem_idx = match_reactions(rx, include_patterns, exclude_patterns)
    notes: List[str] = []
    result = IEMResult(iem=iem, iem_reactions=[rx[i] for i in iem_idx], vmax_healthy=float("nan"), vmax_disease=float("nan"),
                       wb_objective_disease_feasible=False, biomarkers=[], n_solves=0, time_s=0.0, notes=notes)
    if not iem_idx:
        notes.append("no reactions matched the IEM patterns")
        result.status = "no_reactions"
        return result
    # Some source blocks add demands beyond the biomarkers selected for scoring.
    required_mets = list(dict.fromkeys(list(demand_metabolites or []) +
                                      [rid[3:] for rid, _ in biomarkers if rid.startswith("DM_")]))
    for met in required_mets:
        rid = f"DM_{met}"
        if rid not in hw.rxn_pos:
            if met in hw.met_pos:
                hw.add_demand(met)
            else:
                notes.append(f"metabolite for {rid} not in model")
        elif hw.rxn_pos[rid] >= hw.wbm.n_rxns:
            # A previous IEM allocated this column, then disabled it on exit.
            hw.set_bounds([hw.rxn_pos[rid]], lb=[0.0], ub=[1000.0])
    saved_lb = hw.lb[iem_idx].copy(); saved_ub = hw.ub[iem_idx].copy()

    # 1. The original auxiliary variable constrains sum(IEM) to +/-100000.
    # An objective alone omits those bounds and changes the healthy reference.
    row = hw.add_row({i: 1.0 for i in iem_idx}, -1e5, 1e5)
    hw.set_objective({i: 1.0 for i in iem_idx}, "max")
    st, vmax, _, _ = hw.solve()
    result.vmax_healthy = vmax
    if st != "Optimal" or not np.isfinite(vmax) or abs(vmax) <= tol:
        result.status = "inactive" if st == "Optimal" and np.isfinite(vmax) and abs(vmax) <= tol else "healthy_solve_failed"
        notes.append(f"healthy IEM-flux maximisation: status {st}, vmax {vmax:.3g} -> IEM not simulated (as in checkIEM_WBM)")
        return result
    # 2. healthy constraint row: sum(IEM) >= min_flux_healthy * vmax (truncated to 6 decimals)
    lo = math.floor(min_flux_healthy * vmax * 1e6) / 1e6 if vmax > 0 else math.ceil(min_flux_healthy * vmax * 1e6) / 1e6
    hw.set_row_bounds(row, lo, 1e5)

    def healthy_state():
        hw.set_bounds(iem_idx, lb=saved_lb, ub=saved_ub)
        hw.set_row_bounds(row, lo, 1e5)

    def disease_state():
        hw.set_bounds(iem_idx, lb=[0.0] * len(iem_idx), ub=[0.0] * len(iem_idx))
        hw.set_row_bounds(row, 0.0, 1e5)

    # disease: summed IEM flux (should be 0) and whole-body objective feasibility
    disease_state()
    st, vd, _, _ = hw.solve()
    result.vmax_disease = vd if st == "Optimal" and np.isfinite(vd) else float("nan")
    if st != "Optimal" or not np.isfinite(vd):
        result.status = "disease_infeasible" if st == "Infeasible" else "disease_solve_failed"
        notes.append(f"disease IEM-flux solve unavailable ({st}); biomarkers not simulated")
        return result
    wb = hw.rxn_pos.get("Whole_body_objective_rxn")
    if wb is not None:
        hw.set_objective({wb: 1.0}, "max")
        st, fw, _, _ = hw.solve()
        result.wb_objective_disease_feasible = bool(st == "Optimal" and np.isfinite(fw))
        if not result.wb_objective_disease_feasible:
            result.status = "wb_infeasible" if st == "Infeasible" else "wb_solve_failed"
            notes.append(f"whole-body objective infeasible in disease state ({st}); biomarkers not simulated")
            healthy_state(); hw.set_row_bounds(row, -np.inf, np.inf)
            result.n_solves = hw.n_solves - n0; result.time_s = time.time() - t0
            return result
    else:
        result.status = "missing_wb_objective"
        notes.append("Whole_body_objective_rxn is absent; whole-body IEM protocol cannot be checked")
        return result
    # 3. biomarkers
    for rid, label in biomarkers:
        if rid not in hw.rxn_pos:
            result.biomarkers.append(BiomarkerResult(rid, _expected(label), float("nan"), float("nan"), "NA", None, "absent", "absent"))
            continue
        j = hw.rxn_pos[rid]
        old_ub = hw.ub[j]
        hw.set_bounds([j], ub=[1e5])
        hw.set_objective({j: 1.0}, "max")
        healthy_state(); sth, fh, _, _ = hw.solve()
        disease_state(); std, fd, _, _ = hw.solve()
        valid_h = sth == "Optimal" and np.isfinite(fh)
        valid_d = std == "Optimal" and np.isfinite(fd)
        fh = (0.0 if abs(fh) <= tol else fh) if valid_h else float("nan")
        fd = (0.0 if abs(fd) <= tol else fd) if valid_d else float("nan")
        hw.set_bounds([j], ub=[old_ub])
        exp = _expected(label)
        if valid_h and valid_d:
            diff = fd - fh
            pred = "Increased" if diff > tol else "Decreased" if diff < -tol else "Unchanged"
            correct = None if exp == "Unchanged" else (pred == exp)
        else:
            pred, correct = "NA", None
            notes.append(f"{rid}: biomarker prediction unavailable (healthy {sth}, disease {std})")
        result.biomarkers.append(BiomarkerResult(rid, exp, fh, fd, pred, correct, sth, std))
        if verbose:
            print(f"    {iem:8s} {rid:28s} healthy={fh:11.4g} disease={fd:11.4g} pred={pred:9s} expected={exp:9s} {'OK' if correct else ('--' if correct is None else 'X')}", flush=True)
    healthy_state(); hw.set_row_bounds(row, -np.inf, np.inf)   # leave the model as found (row inert)
    result.n_solves = hw.n_solves - n0; result.time_s = time.time() - t0
    result.status = "partial" if any(b.status_healthy not in {"Optimal", "absent"} or
                                    b.status_disease not in {"Optimal", "absent"} or
                                    (b.status_healthy != "absent" and (not np.isfinite(b.healthy) or not np.isfinite(b.disease)))
                                    for b in result.biomarkers) else "complete"
    return result


def _expected(label: str) -> str:
    l = label.lower()
    if "incre" in l:
        return "Increased"
    if "decre" in l:
        return "Decreased"
    return "Unchanged"
