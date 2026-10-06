"""Gurobi backend for the whole-body IEM protocol, with the same interface as wbm_iem.HighsWBM.

It holds one persistent gurobipy model, so bound and objective changes keep the last basis for warm starts.
Added rows (the healthy pin) use the Toolbox's own formulation: an auxiliary variable s with
sum(IEM fluxes) - s = 0, whose bounds are the row bounds. Model statuses are reported with HiGHS's names
("Optimal", "Infeasible", ...) so the protocol code is shared unchanged.

Needs a full Gurobi licence: the size-limited licence that ships with gurobipy refuses whole-body models.
"""
from __future__ import annotations

import time
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from .wbm import WBM, stacked_constraints
from .wbm_iem import HighsWBM


class GurobiWBM:
    STATUS = {2: "Optimal", 3: "Infeasible", 4: "Primal infeasible or unbounded", 5: "Unbounded",
              9: "Time limit reached", 12: "Numeric error", 13: "Suboptimal", 11: "Interrupted"}

    def __init__(self, wbm: WBM, feas_tol: float = 1e-7, opt_tol: float = 1e-7, threads: int = 0,
                 time_limit: float = 1800.0):
        import gurobipy as gp
        from gurobipy import GRB
        self.gp, self.GRB = gp, GRB
        self.wbm = wbm
        self.n = wbm.n_rxns
        A, rl, ru = stacked_constraints(wbm)
        self.m0 = A.shape[0]
        env = gp.Env(empty=True)
        env.setParam("OutputFlag", 0)
        env.start()
        self.env = env
        model = gp.Model("wbm", env=env)
        model.Params.FeasibilityTol = feas_tol
        model.Params.OptimalityTol = opt_tol
        model.Params.Threads = threads
        model.Params.TimeLimit = time_limit
        x = model.addMVar(self.n, lb=self._inf(wbm.lb), ub=self._inf(wbm.ub), name="v")
        eq = rl == ru
        le = np.isneginf(rl) & np.isfinite(ru)
        ge = np.isfinite(rl) & np.isposinf(ru)
        if not np.all(eq | le | ge):
            raise ValueError("ranged or free rows in [S; C] are not expected")
        sense = np.where(eq, "=", np.where(le, "<", ">"))
        rhs = np.where(ge, rl, ru)
        mc = model.addMConstr(A.tocsr(), x, sense, rhs)
        model.ModelSense = GRB.MAXIMIZE
        model.update()
        self.model = model
        self.cols = x.tolist()
        self.rows = mc.tolist()
        self.extra = []   # (constraint, auxiliary variable) per added row
        self.threads = threads
        self.per_solve_time_limit = time_limit
        self.rxn_pos = {r: i for i, r in enumerate(wbm.rxns)}
        self.met_pos = {m: i for i, m in enumerate(wbm.mets)}
        self.lb = wbm.lb.copy(); self.ub = wbm.ub.copy()
        self.n_extra_rows = 0
        self.n_solves = 0
        self.solve_time = 0.0
        self.method = "ipm"
        self._obj_cols = []
        self.last_info = {}

    def _inf(self, v):
        v = np.array(v, dtype=float)
        return np.clip(v, -self.GRB.INFINITY, self.GRB.INFINITY)

    # --- structure -----------------------------------------------------------------
    def add_demand(self, met_id: str, lb: float = 0.0, ub: float = 1000.0) -> str:
        rid = f"DM_{met_id}"
        if rid in self.rxn_pos:
            return rid
        mi = self.met_pos[met_id]
        var = self.model.addVar(lb=lb, ub=ub, obj=0.0, name=rid, column=self.gp.Column([-1.0], [self.rows[mi]]))
        self.model.update()
        self.cols.append(var)
        self.rxn_pos[rid] = self.n
        self.lb = np.append(self.lb, lb); self.ub = np.append(self.ub, ub)
        self.n += 1
        return rid

    def add_row(self, coefs: Dict[int, float], lo: float, hi: float) -> int:
        lo, hi = self._inf([lo, hi])
        s = self.model.addVar(lb=lo, ub=hi, name=f"aux{self.n_extra_rows}")
        expr = self.gp.LinExpr([coefs[i] for i in sorted(coefs)], [self.cols[i] for i in sorted(coefs)])
        c = self.model.addLConstr(expr - s, self.GRB.EQUAL, 0.0)
        self.model.update()
        self.extra.append((c, s))
        row = self.m0 + self.n_extra_rows
        self.n_extra_rows += 1
        return row

    def set_row_bounds(self, row: int, lo: float, hi: float) -> None:
        _, s = self.extra[row - self.m0]
        lo, hi = self._inf([lo, hi])
        s.LB = lo; s.UB = hi

    def set_bounds(self, cols: Sequence[int], lb: Optional[Sequence[float]] = None, ub: Optional[Sequence[float]] = None) -> None:
        cols = np.array(list(cols), dtype=np.int64)
        if len(cols) == 0:
            return
        lbv = self.lb[cols] if lb is None else np.array(lb, dtype=float)
        ubv = self.ub[cols] if ub is None else np.array(ub, dtype=float)
        vars_ = [self.cols[i] for i in cols]
        self.model.setAttr("LB", vars_, self._inf(lbv).tolist())
        self.model.setAttr("UB", vars_, self._inf(ubv).tolist())
        self.lb[cols] = lbv; self.ub[cols] = ubv

    def set_objective(self, coefs: Dict[int, float], sense: str = "max") -> None:
        if self._obj_cols:
            self.model.setAttr("Obj", [self.cols[i] for i in self._obj_cols], [0.0] * len(self._obj_cols))
        cols = sorted(coefs)
        self.model.setAttr("Obj", [self.cols[i] for i in cols], [float(coefs[i]) for i in cols])
        self._obj_cols = list(cols)
        self._objective = dict(coefs)
        self._objective_sense = sense
        self.model.ModelSense = self.GRB.MAXIMIZE if sense == "max" else self.GRB.MINIMIZE

    temporary_state = HighsWBM.temporary_state

    # --- solve ---------------------------------------------------------------------
    def solve(self) -> Tuple[str, float, np.ndarray, float]:
        # "ipm": barrier with crossover (Gurobi's default crossover), as the HiGHS protocol runs.
        # "primal"/"dual": simplex from the basis of the previous solve (LPWarmStart default).
        # "concurrent": primal and dual simplex (warm) and barrier in parallel; the first to finish is used.
        self.model.Params.Method = {"ipm": 2, "primal": 0, "dual": 1, "concurrent": 3}[self.method]
        self.model.Params.TimeLimit = self.per_solve_time_limit
        t = time.time()
        self.model.optimize()
        dt = time.time() - t
        self.n_solves += 1; self.solve_time += dt
        status = self.model.Status
        st = self.STATUS.get(status, f"Gurobi status {status}")
        if status == self.GRB.OPTIMAL:
            obj = float(self.model.ObjVal)
            x = np.array(self.model.getAttr("X", self.cols))
        else:
            obj = float("nan"); x = np.full(self.n, np.nan)
        self.last_info = {"method": self.method, "simplex_iterations": int(self.model.IterCount),
                          "ipm_iterations": int(self.model.BarIterCount)}
        return st, obj, x, dt
