"""Whole-body metabolic models (Harvey/Harvetta) in Python: load from the COBRA Toolbox .mat
format including coupling constraints (C, d, dsense), build the LP directly, solve with open
solvers, and certify solutions by recomputing residuals in double precision.

The LP is
    min/max  c^T v
    s.t.     S v  (csense)  b        equality rows (metabolite balances)
             C v  (dsense)  d        coupling rows (organ-level constraints)
             lb <= v <= ub
Nothing goes through cobrapy here: at 81k variables and 160k rows the point is to see what the
solvers do, without a modelling layer in between.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy.io as sio
import scipy.sparse as sp


@dataclass
class WBM:
    name: str
    rxns: np.ndarray
    mets: np.ndarray
    S: sp.csc_matrix
    b: np.ndarray
    csense: np.ndarray
    C: sp.csc_matrix
    d: np.ndarray
    dsense: np.ndarray
    ctrs: np.ndarray
    lb: np.ndarray
    ub: np.ndarray
    c: np.ndarray
    osense: str
    meta: Dict[str, str] = field(default_factory=dict)

    @property
    def n_rxns(self) -> int:
        return len(self.rxns)

    def rxn_index(self, rid: str) -> int:
        idx = np.where(self.rxns == rid)[0]
        if len(idx) == 0:
            raise KeyError(rid)
        return int(idx[0])

    def find_rxns(self, pattern: str) -> List[int]:
        return [i for i, r in enumerate(self.rxns) if pattern in r]


def _cellstr(a) -> np.ndarray:
    return np.array([str(x[0]) if hasattr(x, "__len__") and len(x) else "" for x in a[:, 0]])


def load_wbm(path: str, varname: Optional[str] = None) -> WBM:
    t = time.time()
    d = sio.loadmat(path, squeeze_me=False, struct_as_record=False)
    keys = [k for k in d if not k.startswith("__")]
    m = d[varname or keys[0]][0, 0]
    S = sp.csc_matrix(m.S)
    C = sp.csc_matrix(m.C) if hasattr(m, "C") else sp.csc_matrix((0, S.shape[1]))
    wbm = WBM(
        name=str(m.modelID[0]) if hasattr(m, "modelID") else keys[0],
        rxns=_cellstr(m.rxns), mets=_cellstr(m.mets), S=S,
        b=np.asarray(m.b, dtype=float).ravel(), csense=np.asarray(m.csense).astype(str).ravel(),
        C=C, d=np.asarray(m.d, dtype=float).ravel() if hasattr(m, "d") else np.zeros(C.shape[0]),
        dsense=np.asarray(m.dsense).astype(str).ravel() if hasattr(m, "dsense") else np.array([]),
        ctrs=_cellstr(m.ctrs) if hasattr(m, "ctrs") else np.array([]),
        lb=np.asarray(m.lb, dtype=float).ravel(), ub=np.asarray(m.ub, dtype=float).ravel(),
        c=np.asarray(m.c, dtype=float).ravel(),
        osense=str(m.osenseStr[0]) if hasattr(m, "osenseStr") else "max",
        meta={"file": path, "load_s": f"{time.time()-t:.1f}", "version": str(getattr(m, "version", [""])[0]),
              "status": str(getattr(m, "status", [""])[0]), "sex": str(getattr(m, "sex", [""])[0])},
    )
    return wbm


def stacked_constraints(wbm: WBM) -> Tuple[sp.csc_matrix, np.ndarray, np.ndarray]:
    """Return A, row_lower, row_upper for [S; C] with senses converted to ranges."""
    A = sp.vstack([wbm.S, wbm.C]).tocsc()
    inf = np.inf
    rl = np.empty(A.shape[0]); ru = np.empty(A.shape[0])
    for i, (s, bi) in enumerate(zip(wbm.csense, wbm.b)):
        rl[i], ru[i] = {"E": (bi, bi), "L": (-inf, bi), "G": (bi, inf)}[s]
    off = wbm.S.shape[0]
    for i, (s, di) in enumerate(zip(wbm.dsense, wbm.d)):
        rl[off + i], ru[off + i] = {"E": (di, di), "L": (-inf, di), "G": (di, inf)}[s]
    return A, rl, ru


@dataclass
class SolveResult:
    solver: str
    method: str
    status: str
    objective: float
    time_s: float
    x: Optional[np.ndarray] = None
    certificate: Dict[str, float] = field(default_factory=dict)
    info: Dict[str, str] = field(default_factory=dict)


def certify(wbm: WBM, x: np.ndarray) -> Dict[str, float]:
    """Recompute constraint residuals in double precision — a solver-independent check."""
    Sx = wbm.S @ x
    eq_res = np.abs(Sx - wbm.b)
    out = {"max_abs_S_residual": float(eq_res.max()), "n_S_rows_over_1e-6": int((eq_res > 1e-6).sum())}
    if wbm.C.shape[0]:
        Cx = wbm.C @ x
        viol = np.zeros_like(Cx)
        L = wbm.dsense == "L"; G = wbm.dsense == "G"; E = wbm.dsense == "E"
        viol[L] = np.maximum(Cx[L] - wbm.d[L], 0)
        viol[G] = np.maximum(wbm.d[G] - Cx[G], 0)
        viol[E] = np.abs(Cx[E] - wbm.d[E])
        out["max_coupling_violation"] = float(viol.max())
        out["n_coupling_rows_over_1e-6"] = int((viol > 1e-6).sum())
    bviol = np.maximum(np.maximum(wbm.lb - x, x - wbm.ub), 0)
    out["max_bound_violation"] = float(bviol.max())
    out["n_bounds_over_1e-6"] = int((bviol > 1e-6).sum())
    return out


def solve_highs(wbm: WBM, objective: Optional[np.ndarray] = None, sense: str = "max", method: str = "simplex",
                presolve: str = "choose", feas_tol: float = 1e-7, opt_tol: float = 1e-7, time_limit: float = 3600.0,
                threads: int = 1, extra_options: Optional[Dict[str, object]] = None) -> SolveResult:
    import highspy
    A, rl, ru = stacked_constraints(wbm)
    n, mrows = A.shape[1], A.shape[0]
    c = wbm.c.copy() if objective is None else np.asarray(objective, dtype=float)
    lp = highspy.HighsLp()
    lp.num_col_ = n; lp.num_row_ = mrows
    lp.col_cost_ = c
    lp.col_lower_ = wbm.lb.copy(); lp.col_upper_ = wbm.ub.copy()
    lp.row_lower_ = rl; lp.row_upper_ = ru
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A.indptr.astype(np.int32)
    lp.a_matrix_.index_ = A.indices.astype(np.int32)
    lp.a_matrix_.value_ = A.data.astype(float)
    lp.sense_ = highspy.ObjSense.kMaximize if sense == "max" else highspy.ObjSense.kMinimize
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("solver", method)                # 'simplex' | 'ipm' | 'pdlp'
    h.setOptionValue("presolve", presolve)
    h.setOptionValue("primal_feasibility_tolerance", feas_tol)
    h.setOptionValue("dual_feasibility_tolerance", opt_tol)
    h.setOptionValue("time_limit", time_limit)
    h.setOptionValue("threads", threads)
    for k, v in (extra_options or {}).items():
        h.setOptionValue(k, v)
    h.passModel(lp)
    t = time.time()
    h.run()
    dt = time.time() - t
    status = h.modelStatusToString(h.getModelStatus())
    info = h.getInfo()
    sol = h.getSolution()
    x = np.array(sol.col_value) if sol.value_valid else None
    obj = float(info.objective_function_value) if x is not None else float("nan")
    res = SolveResult(solver=f"highs {highspy.HIGHS_VERSION_MAJOR}.{highspy.HIGHS_VERSION_MINOR}.{highspy.HIGHS_VERSION_PATCH}" if hasattr(highspy, "HIGHS_VERSION_MAJOR") else "highs",
                      method=method, status=status, objective=obj, time_s=dt, x=x,
                      info={"simplex_iterations": str(info.simplex_iteration_count), "ipm_iterations": str(info.ipm_iteration_count),
                            "max_primal_infeasibility": str(info.max_primal_infeasibility), "max_dual_infeasibility": str(info.max_dual_infeasibility),
                            "presolve": presolve, "feas_tol": str(feas_tol)})
    if x is not None:
        res.certificate = certify(wbm, x)
    return res


def solve_glpk(wbm: WBM, objective: Optional[np.ndarray] = None, sense: str = "max", method: str = "simplex",
               presolve: bool = True, time_limit_ms: int = 3_600_000, feas_tol: float = 1e-7) -> SolveResult:
    import swiglpk as glp
    A, rl, ru = stacked_constraints(wbm)
    n, mrows = A.shape[1], A.shape[0]
    c = wbm.c.copy() if objective is None else np.asarray(objective, dtype=float)
    P = glp.glp_create_prob()
    glp.glp_set_obj_dir(P, glp.GLP_MAX if sense == "max" else glp.GLP_MIN)
    glp.glp_add_rows(P, mrows); glp.glp_add_cols(P, n)
    for i in range(mrows):
        lo, up = rl[i], ru[i]
        if np.isfinite(lo) and np.isfinite(up):
            glp.glp_set_row_bnds(P, i + 1, glp.GLP_FX if lo == up else glp.GLP_DB, float(lo), float(up))
        elif np.isfinite(lo):
            glp.glp_set_row_bnds(P, i + 1, glp.GLP_LO, float(lo), 0.0)
        elif np.isfinite(up):
            glp.glp_set_row_bnds(P, i + 1, glp.GLP_UP, 0.0, float(up))
        else:
            glp.glp_set_row_bnds(P, i + 1, glp.GLP_FR, 0.0, 0.0)
    for j in range(n):
        lo, up = wbm.lb[j], wbm.ub[j]
        glp.glp_set_col_bnds(P, j + 1, glp.GLP_FX if lo == up else glp.GLP_DB, float(lo), float(up))
        if c[j] != 0:
            glp.glp_set_obj_coef(P, j + 1, float(c[j]))
    coo = A.tocoo()
    nnz = coo.nnz
    ia = glp.intArray(nnz + 1); ja = glp.intArray(nnz + 1); ar = glp.doubleArray(nnz + 1)
    for k, (i, j, v) in enumerate(zip(coo.row, coo.col, coo.data), start=1):
        ia[k] = int(i) + 1; ja[k] = int(j) + 1; ar[k] = float(v)
    glp.glp_load_matrix(P, nnz, ia, ja, ar)
    t = time.time()
    if method == "interior":
        parm = glp.glp_iptcp(); glp.glp_init_iptcp(parm); parm.msg_lev = glp.GLP_MSG_OFF
        ret = glp.glp_interior(P, parm)
        st = glp.glp_ipt_status(P)
        status = {glp.GLP_OPT: "Optimal", glp.GLP_INFEAS: "Infeasible", glp.GLP_NOFEAS: "NoFeasible", glp.GLP_UNDEF: "Undefined"}.get(st, str(st))
        x = np.array([glp.glp_ipt_col_prim(P, j + 1) for j in range(n)]) if st == glp.GLP_OPT else None
        obj = glp.glp_ipt_obj_val(P)
    else:
        parm = glp.glp_smcp(); glp.glp_init_smcp(parm); parm.msg_lev = glp.GLP_MSG_OFF
        parm.presolve = glp.GLP_ON if presolve else glp.GLP_OFF
        parm.tm_lim = time_limit_ms; parm.tol_bnd = feas_tol; parm.tol_dj = feas_tol
        parm.meth = glp.GLP_DUALP if method == "dual" else glp.GLP_PRIMAL
        ret = glp.glp_simplex(P, parm)
        st = glp.glp_get_status(P)
        status = {glp.GLP_OPT: "Optimal", glp.GLP_FEAS: "Feasible", glp.GLP_INFEAS: "Infeasible", glp.GLP_NOFEAS: "NoFeasible",
                  glp.GLP_UNBND: "Unbounded", glp.GLP_UNDEF: "Undefined"}.get(st, str(st))
        x = np.array([glp.glp_get_col_prim(P, j + 1) for j in range(n)]) if st in (glp.GLP_OPT, glp.GLP_FEAS) else None
        obj = glp.glp_get_obj_val(P)
    dt = time.time() - t
    glp.glp_delete_prob(P)
    res = SolveResult(solver=f"glpk {glp.glp_version()}", method=method, status=status, objective=float(obj), time_s=dt, x=x,
                      info={"return_code": str(ret), "presolve": str(presolve)})
    if x is not None:
        res.certificate = certify(wbm, x)
    return res


def unit_objective(wbm: WBM, rxn_id: str) -> np.ndarray:
    c = np.zeros(wbm.n_rxns)
    c[wbm.rxn_index(rxn_id)] = 1.0
    return c
