"""Minimal gap-filling of a draft model against a defined medium, with a universal reaction database.

MILP (HiGHS via scipy.optimize.milp): add the smallest weighted set of reactions from the universe
so that the biomass reaction can carry at least `min_growth` on the given medium.  Exchange
reactions are never added (the medium is the medium); every added reaction is recorded with a
'gapfill' annotation so that downstream analyses can tell them apart from annotated reactions.

This is the classic parsimonious gap-fill (Reed et al. 2006; CarveMe's `gapfill`), reimplemented here
so it runs with open-source solvers and so that the added set is explicit in the benchmark card.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import cobra
import numpy as np
import scipy.sparse as sp
from scipy.optimize import LinearConstraint

from .media import Medium, apply_medium


@dataclass
class GapfillResult:
    added_reactions: List[str]
    growth_before: float
    growth_after: float
    objective: float
    status: str
    n_candidates: int
    seconds: float
    weights_used: Dict[str, float] = field(default_factory=dict)


def gapfill(model: cobra.Model, universe: cobra.Model, medium: Medium, extra_uptakes: Optional[Dict[str, float]] = None,
            biomass_id: str = "Growth", min_growth: float = 0.05, transport_weight: float = 1.5,
            time_limit_s: float = 600.0, mip_gap: float = 0.0, forbid_o2_production: bool = True,
            verbose: bool = True) -> GapfillResult:
    """Return the minimal reaction set from `universe` that lets `model` grow on `medium` (+ extra uptakes).

    The model is not modified; use `apply_gapfill` to add the reactions.
    """
    t0 = time.time()
    m = model.copy()
    apply_medium(m, medium, close_all=True)
    for ex_id, lb in (extra_uptakes or {}).items():
        if ex_id in m.reactions:
            m.reactions.get_by_id(ex_id).lower_bound = lb
    growth_before = _nan0(m.slim_optimize())
    if growth_before >= min_growth:
        return GapfillResult([], growth_before, growth_before, 0.0, "not_needed", 0, time.time() - t0)

    model_rxn_ids = {r.id for r in m.reactions}
    cands = [r for r in universe.reactions if r.id not in model_rxn_ids and not r.id.startswith("EX_")
             and not r.id.startswith("sink_") and not r.id.startswith("DM_")]
    # metabolite index over model + candidate metabolites
    mets: Dict[str, int] = {met.id: i for i, met in enumerate(m.metabolites)}
    for r in cands:
        for met in r.metabolites:
            if met.id not in mets:
                mets[met.id] = len(mets)
    n_m = len(m.reactions); n_c = len(cands); n_v = n_m + n_c
    rows, cols, vals = [], [], []
    lb = np.zeros(n_v + n_c); ub = np.zeros(n_v + n_c)
    for j, r in enumerate(m.reactions):
        lb[j], ub[j] = r.lower_bound, r.upper_bound
        for met, coef in r.metabolites.items():
            rows.append(mets[met.id]); cols.append(j); vals.append(coef)
    big = 100.0   # candidate flux bound during gap-filling (tight big-M limits integrality leaks)
    for k, r in enumerate(cands):
        j = n_m + k
        rlb, rub = r.lower_bound, r.upper_bound
        if forbid_o2_production:
            # no candidate reaction may run in a direction that releases molecular oxygen: the universe marks
            # some oxygenases reversible, and a parsimonious fill otherwise 'generates' O2 in anaerobic media
            o2 = [coef for met, coef in r.metabolites.items() if met.id.startswith("o2_")]
            if any(c > 0 for c in o2):
                rub = 0.0
            if any(c < 0 for c in o2):
                rlb = min(rlb, 0.0) if rlb < 0 else 0.0
                rlb = 0.0
        lb[j], ub[j] = (-big if rlb < 0 else 0.0), (big if rub > 0 else 0.0)
        for met, coef in r.metabolites.items():
            rows.append(mets[met.id]); cols.append(j); vals.append(coef)
    # binaries y_k in [0,1]
    lb[n_v:] = 0.0; ub[n_v:] = 1.0
    S = sp.csr_matrix((vals, (rows, cols)), shape=(len(mets), n_v + n_c))
    constraints = [LinearConstraint(S, np.zeros(len(mets)), np.zeros(len(mets)))]
    # v_k - ub_k * y_k <= 0 ; v_k - lb_k * y_k >= 0
    r2, c2, v2 = [], [], []
    lo2, hi2 = [], []
    row = 0
    for k in range(n_c):
        j = n_m + k
        if ub[j] > 0:
            r2 += [row, row]; c2 += [j, n_v + k]; v2 += [1.0, -ub[j]]; lo2.append(-np.inf); hi2.append(0.0); row += 1
        if lb[j] < 0:
            r2 += [row, row]; c2 += [j, n_v + k]; v2 += [1.0, -lb[j]]; lo2.append(0.0); hi2.append(np.inf); row += 1
    if row:
        A2 = sp.csr_matrix((v2, (r2, c2)), shape=(row, n_v + n_c))
        constraints.append(LinearConstraint(A2, np.array(lo2), np.array(hi2)))
    # growth >= min_growth
    bidx = m.reactions.index(m.reactions.get_by_id(biomass_id))
    g = sp.csr_matrix(([1.0], ([0], [bidx])), shape=(1, n_v + n_c))
    constraints.append(LinearConstraint(g, np.array([min_growth]), np.array([np.inf])))
    # objective: weighted count of added reactions (transport/periplasmic reactions slightly dearer)
    c = np.zeros(n_v + n_c)
    weights = {}
    for k, r in enumerate(cands):
        w = transport_weight if _is_transport(r) else 1.0
        c[n_v + k] = w; weights[r.id] = w
    if verbose:
        print(f"  gapfill MILP: {n_v + n_c} vars ({n_c} binaries), {len(mets)} metabolites; solving...", flush=True)
    x, status = _solve_highs(c, constraints, lb, ub, n_v, time_limit_s, mip_gap)
    if x is None:
        return GapfillResult([], growth_before, growth_before, float("nan"), f"failed: {status}", n_c, time.time() - t0)
    y = x[n_v:]
    added = [cands[k].id for k in range(n_c) if y[k] > 0.5]
    res_fun = float(np.dot(c, x))
    # verify in cobra
    m2 = model.copy()
    m2.add_reactions([universe.reactions.get_by_id(a).copy() for a in added])
    apply_medium(m2, medium, close_all=True)
    for ex_id, lbv in (extra_uptakes or {}).items():
        if ex_id in m2.reactions:
            m2.reactions.get_by_id(ex_id).lower_bound = lbv
    growth_after = _nan0(m2.slim_optimize())
    return GapfillResult(added, growth_before, growth_after, res_fun, status, n_c, time.time() - t0,
                         {a: weights[a] for a in added})


def _solve_highs(c, constraints, lb, ub, n_v, time_limit_s, mip_gap):
    """Solve min c'x s.t. constraints, bounds, x[n_v:] binary, with HiGHS and tight integrality tolerances."""
    import highspy
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", float(time_limit_s))
    h.setOptionValue("mip_rel_gap", float(mip_gap))
    h.setOptionValue("mip_feasibility_tolerance", 1e-9)
    h.setOptionValue("primal_feasibility_tolerance", 1e-9)
    n = len(c)
    inf = h.getInfinity()
    lbv = np.where(np.isfinite(lb), lb, -inf); ubv = np.where(np.isfinite(ub), ub, inf)
    h.addVars(n, lbv, ubv)
    h.changeColsCost(n, np.arange(n, dtype=np.int32), np.asarray(c, dtype=float))
    A = sp.vstack([cc.A if sp.issparse(cc.A) else sp.csr_matrix(cc.A) for cc in constraints]).tocsr()
    rlo = np.concatenate([np.atleast_1d(cc.lb) for cc in constraints]).astype(float)
    rhi = np.concatenate([np.atleast_1d(cc.ub) for cc in constraints]).astype(float)
    rlo = np.where(np.isfinite(rlo), rlo, -inf); rhi = np.where(np.isfinite(rhi), rhi, inf)
    h.addRows(A.shape[0], rlo, rhi, A.nnz, A.indptr.astype(np.int32), A.indices.astype(np.int32), A.data.astype(float))
    idx = np.arange(n_v, n, dtype=np.int32)
    h.changeColsIntegrality(len(idx), idx, np.array([highspy.HighsVarType.kInteger] * len(idx)))
    h.run()
    st = h.getModelStatus()
    status = h.modelStatusToString(st)
    sol = h.getSolution()
    x = np.array(sol.col_value) if len(sol.col_value) == n else None
    if x is None or status.lower() in ("infeasible", "unbounded"):
        return None, status
    return x, status


def apply_gapfill(model: cobra.Model, universe: cobra.Model, added: List[str], note: str = "") -> cobra.Model:
    m = model.copy()
    rx = []
    for a in added:
        r = universe.reactions.get_by_id(a).copy()
        r.annotation = dict(r.annotation); r.annotation["gapfill"] = note or "added by gembench.gapfill"
        r.gene_reaction_rule = ""
        rx.append(r)
    m.add_reactions(rx)
    return m


def _is_transport(r: cobra.Reaction) -> bool:
    comps = {met.compartment for met in r.metabolites}
    return len(comps) > 1


def _nan0(x) -> float:
    return 0.0 if (x is None or (isinstance(x, float) and np.isnan(x))) else float(x)
