"""Minimal gap-filling of a draft model against a defined medium, with a universal reaction database.

MILP (HiGHS): add the smallest weighted set of reactions from the universe
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
    rejected_energy_cycles: List[List[str]] = field(default_factory=list)   # added sets cut off for creating an EGC


def gapfill(model: cobra.Model, universe: cobra.Model, medium: Medium, extra_uptakes: Optional[Dict[str, float]] = None,
            biomass_id: str = "Growth", min_growth: float = 0.05, transport_weight: float = 1.5,
            time_limit_s: float = 600.0, mip_gap: float = 0.0, forbid_o2_production: bool = True,
            forbid_energy_cycles: bool = True, max_rounds: int = 6, verbose: bool = True) -> GapfillResult:
    """Return the minimal reaction set from `universe` that lets `model` grow on `medium` (+ extra uptakes).

    The model is not modified; use `apply_gapfill` to add the reactions.
    With `forbid_energy_cycles`, a solution whose added set creates an energy-generating cycle (ATP, NAD(P)H, quinol or a
    proton gradient from nothing with every boundary reaction closed; gembench.checks.energy_from_nothing) is cut off with a
    no-good constraint and the MILP is re-solved, up to `max_rounds` times; the rejected sets are reported.
    """
    if not np.isfinite(min_growth) or min_growth <= 0:
        raise ValueError("min_growth must be finite and positive")
    if not np.isfinite(transport_weight) or transport_weight <= 0:
        raise ValueError("transport_weight must be finite and positive")
    if max_rounds < 1 or not np.isfinite(time_limit_s) or time_limit_s <= 0:
        raise ValueError("max_rounds and time_limit_s must be positive, with a finite time limit")
    t0 = time.time()
    m = model.copy()
    m.objective = m.reactions.get_by_id(biomass_id)
    m.objective_direction = "max"
    apply_medium(m, medium, close_all=True)
    for ex_id, lb in (extra_uptakes or {}).items():
        if ex_id in m.reactions:
            m.reactions.get_by_id(ex_id).lower_bound = lb
    growth_before = _growth(m, allow_infeasible=True)
    if forbid_energy_cycles:
        from .checks import energy_from_nothing
        egc = {k: v for k, v in energy_from_nothing(m).items() if v > 1e-6}
        if egc:
            return GapfillResult([], growth_before, growth_before, float("nan"),
                                 f"failed: input model has energy-generating cycles {egc}", 0, time.time() - t0)
    if growth_before >= min_growth:
        return GapfillResult([], growth_before, growth_before, 0.0, "not_needed", 0, time.time() - t0)

    model_rxn_ids = {r.id for r in m.reactions}
    cands = [r for r in universe.reactions if r.id not in model_rxn_ids and not r.boundary and r.metabolites
             and not r.id.startswith(("EX_", "sink_", "DM_"))
             and _candidate_bounds(r, forbid_o2_production) is not None]
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
    candidate_bounds = []
    for k, r in enumerate(cands):
        j = n_m + k
        rlb, rub = _candidate_bounds(r, forbid_o2_production)
        if not np.isfinite(rlb) or not np.isfinite(rub):
            raise ValueError(f"Gapfill candidate {r.id} needs finite flux bounds")
        candidate_bounds.append((rlb, rub))
        # Unselected candidates carry zero; selected candidates keep their actual
        # capacity and any mandatory flux, enforced by the indicator constraints.
        lb[j], ub[j] = min(0.0, rlb), max(0.0, rub)
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
        rlb, rub = candidate_bounds[k]
        r2 += [row, row]; c2 += [j, n_v + k]; v2 += [1.0, -rub]; lo2.append(-np.inf); hi2.append(0.0); row += 1
        r2 += [row, row]; c2 += [j, n_v + k]; v2 += [1.0, -rlb]; lo2.append(0.0); hi2.append(np.inf); row += 1
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
    rejected: List[List[str]] = []
    for _round in range(max_rounds):
        remaining = time_limit_s - (time.time() - t0)
        if remaining <= 0:
            return GapfillResult([], growth_before, growth_before, float("nan"), "failed: time limit", n_c,
                                 time.time() - t0, rejected_energy_cycles=rejected)
        x, status = _solve_highs(c, constraints, lb, ub, n_v, remaining, mip_gap)
        if x is None:
            return GapfillResult([], growth_before, growth_before, float("nan"), f"failed: {status}", n_c, time.time() - t0,
                                 rejected_energy_cycles=rejected)
        y = x[n_v:]
        added = [cands[k].id for k in range(n_c) if y[k] > 0.5]
        res_fun = float(np.dot(c, x))
        # verify in cobra
        m2 = apply_gapfill(model, universe, added, forbid_o2_production=forbid_o2_production)
        m2.objective = m2.reactions.get_by_id(biomass_id)
        m2.objective_direction = "max"
        if forbid_energy_cycles and added:
            from .checks import energy_from_nothing
            egc = {k: v for k, v in energy_from_nothing(m2).items() if v > 1e-6}
            if egc:
                rejected.append(added)
                if verbose:
                    print(f"  gapfill: added set {added} creates an energy-generating cycle {egc}; cutting it off", flush=True)
                idx = [n_v + k for k in range(n_c) if y[k] > 0.5]
                cut = sp.csr_matrix(([1.0] * len(idx), ([0] * len(idx), idx)), shape=(1, n_v + n_c))
                constraints.append(LinearConstraint(cut, np.array([-np.inf]), np.array([len(idx) - 1.0])))
                continue
        apply_medium(m2, medium, close_all=True)
        for ex_id, lbv in (extra_uptakes or {}).items():
            if ex_id in m2.reactions:
                m2.reactions.get_by_id(ex_id).lower_bound = lbv
        growth_after = _growth(m2, allow_infeasible=True)
        if growth_after < min_growth - min(1e-7, min_growth * 1e-6):
            return GapfillResult([], growth_before, growth_after, float("nan"),
                                 "failed: added reactions do not meet min_growth on independent verification",
                                 n_c, time.time() - t0, rejected_energy_cycles=rejected)
        return GapfillResult(added, growth_before, growth_after, res_fun, status, n_c, time.time() - t0,
                             {a: weights[a] for a in added}, rejected_energy_cycles=rejected)
    return GapfillResult([], growth_before, growth_before, float("nan"), "failed: every solution creates an energy-generating cycle",
                         n_c, time.time() - t0, rejected_energy_cycles=rejected)


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
    x = np.asarray(sol.col_value, dtype=float) if len(sol.col_value) == n else None
    info = h.getInfo()
    if (x is None or not sol.value_valid
            or info.primal_solution_status != highspy.SolutionStatus.kSolutionStatusFeasible
            or not np.all(np.isfinite(x))):
        return None, status
    # A time-limited run can have a valid incumbent, but a vector returned by the
    # API alone is not a feasibility certificate. Check bounds, rows and binaries.
    tol = 1e-7
    ax = A @ x
    if (np.any(x < lb - tol) or np.any(x > ub + tol)
            or np.any(ax < rlo - tol) or np.any(ax > rhi + tol)
            or np.any(np.abs(x[n_v:] - np.rint(x[n_v:])) > tol)):
        return None, f"{status}: failed independent feasibility check"
    return x, status


def apply_gapfill(model: cobra.Model, universe: cobra.Model, added: List[str], note: str = "",
                  forbid_o2_production: bool = True) -> cobra.Model:
    """Copy the model and add reactions with the same oxygen policy as `gapfill`.

    Pass `forbid_o2_production=False` here too if that policy was disabled during
    gap-filling. Existing reaction capacities are retained in either case.
    """
    m = model.copy()
    rx = []
    for a in added:
        r = universe.reactions.get_by_id(a).copy()
        bounds = _candidate_bounds(r, forbid_o2_production)
        if bounds is None:
            raise ValueError(f"Gapfill reaction {a} has no allowed flux under the oxygen policy")
        r.bounds = bounds
        r.annotation = dict(r.annotation); r.annotation["gapfill"] = note or "added by gembench.gapfill"
        r.gene_reaction_rule = ""
        rx.append(r)
    m.add_reactions(rx)
    return m


def _is_transport(r: cobra.Reaction) -> bool:
    comps = {met.compartment for met in r.metabolites}
    return len(comps) > 1


def _candidate_bounds(r: cobra.Reaction, forbid_o2_production: bool):
    lower, upper = r.bounds
    if forbid_o2_production:
        # Count net molecular O2 across compartments. Transport moves oxygen but
        # does not produce it and must retain its allowed direction(s).
        net_o2 = sum(coef for met, coef in r.metabolites.items() if met.id.startswith("o2_"))
        if net_o2 > 0:
            upper = min(upper, 0.0)
        elif net_o2 < 0:
            lower = max(lower, 0.0)
    if lower > upper or (lower == 0.0 and upper == 0.0):
        return None
    return lower, upper


def _growth(model: cobra.Model, allow_infeasible: bool = False) -> float:
    value = model.slim_optimize()
    if allow_infeasible and model.solver.status == "infeasible":
        return 0.0
    if model.solver.status != "optimal" or value is None or not np.isfinite(value):
        raise RuntimeError(f"Growth optimization failed: solver status {model.solver.status}")
    return float(value)
