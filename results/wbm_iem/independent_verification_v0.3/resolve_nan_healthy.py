"""Independent re-solve of the healthy biomarker optima for which the MATLAB/Gurobi run printed NaN (claim 7 support).

Nothing from gembench is imported. The LP is assembled here from the model file and the MATLAB global bounds
(results/wbm_iem/matlab_reference/step2/matlab_global_bounds_Harvey_1_03d.mat), following checkIEM_WBM.m
(external/cobratoolbox/src/analysis/wholeBody/PSCMToolbox/checkIEM_WBM.m) rather than the Python port:
  * demand reactions DM_<met> for the metabolites the runIEM_HH.m block adds (bounds 0..1000, addDemandReaction default);
  * an auxiliary column a with a row  sum(IEM reaction fluxes) - a = 0,  -1e5 <= a <= 1e5;
  * vmax = max a;  healthy: a >= fix(vmax * 1e6) / 1e6 (truncation);
  * healthy biomarker value = max v_biomarker with its upper bound raised to 1e5.
IEM reactions: reactions whose id contains any of the block's patterns (runIEM_HH.m, strfind), minus 'Micro_' ones.
Each optimum is checked with my own code: row and bound violations of the returned vector against my own matrix,
the objective recomputed from the vector, and a weak-duality (Lagrangian) upper bound computed from the solver's row
duals:  for any y,  max c'x <= sum_i max(y_i*rl_i, y_i*ru_i) + sum_j max(z_j*l_j, z_j*u_j),  z = c - A'y.

Usage: nice -n 19 python3 resolve_nan_healthy.py EF ASNSD GMT PHOX1 SUCLA BTD
Appends one JSON object per IEM to resolve_nan_healthy.jsonl next to this script.
"""
import json
import math
import os
import sys
import time

import numpy as np
import scipy.io as sio
import scipy.sparse as sp

ROOT = "/home/claude/mma"
HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = os.path.join(ROOT, "external/COBRA.models/mat/Harvey_1_03d.mat")
GLOBAL = os.path.join(ROOT, "results/wbm_iem/matlab_reference/step2/matlab_global_bounds_Harvey_1_03d.mat")
OUT = os.path.join(HERE, "resolve_nan_healthy.jsonl")

# From runIEM_HH_ref.m (executed copy): patterns, demand metabolites, and the biomarkers whose MATLAB healthy value is NaN.
BLOCKS = {
    "EF": (["_HMR_8761", "_HMR_9800", "_KHK", "_KHK2", "_KHK3"], [], ["EX_fru[u]"]),
    "ASNSD": (["_ASNS1"], ["asn_L[bc]", "asn_L[csf]", "gln_L[bc]", "gln_L[csf]"], ["DM_asn_L[bc]", "DM_asn_L[csf]"]),
    "GMT": (["_GACMTRc"], ["creat[bc]", "crtn[bc]"], ["EX_urate[u]"]),
    "PHOX1": (["_AGTix", "_SPTix", "_r0160"], [], ["EX_oxa[u]", "EX_glyclt[u]", "EX_glx[u]"]),
    "SUCLA": (["_ITCOALm", "_MECOALm", "_SUCOASm", "_ITCOAL1m", "_MECOAS1m", "_SUCOAS1m"], ["lac_L[bc]", "pyr[bc]"],
              ["DM_lac_L[bc]", "DM_pyr[bc]"]),
    "BTD": (["_BTND1", "_BTND1n", "_BTNDe", "_BTNDm", "_ACCOACm", "_ACCOAC", "_PCm", "_MCCCrm", "_RE2453M", "_RE2454M",
             "_PPCOACm"], ["acac[bc]", "acetone[bc]", "bhb[bc]", "3ivcrn[bc]"],
            ["DM_acac[bc]", "DM_acetone[bc]", "DM_bhb[bc]", "EX_nh4[u]", "EX_3hpp[u]"]),
}


def cellstr(a):
    return [str(np.asarray(x).ravel()[0]) if np.asarray(x).size else "" for x in np.asarray(a, dtype=object).ravel()]


def load():
    d = sio.loadmat(MODEL, squeeze_me=False, struct_as_record=True)
    m = d["male"][0, 0]
    rxns, mets = cellstr(m["rxns"]), cellstr(m["mets"])
    S = sp.csc_matrix(m["S"], dtype=float)
    C = sp.csc_matrix(m["C"], dtype=float)
    b = np.asarray(m["b"], dtype=float).ravel()
    dvec = np.asarray(m["d"], dtype=float).ravel()
    csense = np.asarray(m["csense"]).astype(str).ravel()
    dsense = np.asarray(m["dsense"]).astype(str).ravel()
    g = sio.loadmat(GLOBAL, squeeze_me=False)
    assert cellstr(g["rxns"]) == rxns
    lb = np.asarray(g["lb"], dtype=float).ravel()
    ub = np.asarray(g["ub"], dtype=float).ravel()
    assert set(csense) == {"E"} and set(dsense) <= {"L", "G"}
    rl = np.concatenate([b, np.where(dsense == "G", dvec, -np.inf)])
    ru = np.concatenate([b, np.where(dsense == "L", dvec, np.inf)])
    return rxns, mets, S, C, rl, ru, lb, ub


def build(rxns, mets, S, C, rl0, ru0, lb0, ub0, patterns, demand_mets):
    n, ms = len(rxns), S.shape[0]
    iem = [j for j, r in enumerate(rxns) if any(p in r for p in patterns) and "Micro_" not in r]
    met_pos = {mm: i for i, mm in enumerate(mets)}
    names = list(rxns)
    cols_extra = []
    for mm in demand_mets:
        rid = "DM_" + mm
        assert rid not in names and mm in met_pos, rid
        cols_extra.append((rid, met_pos[mm]))
        names.append(rid)
    k = len(cols_extra)
    # demand columns: -1 in the metabolite row (S part only)
    D = sp.csc_matrix((-np.ones(k), ([i for _, i in cols_extra], np.arange(k))), shape=(ms, k)) if k else sp.csc_matrix((ms, 0))
    top = sp.hstack([S, D, sp.csc_matrix((ms, 1))])
    mid = sp.hstack([C, sp.csc_matrix((C.shape[0], k + 1))])
    aux_row = sp.csr_matrix((np.r_[np.ones(len(iem)), -1.0], (np.zeros(len(iem) + 1, dtype=int), np.r_[iem, n + k])),
                            shape=(1, n + k + 1))
    A = sp.vstack([top, mid, aux_row]).tocsc()
    rl = np.r_[rl0, 0.0]
    ru = np.r_[ru0, 0.0]
    lb = np.r_[lb0, np.zeros(k), -1e5]
    ub = np.r_[ub0, np.full(k, 1000.0), 1e5]
    names.append("aux_IEM_sum")
    return A, rl, ru, lb, ub, names, [rxns[j] for j in iem]


def solve(A, rl, ru, lb, ub, obj_col):
    import highspy
    n = A.shape[1]
    lp = highspy.HighsLp()
    lp.num_col_, lp.num_row_ = n, A.shape[0]
    c = np.zeros(n)
    c[obj_col] = 1.0
    lp.col_cost_ = c
    lp.col_lower_, lp.col_upper_ = lb.copy(), ub.copy()
    lp.row_lower_, lp.row_upper_ = rl.copy(), ru.copy()
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A.indptr.astype(np.int32)
    lp.a_matrix_.index_ = A.indices.astype(np.int32)
    lp.a_matrix_.value_ = A.data.astype(float)
    lp.sense_ = highspy.ObjSense.kMaximize
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("solver", "ipm")
    h.setOptionValue("run_crossover", "on")
    h.setOptionValue("primal_feasibility_tolerance", 1e-7)
    h.setOptionValue("dual_feasibility_tolerance", 1e-7)
    h.setOptionValue("threads", 1)
    h.setOptionValue("time_limit", 3600.0)
    h.passModel(lp)
    t = time.time()
    h.run()
    dt = time.time() - t
    status = h.modelStatusToString(h.getModelStatus())
    sol = h.getSolution()
    x = np.array(sol.col_value) if sol.value_valid else None
    y = np.array(sol.row_dual) if sol.dual_valid else None
    return status, x, y, c, dt


def certificate(A, rl, ru, lb, ub, c, x, y):
    r = A @ x
    row_v = np.maximum(np.maximum(rl - r, r - ru), 0.0)
    col_v = np.maximum(np.maximum(lb - x, x - ub), 0.0)
    out = {"max_row_violation": float(row_v.max()), "n_rows_over_1e-6": int((row_v > 1e-6).sum()),
           "max_bound_violation": float(col_v.max()), "n_bounds_over_1e-6": int((col_v > 1e-6).sum()),
           "objective_from_vector": float(c @ x)}
    if y is not None:
        # Weak duality holds for ANY y, so components of y that would multiply an infinite row bound are set to 0
        # (projection); the bound stays valid. All column bounds here are finite (|bound| <= 1e6).
        best = None
        for sign in (1.0, -1.0):
            yy = sign * y
            bad = ((yy > 0) & ~np.isfinite(ru)) | ((yy < 0) & ~np.isfinite(rl))
            dropped = float(np.abs(yy[bad]).sum())
            yy = np.where(bad, 0.0, yy)
            z = c - A.T @ yy

            def term(coef, lo, hi):
                v = np.zeros_like(coef)
                pos, neg = coef > 0, coef < 0
                v[pos] = coef[pos] * hi[pos]
                v[neg] = coef[neg] * lo[neg]
                return v
            total = float(term(yy, rl, ru).sum() + term(z, lb, ub).sum())
            if np.isfinite(total) and (best is None or total < best[0]):
                best = (total, sign, int(bad.sum()), dropped)
        if best is not None:
            out["lagrangian_upper_bound"] = best[0]
            out["dual_sign_used"] = best[1]
            out["n_dual_components_projected"] = best[2]
            out["sum_abs_projected_duals"] = best[3]
            out["duality_gap"] = best[0] - out["objective_from_vector"]
    return out


def main():
    iems = sys.argv[1:] or list(BLOCKS)
    rxns, mets, S, C, rl0, ru0, lb0, ub0 = load()
    v03 = {r["iem"]: r for r in json.load(open(os.path.join(ROOT, "results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json")))}
    for iem in iems:
        patterns, demands, targets = BLOCKS[iem]
        A, rl, ru, lb, ub, names, iem_rxns = build(rxns, mets, S, C, rl0, ru0, lb0, ub0, patterns, demands)
        aux = len(names) - 1
        rec = {"iem": iem, "n_iem_rxns": len(iem_rxns),
               "iem_rxns_equal_python_record": sorted(iem_rxns) == sorted(v03[iem]["iem_reactions"])}
        st, x, y, c, dt = solve(A, rl, ru, lb, ub, aux)
        vmax = float(x[aux]) if x is not None else float("nan")
        rec["vmax"] = {"status": st, "value": vmax, "python_v0.3": v03[iem]["vmax_healthy"], "seconds": round(dt, 1),
                       **(certificate(A, rl, ru, lb, ub, c, x, y) if x is not None else {})}
        print(f"{iem}: vmax {st} {vmax!r} (python {v03[iem]['vmax_healthy']!r}) {dt:.0f}s", flush=True)
        lo = math.trunc(vmax * 1e6) / 1e6        # MATLAB fix()
        lb_h = lb.copy()
        lb_h[aux] = lo
        rec["healthy_aux_lb"] = lo
        rec["biomarkers"] = []
        py = {b["reaction"]: b for b in v03[iem]["biomarkers"]}
        for t in targets:
            j = names.index(t)
            ub_t = ub.copy()
            ub_t[j] = 1e5
            st, x, y, c, dt = solve(A, rl, ru, lb_h, ub_t, j)
            cert = certificate(A, rl, ru, lb_h, ub_t, c, x, y) if x is not None else {}
            val = cert.get("objective_from_vector", float("nan"))
            pv = py[t]["healthy"]
            diff = abs(val - pv) if x is not None and pv is not None else None
            rec["biomarkers"].append({"biomarker": t, "status": st, "value": val, "python_v0.3_healthy": pv,
                                      "abs_diff_vs_python": diff, "seconds": round(dt, 1), **cert})
            print(f"   {t}: {st} value {val!r} python {pv!r} |diff| {diff} row {cert.get('max_row_violation')} "
                  f"bound {cert.get('max_bound_violation')} UB {cert.get('lagrangian_upper_bound')} {dt:.0f}s", flush=True)
        with open(OUT, "a") as fh:
            fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
