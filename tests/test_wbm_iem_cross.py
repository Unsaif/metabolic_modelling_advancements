"""Cross-disease readouts (scripts/run_wbm_iem_cross.py) on E. coli core dressed up as a whole-body model.

Checks that own-biomarker readouts reproduce the unchanged protocol code (run_iem), that simplex warm starts
give the same optima and calls as interior point, that the model is left as found, and that the Gurobi
backend agrees with HiGHS (skipped without gurobipy; the size-limited licence is enough for this model).
"""
import copy
import importlib.util
import math
import os

import numpy as np
import pytest
import scipy.sparse as sp

from gembench import wbm_iem as I
from gembench.wbm import WBM

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("run_wbm_iem_cross", os.path.join(ROOT, "scripts", "run_wbm_iem_cross.py"))
X = importlib.util.module_from_spec(spec)
spec.loader.exec_module(X)

PROTOCOL = [
    {"iem": "TPI_def", "call_index": 1, "include_patterns": ["TPI"], "exclude_patterns": [], "bound_tweaks": [],
     "biomarkers": [["DM_dhap_c", "Increased (blood)"], ["EX_ac_e", "Decreased (urine)"], ["EX_lac__D_e", "Increased"]],
     "demand_metabolites": ["g3p_c"]},
    {"iem": "PGI_def", "call_index": 2, "include_patterns": ["PGI"], "exclude_patterns": [],
     "bound_tweaks": [{"pattern": "G6PDH2r", "bound": "lb", "value": 0.0}],
     "biomarkers": [["DM_g6p_c", "Increased"], ["EX_etoh_e", "Increased"], ["DM_f6p_c", "Decreased"]],
     "demand_metabolites": []},
    {"iem": "FUM_def", "call_index": 3, "include_patterns": ["FUM"], "exclude_patterns": [], "bound_tweaks": [],
     "biomarkers": [["DM_fum_c", "Increased"], ["EX_succ_e", "Increased"], ["DM_nothere_c", "Increased"]],
     "demand_metabolites": ["mal__L_c"]},
]
READOUTS = list(dict.fromkeys(r for p in PROTOCOL for r, _ in p["biomarkers"]))


@pytest.fixture(scope="module")
def core_wbm():
    cobra = pytest.importorskip("cobra")
    from cobra.util.array import create_stoichiometric_matrix
    cm = cobra.io.load_model("textbook")
    bio = cm.reactions.get_by_id("Biomass_Ecoli_core")
    bio.id = "Whole_body_objective_rxn"
    cm.repair()
    bio.bounds = (0.5, 0.5)
    rxns = np.array([r.id for r in cm.reactions]); mets = np.array([m.id for m in cm.metabolites])
    n = len(rxns)
    C = sp.lil_matrix((1, n)); C[0, list(rxns).index("PGI")] = 1.0; C[0, list(rxns).index("PFK")] = -20.0
    return WBM(name="core", rxns=rxns, mets=mets, S=sp.csc_matrix(create_stoichiometric_matrix(cm, array_type="dok")),
               b=np.zeros(len(mets)), csense=np.array(["E"] * len(mets)), C=sp.csc_matrix(C), d=np.array([5.0]),
               dsense=np.array(["L"]), ctrs=np.array(["c1"]),
               lb=np.array([r.lower_bound for r in cm.reactions], float),
               ub=np.array([r.upper_bound for r in cm.reactions], float), c=np.zeros(n), osense="max")


def backend(name, wbm):
    m = copy.deepcopy(wbm)
    if name == "gurobi":
        G = pytest.importorskip("gembench.wbm_iem_gurobi")
        pytest.importorskip("gurobipy")
        hw = G.GurobiWBM(m)
    else:
        hw = I.HighsWBM(m)
    sinks = list(dict.fromkeys([r[3:] for r in READOUTS if r.startswith("DM_")] +
                               [x for p in PROTOCOL for x in p.get("demand_metabolites") or []]))
    for met in sinks:
        if met in hw.met_pos:
            hw.add_demand(met, lb=0.0, ub=0.0)
    return hw


def cross(hw, warm, context="protocol"):
    return {p["iem"]: X.run_one(hw, p, READOUTS, context, warm, log=False) for p in PROTOCOL}


def same(a, b, tol=2e-6):
    if a is None or b is None:
        return a is None and b is None
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    return abs(a - b) <= tol


def test_own_readouts_reproduce_protocol_and_state_is_restored(core_wbm):
    ref_hw = backend("highs", core_wbm)
    ref = {}
    for p in PROTOCOL:
        with ref_hw.temporary_state():
            for tw in p["bound_tweaks"]:
                idx = I.match_reactions(ref_hw.wbm.rxns, [tw["pattern"]])
                ref_hw.set_bounds(idx, **{tw["bound"]: [tw["value"]] * len(idx)})
            ref[p["iem"]] = I.run_iem(ref_hw, p["iem"], p["include_patterns"], p["exclude_patterns"],
                                      [tuple(b) for b in p["biomarkers"]], demand_metabolites=p["demand_metabolites"],
                                      verbose=False)
    hw = backend("highs", core_wbm)
    lb0, ub0 = hw.lb.copy(), hw.ub.copy()
    out = cross(hw, "primal")
    np.testing.assert_array_equal(hw.lb, lb0)
    np.testing.assert_array_equal(hw.ub, ub0)
    for iem, rec in out.items():
        assert rec["status"] == "complete"
        assert rec["vmax_healthy"] == pytest.approx(ref[iem].vmax_healthy)
        by = {r["reaction"]: r for r in rec["readouts"]}
        for b in ref[iem].biomarkers:
            r = by[b.reaction]
            assert r["own"] and r["expected"] == b.expected
            assert r["predicted"] == b.predicted
            assert same(r["healthy"], b.healthy) and same(r["disease"], b.disease)
        absent = by["DM_nothere_c"]
        assert absent["status_healthy"] == "absent" and absent["predicted"] == "NA"


@pytest.mark.parametrize("context", ["protocol", "minimal"])
@pytest.mark.parametrize("name,warm", [("highs", "primal"), ("highs", "dual"), ("gurobi", "ipm"), ("gurobi", "primal")])
def test_backends_and_warm_starts_agree(core_wbm, context, name, warm):
    base = cross(backend("highs", core_wbm), "ipm", context)
    out = cross(backend(name, core_wbm), warm, context)
    for iem in base:
        assert out[iem]["status"] == base[iem]["status"]
        for ra, rb in zip(base[iem]["readouts"], out[iem]["readouts"]):
            assert ra["reaction"] == rb["reaction"]
            assert ra["predicted"] == rb["predicted"], (iem, ra["reaction"])
            assert same(ra["healthy"], rb["healthy"]) and same(ra["disease"], rb["disease"]), (iem, ra["reaction"])


def test_minimal_context_keeps_own_sinks_closed(core_wbm):
    hw = backend("highs", core_wbm)
    rec = X.run_one(hw, PROTOCOL[0], READOUTS, "minimal", "ipm", log=False)
    # With its g3p sink closed, the TPI model's v_max differs from the protocol context only if the sink matters;
    # either way the run completes and the model is left as found.
    assert rec["status"] == "complete"
    j = hw.rxn_pos["DM_g3p_c"]
    assert hw.ub[j] == 0.0
