"""Cross-disease readouts (scripts/run_wbm_iem_cross.py) on E. coli core dressed up as a whole-body model.

Checks that own-biomarker readouts reproduce the unchanged protocol code (run_iem), that simplex warm starts
give the same optima and calls as interior point, that the model is left as found, and that the Gurobi
backend agrees with HiGHS (skipped without gurobipy; the size-limited licence is enough for this model).
"""
import copy
import importlib.util
import json
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


def test_rechecks_are_recorded_and_agree(core_wbm):
    hw = backend("highs", core_wbm)
    rec = X.run_one(hw, PROTOCOL[1], READOUTS, "protocol", "primal", log=False, recheck_every=1)
    checks = [r[f"solve_{s}"]["recheck"] for r in rec["readouts"] for s in ("healthy", "disease")
              if "recheck" in (r[f"solve_{s}"] or {})]
    # every warm solve (all but the first present readout of each state) is rechecked when N = 1
    n_present = sum(1 for r in rec["readouts"] if r["status_healthy"] != "absent")
    assert len(checks) == 2 * (n_present - 1)
    assert all(c["status"] == "Optimal" and c["abs_diff"] <= 1e-6 for c in checks)
    # a recheck is a real solve from scratch, not the stored solution handed back
    assert all((c["ipm_iterations"] or 0) > 0 for c in checks)


def test_rechecked_selection_is_deterministic():
    picks = [X.rechecked("HIS", "healthy", f"R{i}", 10) for i in range(1000)]
    assert picks == [X.rechecked("HIS", "healthy", f"R{i}", 10) for i in range(1000)]
    assert 60 <= sum(picks) <= 140
    assert not X.rechecked("HIS", "healthy", "R1", 0)


@pytest.mark.parametrize("context", ["protocol", "minimal"])
@pytest.mark.parametrize("name,warm,first", [("highs", "dual", "ipm"), ("highs", "primal", "warm"), ("highs", "ipm", "ipm"),
                                             ("gurobi", "dual", "ipm"), ("gurobi", "dual", "warm"),
                                             ("gurobi", "concurrent", "ipm")])
def test_readout_major_equals_iem_major(core_wbm, context, name, warm, first):
    base = cross(backend("highs", core_wbm), "ipm", context)
    hw = backend(name, core_wbm)
    lb0, ub0 = hw.lb.copy(), hw.ub.copy()
    setups, entries = X.run_readout_major(hw, PROTOCOL, READOUTS, context, warm, first, False, 0, {}, lambda s, e: None)
    recs = {r["iem"]: r for r in X.assemble(setups, entries, READOUTS, final=True)}
    np.testing.assert_array_equal(hw.lb, lb0)
    np.testing.assert_array_equal(hw.ub, ub0)
    for iem, ref in base.items():
        rec = recs[iem]
        assert rec["status"] == ref["status"]
        assert rec["vmax_healthy"] == pytest.approx(ref["vmax_healthy"])
        assert [r["reaction"] for r in rec["readouts"]] == [r["reaction"] for r in ref["readouts"]]
        for ra, rb in zip(ref["readouts"], rec["readouts"]):
            assert ra["predicted"] == rb["predicted"], (iem, ra["reaction"])
            assert same(ra["healthy"], rb["healthy"]) and same(ra["disease"], rb["disease"]), (iem, ra["reaction"])
            assert ra["own"] == rb["own"] and ra["expected"] == rb["expected"]


def test_readout_major_resumes(core_wbm):
    hw = backend("highs", core_wbm)
    setups, entries = X.run_readout_major(hw, PROTOCOL, READOUTS, "protocol", "dual", "ipm", False, 0, {}, lambda s, e: None)
    first = READOUTS[0]
    done = {k: v for k, v in entries.items() if k[1] == first}
    hw2 = backend("highs", core_wbm)
    calls = []
    _, entries2 = X.run_readout_major(hw2, PROTOCOL, READOUTS, "protocol", "dual", "ipm", False, 0, done,
                                      lambda s, e: calls.append(len(e)))
    assert len(calls) == len(READOUTS) - 1     # the stored readout is not recomputed
    for k, v in entries.items():
        assert entries2[k]["predicted"] == v["predicted"]


def test_shards_merge_to_the_full_run(core_wbm):
    import hashlib
    import json as _json
    spec_m = importlib.util.spec_from_file_location("merge_cross_shards", os.path.join(ROOT, "scripts", "merge_cross_shards.py"))
    M = importlib.util.module_from_spec(spec_m)
    spec_m.loader.exec_module(M)
    full_setups, full_entries = X.run_readout_major(backend("highs", core_wbm), PROTOCOL, READOUTS, "protocol", "dual", "ipm",
                                                    False, 0, {}, lambda s, e: None)
    full = {r["iem"]: r for r in X.assemble(full_setups, full_entries, READOUTS, final=True)}
    n = 3
    panel_sha = hashlib.sha256(_json.dumps(READOUTS).encode()).hexdigest()
    shards = []
    for k in range(1, n + 1):
        sub = READOUTS[k - 1::n]
        setups, entries = X.run_readout_major(backend("highs", core_wbm), PROTOCOL, sub, "protocol", "dual", "ipm",
                                              False, 0, {}, lambda s, e: None)
        recs = X.assemble(setups, entries, sub, final=True)
        prov = {"engine_version": "test", "shard": f"{k}/{n}", "panel_sha256": panel_sha, "n_panel": len(READOUTS),
                "readouts_sha256": hashlib.sha256(_json.dumps(sub).encode()).hexdigest(), "n_readouts": len(sub)}
        recs = _json.loads(_json.dumps(X.json_safe([dict(r, provenance=prov) for r in recs])))
        shards.append((f"shard{k}", recs))
    merged, panel = M.merge(shards)
    assert panel == READOUTS
    for rec in merged:
        ref = full[rec["iem"]]
        assert rec["status"] == ref["status"]
        assert [e["reaction"] for e in rec["readouts"]] == READOUTS
        for a, b in zip(ref["readouts"], rec["readouts"]):
            assert a["predicted"] == b["predicted"]
            assert same(a["healthy"], b["healthy"] if b["healthy"] is not None else float("nan"))
    with pytest.raises(ValueError):
        M.merge(shards[:2])          # a missing shard is refused


@pytest.mark.parametrize("name", ["highs", "gurobi"])
def test_readout_major_rechecks_really_resolve(core_wbm, name):
    hw = backend(name, core_wbm)
    setups, entries = X.run_readout_major(hw, PROTOCOL, READOUTS, "protocol", "dual", "ipm", False, 1, {}, lambda s, e: None)
    checks = [e[f"solve_{st}"]["recheck"] for e in entries.values() for st in ("healthy", "disease")
              if "recheck" in (e[f"solve_{st}"] or {})]
    assert checks and all(c["status"] == "Optimal" and (c["ipm_iterations"] or 0) > 0 and c["abs_diff"] <= 1e-6 for c in checks)


def test_prepare_from_stored_record_rebuilds_the_same_lps(core_wbm):
    hw = backend("highs", core_wbm)
    setups, entries = X.run_readout_major(hw, PROTOCOL, READOUTS, "protocol", "ipm", "ipm", False, 0, {}, lambda s, e: None)
    recs = {r["iem"]: r for r in X.assemble(setups, entries, READOUTS, final=True)}
    hw2 = backend("highs", core_wbm)
    rebuilt = [X.prepare(hw2, p, "protocol", log=False, stored=recs[p["iem"]]) for p in PROTOCOL]
    assert hw2.n_solves == 0
    mover = X.Mover(hw2)
    for s in rebuilt:
        for rid in READOUTS:
            if rid not in hw2.rxn_pos:
                continue
            rcol = hw2.rxn_pos[rid]
            hw2.set_objective({rcol: 1.0}, "max")
            for state in ("healthy", "disease"):
                target, rb = X.target_bounds(s, state, rcol, mover.base_lb, mover.base_ub)
                mover.move(target, s.row, rb)
                hw2.fresh(); hw2.method = "ipm"
                st, f, _, _ = hw2.solve()
                f = 0.0 if abs(f) <= 1e-6 else f
                assert same(f, entries[(s.iem, rid)][state]), (s.iem, rid, state)


def test_recheck_retries_an_infeasible_healthy_pin_with_its_own(core_wbm):
    """scripts/recheck_cross_matrix.py: a pin just above this solver's maximum makes the healthy state infeasible;
    the retry with the pin from this solver's own v_max reproduces the matrix value and leaves the pin as it was."""
    spec_r = importlib.util.spec_from_file_location("recheck_cross_matrix", os.path.join(ROOT, "scripts", "recheck_cross_matrix.py"))
    R = importlib.util.module_from_spec(spec_r)
    spec_r.loader.exec_module(R)
    hw = backend("highs", core_wbm)
    setups, entries = X.run_readout_major(hw, PROTOCOL, READOUTS, "protocol", "ipm", "ipm", False, 0, {}, lambda s, e: None)
    recs = {r["iem"]: r for r in X.assemble(setups, entries, READOUTS, final=True)}
    hw2 = backend("highs", core_wbm)
    rebuilt = {p["iem"]: X.prepare(hw2, p, "protocol", log=False, stored=recs[p["iem"]]) for p in PROTOCOL}
    mover = X.Mover(hw2)
    s = rebuilt["TPI_def"]
    vmax = recs["TPI_def"]["vmax_healthy"]
    assert R.truncated_pin(vmax) == recs["TPI_def"]["healthy_pin"]
    rid = "DM_dhap_c"
    rcol = hw2.rxn_pos[rid]
    hw2.set_objective({rcol: 1.0}, "max")
    s.lo = vmax + 1e-3                                   # a pin from "another solver", above this one's maximum
    status, f, _, _, first = R.solve_one(hw2, mover, s, "healthy", rcol)
    assert status == "Infeasible" and first is None      # no retry without a pin to retry with
    status, f, _, _, first = R.solve_one(hw2, mover, s, "healthy", rcol, alt_lo=R.truncated_pin(vmax))
    assert status == "Optimal" and first["status"] == "Infeasible" and first["pin"] == vmax + 1e-3
    assert first["retry_pin"] == R.truncated_pin(vmax) and s.lo == vmax + 1e-3
    assert same(0.0 if abs(f) <= 1e-6 else f, entries[("TPI_def", rid)]["healthy"])
    status, f, _, _, first = R.solve_one(hw2, mover, s, "disease", rcol, alt_lo=R.truncated_pin(vmax))
    assert status == "Optimal" and first is None         # the disease state does not use the pin
    assert same(0.0 if abs(f) <= 1e-6 else f, entries[("TPI_def", rid)]["disease"])


def _pin_sweep_module():
    spec_s = importlib.util.spec_from_file_location("run_wbm_iem_pin_sweep", os.path.join(ROOT, "scripts", "run_wbm_iem_pin_sweep.py"))
    S = importlib.util.module_from_spec(spec_s)
    spec_s.loader.exec_module(S)
    return S


@pytest.mark.parametrize("name", ["highs", "gurobi"])
def test_pin_sweep_reproduces_the_matrix_at_alpha_1_and_relaxes_monotonically(core_wbm, name):
    """scripts/run_wbm_iem_pin_sweep.py: alpha = 1 rebuilds the matrix's healthy LPs; a lower pin can only raise a
    healthy maximum; every alpha level equals a from-scratch barrier solve of the same LP."""
    S = _pin_sweep_module()
    hw = backend("highs", core_wbm)
    setups, entries = X.run_readout_major(hw, PROTOCOL, READOUTS, "protocol", "ipm", "ipm", False, 0, {}, lambda s, e: None)
    recs = {r["iem"]: r for r in X.assemble(setups, entries, READOUTS, final=True)}
    hw2 = backend(name, core_wbm)
    rebuilt = [X.prepare(hw2, p, "protocol", log=False, stored=recs[p["iem"]]) for p in PROTOCOL]
    rebuilt = [s for s in rebuilt if s.ok]
    mover = X.Mover(hw2)
    alphas = [1.0, 0.5, 0.1]
    for rid in READOUTS:
        out = S.sweep_readout(hw2, mover, rebuilt, rid, alphas, "dual", "ipm", 2)
        for s in rebuilt:
            got = out[s.iem]
            assert [e["alpha"] for e in got] == alphas
            if rid not in hw2.rxn_pos:
                assert all(e["status"] == "absent" for e in got)
                continue
            assert all(e["status"] == "Optimal" for e in got), (s.iem, rid, got)
            assert s.lo == recs[s.iem]["healthy_pin"]                      # the pin is restored
            assert same(got[0]["value"], entries[(s.iem, rid)]["healthy"]), (s.iem, rid)
            vals = [e["value"] for e in got]
            assert all(b >= a - 2e-6 for a, b in zip(vals, vals[1:])), (s.iem, rid, vals)
            assert got[1]["pin"] == math.floor(0.5 * recs[s.iem]["vmax_healthy"] * 1e6) / 1e6
            for e in got:
                if "recheck" in e:
                    assert e["recheck"]["status"] == "Optimal" and e["recheck"]["abs_diff"] <= 2e-6
    # Each level against an independent from-scratch solve of the same LP.
    hw3 = backend("highs", core_wbm)
    fresh = [X.prepare(hw3, p, "protocol", log=False, stored=recs[p["iem"]]) for p in PROTOCOL]
    fresh = [s for s in fresh if s.ok]
    mover3 = X.Mover(hw3)
    rid = "EX_ac_e"
    out = S.sweep_readout(hw2, mover, rebuilt, rid, alphas, "dual", "ipm", 0)
    rcol = hw3.rxn_pos[rid]
    hw3.set_objective({rcol: 1.0}, "max")
    for s in fresh:
        for e in out[s.iem]:
            s.lo = e["pin"]
            target, rb = X.target_bounds(s, "healthy", rcol, mover3.base_lb, mover3.base_ub)
            mover3.move(target, s.row, rb)
            hw3.fresh(); hw3.method = "ipm"
            st, f, _, _ = hw3.solve()
            assert st == "Optimal" and same(0.0 if abs(f) <= 1e-6 else f, e["value"]), (s.iem, e)


def test_pin_sweep_merge_checks_completeness(tmp_path):
    S = _pin_sweep_module()
    prov = {"n_panel": 3, "alphas": [1.0], "threads": 1}
    a = {"provenance": dict(prov, shard="1/2", n_readouts=2, readouts_sha256="x"), "setups": {"A": {}},
         "readouts": {"r1": {"A": []}, "r3": {"A": []}}}
    b = {"provenance": dict(prov, shard="2/2", n_readouts=1, readouts_sha256="y"), "setups": {"A": {}},
         "readouts": {"r2": {"A": []}}}
    pa, pb = tmp_path / "a.json", tmp_path / "b.json"
    pa.write_text(json.dumps(a)); pb.write_text(json.dumps(b))
    merged = S.merge([str(pa), str(pb)], str(tmp_path / "m.json"))
    assert sorted(merged["readouts"]) == ["r1", "r2", "r3"] and merged["provenance"]["shards"] == 2
    b["readouts"] = {}
    pb.write_text(json.dumps(b))
    with pytest.raises(ValueError):
        S.merge([str(pa), str(pb)], str(tmp_path / "m.json"))
    b["readouts"] = {"r2": {"A": []}}
    b["provenance"]["alphas"] = [0.5]
    pb.write_text(json.dumps(b))
    with pytest.raises(ValueError):
        S.merge([str(pa), str(pb)], str(tmp_path / "m.json"))


def test_pin_for_truncates_like_the_protocol():
    S = _pin_sweep_module()
    assert S.pin_for(1234.5678919, 1.0, 1234.567891) == 1234.567891     # alpha 1: the matrix's pin
    assert S.pin_for(1234.5678919, 0.5) == math.floor(617.28394595 * 1e6) / 1e6
    assert S.pin_for(0.0042433602, 0.001) == 0.000004
    assert S.pin_for(-10.0000017, 0.5) == math.ceil(-5.00000085 * 1e6) / 1e6
