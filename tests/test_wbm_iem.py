"""Small solvable networks test numerical semantics without downloading Harvey."""
import numpy as np
import pytest
import scipy.sparse as sp
import subprocess
import sys

from gembench.wbm import WBM
from gembench.wbm_iem import HighsWBM, run_iem


def toy_model(capacity=10.0, product_outlet=True):
    # supply -> A --IEM--> B, with separate A and B outlets.
    return WBM(
        name="toy", rxns=np.array(["supply", "organ_IEM", "EX_a", "EX_b", "Whole_body_objective_rxn"]),
        mets=np.array(["a[bc]", "b[bc]"]),
        S=sp.csc_matrix([[1, -1, -1, 0, 0], [0, 1, 0, -1, 0]]),
        b=np.zeros(2), csense=np.array(["E", "E"]), C=sp.csc_matrix((0, 5)),
        d=np.array([]), dsense=np.array([]), ctrs=np.array([]),
        lb=np.array([0., 0., 0., 0., 1.]),
        ub=np.array([capacity, capacity, capacity, capacity if product_outlet else 0., 1.]),
        c=np.array([0., 0., 0., 0., 1.]), osense="max",
    )


def run(hw, biomarkers=(("EX_a", "Increased"),)):
    return run_iem(hw, "toy", ["_IEM"], [], biomarkers, verbose=False)


def test_knockout_prediction_and_original_bounds_restored():
    model = toy_model()
    hw = HighsWBM(model)
    result = run(hw, [("EX_a", "Increased"), ("EX_b", "Decreased")])
    assert result.vmax_healthy == pytest.approx(10)
    assert result.wb_objective_disease_feasible
    assert [b.predicted for b in result.biomarkers] == ["Increased", "Decreased"]
    assert all(b.correct for b in result.biomarkers)
    assert result.n_solves == 7
    assert result.status == "complete"
    np.testing.assert_array_equal(hw.lb, model.lb)
    np.testing.assert_array_equal(hw.ub, model.ub)


def test_joint_objective_has_toolbox_auxiliary_variable_cap():
    result = run(HighsWBM(toy_model(capacity=200000)))
    assert result.vmax_healthy == pytest.approx(100000)


@pytest.mark.parametrize("status,value", [("Time limit reached", 10.), ("Infeasible", np.nan), ("Optimal", np.nan)])
def test_unavailable_optimum_is_not_a_zero_flux_prediction(monkeypatch, status, value):
    hw = HighsWBM(toy_model())
    solve = hw.solve

    def interrupted_solve():
        answer = solve()
        if hw.n_solves == 4:  # healthy biomarker solve
            return status, value, answer[2], answer[3]
        return answer

    monkeypatch.setattr(hw, "solve", interrupted_solve)
    result = run(hw)
    bm = result.biomarkers[0]
    assert np.isnan(bm.healthy)
    assert bm.predicted == "NA"
    assert bm.correct is None
    assert result.notes
    assert result.status == "partial"


def test_added_demands_do_not_leak_between_iems_and_can_be_reused():
    hw = HighsWBM(toy_model(product_outlet=False))
    first = run(hw, [("DM_b[bc]", "Decreased")])
    assert first.vmax_healthy == pytest.approx(10)
    assert hw.ub[hw.rxn_pos["DM_b[bc]"]] == 0
    without_demand = run(hw)
    assert without_demand.vmax_healthy == pytest.approx(0)
    assert not without_demand.biomarkers
    again = run(hw, [("DM_b[bc]", "Decreased")])
    assert again.vmax_healthy == pytest.approx(first.vmax_healthy)
    assert again.biomarkers[0].healthy == pytest.approx(first.biomarkers[0].healthy)
    assert hw.ub[hw.rxn_pos["DM_b[bc]"]] == 0


def test_exception_restores_knockout_bounds_objective_and_demand_sinks(monkeypatch):
    model = toy_model(product_outlet=False)
    hw = HighsWBM(model)
    hw.set_objective({2: 1.0})
    solve = hw.solve

    def failing_solve():
        if hw.n_solves == 4:
            raise RuntimeError("simulated solver failure")
        return solve()

    monkeypatch.setattr(hw, "solve", failing_solve)
    with pytest.raises(RuntimeError, match="simulated solver failure"):
        run(hw, [("DM_b[bc]", "Decreased")])
    np.testing.assert_array_equal(hw.lb[:model.n_rxns], model.lb)
    np.testing.assert_array_equal(hw.ub[:model.n_rxns], model.ub)
    assert hw.ub[-1] == 0
    # The original objective remains solvable; all temporary rows are inert.
    assert solve()[1] == pytest.approx(10)


def test_failed_initial_solve_still_reports_work_and_restores_model(monkeypatch):
    hw = HighsWBM(toy_model())
    solve = hw.solve

    def timeout():
        st, obj, x, dt = solve()
        return "Time limit reached", obj, x, dt

    monkeypatch.setattr(hw, "solve", timeout)
    result = run(hw)
    assert result.n_solves == 1
    assert not result.biomarkers
    assert hw._objective == {}


def test_source_demands_not_used_as_biomarkers_are_still_applied_temporarily():
    hw = HighsWBM(toy_model(product_outlet=False))
    result = run_iem(hw, "toy", ["_IEM"], [], [("EX_a", "Increased")],
                     verbose=False, demand_metabolites=["b[bc]"])
    assert result.vmax_healthy == pytest.approx(10)
    assert result.biomarkers[0].correct
    assert hw.ub[hw.rxn_pos["DM_b[bc]"]] == 0


def test_existing_global_scheduler_is_reused_and_explicit_conflict_is_reported():
    # Isolate scheduler setup so the test does not reset another solver's state.
    script = f"""
import highspy
import runpy
from gembench.wbm_iem import HighsWBM
toy_model = runpy.run_path({__file__!r})['toy_model']
previous = highspy.Highs()
previous.setOptionValue('output_flag', False)
previous.setOptionValue('threads', 2)
previous.addVar(0, 1)
assert previous.run() == highspy.HighsStatus.kOk
hw = HighsWBM(toy_model())
assert hw.solve()[0] == 'Optimal'
conflicting = HighsWBM(toy_model(), threads=1)
try:
    conflicting.solve()
except RuntimeError as error:
    assert 'scheduler' in str(error)
else:
    raise AssertionError('An explicit incompatible thread count must be reported')
assert previous.run() == highspy.HighsStatus.kOk
"""
    completed = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, timeout=30)
    assert completed.returncode == 0, completed.stdout + completed.stderr
