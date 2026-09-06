import json
import sys

import numpy as np
import pytest
import scipy.io as sio

from gembench.wbm_iem import HighsWBM
from scripts import run_wbm_iem as runner
from test_wbm_iem import toy_model


def test_resume_rejects_unfingerprinted_and_changed_experiments(tmp_path):
    path = tmp_path / "results.json"
    path.write_text('[{"iem": "CPS1"}]')
    with pytest.raises(ValueError, match="Choose a new"):
        runner.load_resume(path, "new")
    path.write_text('[{"run_fingerprint": "old"}]')
    with pytest.raises(ValueError, match="different"):
        runner.load_resume(path, "new")
    assert runner.load_resume(path, "old") == [{"run_fingerprint": "old"}]


def test_stored_timeout_with_correctness_flag_is_not_scored():
    bad = {"correct": True, "status_healthy": "Time limit reached", "status_disease": "Optimal", "healthy": 0., "disease": 1.}
    assert runner.scored_biomarkers([{"biomarkers": [bad]}]) == []


def test_runner_retries_unknown_results_and_writes_strict_json(tmp_path, monkeypatch):
    model = toy_model()
    model_file = tmp_path / "toy.mat"
    fields = {name: getattr(model, name) for name in ("S", "b", "csense", "C", "d", "dsense", "lb", "ub", "c")}
    fields.update(rxns=np.asarray(model.rxns, dtype=object)[:, None],
                  mets=np.asarray(model.mets, dtype=object)[:, None],
                  ctrs=np.empty((0, 1), dtype=object))
    sio.savemat(model_file, {"toy": fields})
    protocol = tmp_path / "protocol.json"
    protocol.write_text(json.dumps([{"iem": "toy", "call_index": 1, "include_patterns": ["_IEM"],
                                     "exclude_patterns": [], "bound_tweaks": [], "biomarkers": [["EX_a", "Increased"]]}]))
    monkeypatch.setattr(runner, "OUT", str(tmp_path))
    monkeypatch.setattr(sys, "argv", ["run_wbm_iem.py", "toy", "--model-file", str(model_file), "--protocol", str(protocol)])
    original = HighsWBM.solve

    def interrupted_solve(self):
        answer = original(self)
        if self.n_solves == 4:
            return "Time limit reached", answer[1], answer[2], answer[3]
        return answer

    monkeypatch.setattr(HighsWBM, "solve", interrupted_solve)
    runner.main()
    output = tmp_path / "toy_iem_results_v0.2.json"
    first = json.loads(output.read_text())
    assert len(first) == 1
    assert first[0]["status"] == "partial"
    assert first[0]["biomarkers"][0]["healthy"] is None
    assert "NaN" not in output.read_text()

    monkeypatch.setattr(HighsWBM, "solve", original)
    runner.main()
    second = json.loads(output.read_text())
    assert len(second) == 1
    assert second[0]["status"] == "complete"
    assert second[0]["biomarkers"][0]["correct"] is True
    runner.main()
    summary = json.loads((tmp_path / "toy_iem_summary_v0.2.json").read_text())
    assert summary["session_solves"] == 0
    assert summary["n_biomarkers_scored"] == summary["n_biomarkers_expected_in_attempted_iems"] == 1
    assert summary["total_solves_in_recorded_attempts"] == 5
