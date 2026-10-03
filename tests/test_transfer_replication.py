"""Tests for the panel B replication wrapper (scripts/transfer_replication.py)."""
import importlib.util
import json
import os

import cobra
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load():
    spec = importlib.util.spec_from_file_location("transfer_replication", os.path.join(ROOT, "scripts", "transfer_replication.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


TR = _load()


def _rxn(stoich):
    r = cobra.Reaction("R")
    r.add_metabolites({cobra.Metabolite(m, compartment=m.rsplit("_", 1)[1]): v for m, v in stoich.items()})
    return r


@pytest.mark.parametrize("stoich, expected", [
    ({"zn2_p": -1, "atp_c": -1, "h2o_c": -1, "zn2_c": 1, "adp_c": 1, "pi_c": 1, "h_c": 1}, True),    # zinc ABC transporter
    ({"cu2_e": -1, "cu2_c": 1}, True),                                                                # copper uptake
    ({"pi_p": -1, "atp_c": -1, "h2o_c": -1, "pi_c": 2, "adp_c": 1, "h_c": 1}, True),                  # phosphate ABC
    ({"glc__D_p": -1, "atp_c": -1, "h2o_c": -1, "glc__D_c": 1, "adp_c": 1, "pi_c": 1, "h_c": 1}, False),  # sugar ABC
    ({"glc__D_e": -1, "na1_e": -1, "glc__D_c": 1, "na1_c": 1}, False),                                # Na+/glucose symport
    ({"h_p": -1, "h_c": 1}, False),                                                                   # protons only
    ({"pyr_c": -1, "acald_c": 1, "co2_c": 1}, False),                                                 # not a transport
])
def test_inorganic_ion_transport_classification(stoich, expected):
    assert TR.is_inorganic_ion_transport(_rxn(stoich)) is expected


def test_replication_arms_are_the_frozen_v1_arms_plus_one():
    v1 = json.load(open(os.path.join(ROOT, "data", "studies", "transfer_v1", "arms.json")))["arms"]
    rep = json.load(open(os.path.join(ROOT, "data", "studies", "transfer_v1_replication", "arms.json")))["arms"]
    assert set(rep) - set(v1) == {"M_noIonR6"}
    assert all(rep[k] == v1[k] for k in v1)
    assert rep["M_noIonR6"]["transforms"][:-1] == v1["UNQ"]["transforms"]


def test_pool_rejects_duplicate_organisms(tmp_path):
    a = [{"org": "x", "A": "B0", "B": "UNQ", "union_delta": 0.1}]
    (tmp_path / "a.json").write_text(json.dumps(a))
    (tmp_path / "b.json").write_text(json.dumps(a))
    with pytest.raises(SystemExit):
        TR.cmd_pool(f"{tmp_path / 'a.json'},{tmp_path / 'b.json'}", str(tmp_path / "out.json"))
    b = [{"org": "y", "A": "B0", "B": "UNQ", "union_delta": 0.2}]
    (tmp_path / "b.json").write_text(json.dumps(b))
    TR.cmd_pool(f"{tmp_path / 'a.json'},{tmp_path / 'b.json'}", str(tmp_path / "out.json"))
    assert [r["org"] for r in json.load(open(tmp_path / "out.json"))] == ["x", "y"]
