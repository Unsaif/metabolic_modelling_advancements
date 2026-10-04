"""Tests for the port of the COBRA Toolbox whole-body constraint setup (gembench/wbm_constraints.py)."""
import json
import os

import numpy as np
import pytest
import scipy.sparse as sp

from gembench import wbm_constraints as C
from gembench import wbm_iem as I

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LITERATURE_RXNS = ["Muscle_EX_glc_D(e)_[bc]", "Muscle_EX_ala_L(e)_[bc]", "RBC_EX_glc_D(e)_[bc]",
                   "Brain_EX_glc_D(e)_[csf]", "Brain_EX_o2(e)_[csf]", "Liver_EX_ala_L(e)_[bc]",
                   "Liver_EX_glc_D(e)_[bc]", "Adipocytes_EX_glyc(e)_[bc]"]


def model(spec):
    """spec: {rxn: ({met: coef}, lb, ub)} -> rxns, mets, S, lb, ub. Literature reactions are added if absent."""
    spec = dict(spec)
    for r in LITERATURE_RXNS:
        if r not in spec:
            organ = r.split("_EX_")[0]
            met = r.split("_EX_")[1].split("(e)")[0]
            comp = "[csf]" if r.endswith("[csf]") else "[bc]"
            spec[r] = ({f"{organ}_{met}[e]": -1, f"lit_{met}{comp}": 1}, -1e6, 1e6)
    rxns = list(spec)
    mets = sorted({m for st, _, _ in spec.values() for m in st})
    S = sp.lil_matrix((len(mets), len(rxns)))
    for j, r in enumerate(rxns):
        for met, v in spec[r][0].items():
            S[mets.index(met), j] = v
    lb = np.array([spec[r][1] for r in rxns], dtype=float)
    ub = np.array([spec[r][2] for r in rxns], dtype=float)
    return rxns, mets, S.tocsc(), lb, ub


def inputs(blood=(), csf=(), urine=()):
    def rows(items):
        return [{"vmh": v, "min": lo, "max": hi, "unit": "uM", "hmdb": "HMDB"} for v, lo, hi in items]
    flow = {"Muscle": 0.15, "Liver": 0.1, "Kidney": 0.2, "Brain": 0.1, "Scord": 0.01, "Adipocytes": 0.05, "Heart": None}
    return {"hmdb": {"blood": rows(blood), "csf": rows(csf), "urine": rows(urine)},
            "blood_flow": {"organs": {k: {"male": v, "female": v} for k, v in flow.items()}},
            "organ_weights": {s: {"body_weight_g": 70000.0, "organs": [{"name": "Brain", "weight_g": 1400.0, "fraction": 0.02}]}
                              for s in ("male", "female")},
            "diet_eu_average": [], "agora_essential": [], "provenance": {"toolbox_commit": "test"}}


def apply(spec, inp):
    rxns, mets, S, lb, ub = model(spec)
    params = C.default_parameters(inp, "male")
    lb2, ub2, log = C.physiological_constraints(rxns, mets, S, lb, ub, params, inp)
    return dict(zip(rxns, lb2)), dict(zip(rxns, ub2)), log


def ex(organ, met, comp="[bc]"):
    return {f"{organ}_{met}[e]": -1, f"{met}{comp}": 1}


def test_plasma_flows_use_listed_fraction_one_percent_default_and_brain_plus_spinal_cord_for_bbb():
    pf = C.plasma_flow_rates(C.default_parameters(inputs(), "male"))
    assert pf["Muscle"] == pytest.approx(0.15 * 5360 * 0.6)
    assert pf["Heart"] == pytest.approx(0.01 * 5360 * 0.6)        # empty cell -> 1 percent
    assert pf["RBC"] == pytest.approx(0.01 * 5360 * 0.6)          # absent organ -> 1 percent
    assert pf["BBB"] == pytest.approx(0.11 * 5360 * 0.6)


def test_kidney_filtration_and_reabsorption_from_gfr():
    inp = inputs(blood=[("glc_D", 3900, 6100), ("urea", 2500, 7000), ("lowx", 0, 0.0005)])
    spec = {"Kidney_EX_glc_D(e)_[bc]": (ex("Kidney", "glc_D"), -1e6, 1e6),
            "Kidney_EX_urea(e)_[bc]": (ex("Kidney", "urea"), -1e6, 1e6),
            "Kidney_EX_lowx(e)_[bc]": (ex("Kidney", "lowx"), -1e6, 1e6),
            "Kidney_EX_nodata(e)_[bc]": (ex("Kidney", "nodata"), -1e6, 1e6),
            "Kidney_EX_y(e)_[bcK]": (ex("Kidney", "y", "[bcK]"), -1e6, 1e6)}
    lb, ub, _ = apply(spec, inp)
    day = 90 * 60 * 24 / 1000          # GFR in L/day
    assert ub["Kidney_EX_glc_D(e)_[bc]"] == pytest.approx(-3.9 * day)   # forced filtration of the minimum
    assert lb["Kidney_EX_glc_D(e)_[bc]"] == pytest.approx(-6.1 * day)
    assert ub["Kidney_EX_urea(e)_[bc]"] == 0                             # listed: no forced filtration
    assert lb["Kidney_EX_urea(e)_[bc]"] == pytest.approx(-7.0 * day)
    assert ub["Kidney_EX_lowx(e)_[bc]"] == 0                             # maximum below 50 uM
    assert ub["Kidney_EX_nodata(e)_[bc]"] == 0 and lb["Kidney_EX_nodata(e)_[bc]"] == pytest.approx(-0.02 * day)
    assert lb["Kidney_EX_y(e)_[bcK]"] == -1e6                            # no '[bc]' in the name: untouched


def test_blood_uptake_limits_only_open_uptakes_and_floors_tiny_concentrations():
    inp = inputs(blood=[("glc_D", 3900, 6100), ("tiny", 0, 0.0005), ("o2", 0, 100)])
    spec = {"Muscle_EX_x(e)_[bc]": (ex("Muscle", "glc_D"), -1e6, 1e6),
            "Muscle_EX_y(e)_[bc]": (ex("Muscle", "glc_D"), 0, 1e6),
            "Muscle_EX_tiny(e)_[bc]": (ex("Muscle", "tiny"), -1e6, 1e6),
            "Muscle_EX_none(e)_[bc]": (ex("Muscle", "none"), -1e6, 1e6),
            "Muscle_EX_o2(e)_[bc]": ({"Muscle_o2[e]": -1, "RBC_o2[bc]": 1}, -1e6, 1e6),
            "Muscle_EX_h2o(e)_[bc]": (ex("Muscle", "h2o"), -1e6, 1e6)}
    lb, ub, _ = apply(spec, inp)
    pf_day = 0.15 * 5360 * 0.6 * 60 * 24 / 1000
    assert lb["Muscle_EX_x(e)_[bc]"] == pytest.approx(-6.1 * pf_day)
    assert lb["Muscle_EX_y(e)_[bc]"] == 0                       # no uptake allowed: unchanged
    assert lb["Muscle_EX_tiny(e)_[bc]"] == pytest.approx(-1e-6 * pf_day)   # floor of 1 nM
    assert lb["Muscle_EX_none(e)_[bc]"] == pytest.approx(-0.02 * pf_day)   # default 20 uM
    assert lb["Muscle_EX_o2(e)_[bc]"] == pytest.approx(-0.1 * pf_day)      # RBC_ fallback
    assert lb["Muscle_EX_h2o(e)_[bc]"] == -1e6                  # water excluded
    assert ub["Muscle_EX_x(e)_[bc]"] == 1e6                     # secretion is never limited


def test_bbb_uptake_uses_brain_and_spinal_cord_flow_and_listed_exceptions():
    inp = inputs(blood=[("glc_D", 3900, 6100)])
    spec = {"BBB_GLC_D[CSF]upt": ({"glc_D[csf]": -1, "glc_D[bc]": 1}, -1e6, 0),
            "BBB_NH4[CSF]upt": ({"nh4[csf]": -1, "nh4[bc]": 1}, -1e6, 0),
            "BBB_H2O[CSF]upt": ({"h2o[csf]": -1, "h2o[bc]": 1}, -1e6, 0)}
    lb, _, _ = apply(spec, inp)
    assert lb["BBB_GLC_D[CSF]upt"] == pytest.approx(-6.1 * 0.11 * 5360 * 0.6 * 1.44)
    assert lb["BBB_NH4[CSF]upt"] == -1e6
    assert lb["BBB_H2O[CSF]upt"] == -1e6


def test_csf_export_bounds_and_legacy_flow_variant():
    inp = inputs(csf=[("x", 10, 100), ("na1", 100000, 150000), ("lo", 1, 100)])
    spec = {"BBB_X[CSF]exp": ({"x[csf]": -1, "x[bc]": 1}, -1e6, 1e6),
            "BBB_NA1[CSF]exp": ({"na1[csf]": -1, "na1[bc]": 1}, -1e6, 1e6),
            "BBB_LO[CSF]exp": ({"lo[csf]": -1, "lo[bc]": 1}, -1e6, 1e6),
            "BBB_NONE[CSF]exp": ({"none[csf]": -1, "none[bc]": 1}, -1e6, 1e6)}
    lb, ub, _ = apply(spec, inp)
    day = 0.52 * 60 * 24 / 1000
    assert lb["BBB_X[CSF]exp"] == pytest.approx(0.010 * day) and ub["BBB_X[CSF]exp"] == pytest.approx(0.1 * day)
    assert lb["BBB_NA1[CSF]exp"] == 0 and ub["BBB_NA1[CSF]exp"] == pytest.approx(150 * day)
    assert lb["BBB_LO[CSF]exp"] == 0                         # minimum below 5 uM
    assert lb["BBB_NONE[CSF]exp"] == 0 and ub["BBB_NONE[CSF]exp"] == pytest.approx(0.02 * day)
    rxns, mets, S, lb0, ub0 = model(spec)
    params = C.default_parameters(inp, "male")
    params.csf_export_ub_flow_rate = params.csf_flow_rate
    _, ub_legacy, _ = C.physiological_constraints(rxns, mets, S, lb0, ub0, params, inp)
    assert dict(zip(rxns, ub_legacy))["BBB_X[CSF]exp"] == pytest.approx(0.1 * 0.35 * 1.44)


def test_urine_excretion_only_where_secretion_is_open():
    inp = inputs(urine=[("x", 10, 100), ("k", 30000, 90000), ("low", 1, 100)])
    spec = {"EX_x[u]": ({"x[u]": -1}, 0, 1e6), "EX_closed[u]": ({"x[u]": -1}, 0, 0),
            "EX_k[u]": ({"k[u]": -1}, 0, 1e6), "EX_low[u]": ({"low[u]": -1}, 0, 1e6),
            "EX_none[u]": ({"none[u]": -1}, 0, 1e6)}
    lb, ub, _ = apply(spec, inp)
    cr_min, cr_max = 0.5 * 10 / 113.1179, 1.2 * 10 / 113.1179
    assert lb["EX_x[u]"] == pytest.approx(0.010 * cr_min * 2) and ub["EX_x[u]"] == pytest.approx(0.1 * cr_max * 2)
    assert ub["EX_closed[u]"] == 0 and lb["EX_closed[u]"] == 0
    assert lb["EX_k[u]"] == 0 and ub["EX_k[u]"] == pytest.approx(90 * cr_max * 2)   # listed: no forced excretion
    assert lb["EX_low[u]"] == 0
    assert lb["EX_none[u]"] == 0 and ub["EX_none[u]"] == pytest.approx(0.02 * cr_max * 2)


def test_literature_constraints_respect_existing_ranges():
    inp = inputs(blood=[("glc_D", 3900, 6100), ("ala_L", 100, 500)])
    spec = {"Muscle_EX_glc_D(e)_[bc]": (ex("Muscle", "glc_D"), -1e6, 1e6),
            "Liver_EX_glc_D(e)_[bc]": (ex("Liver", "glc_D"), -1e6, 1e6),
            "Liver_EX_ala_L(e)_[bc]": (ex("Liver", "ala_L"), -1e6, 1e6),
            "Heart_EX_co2(e)_[bc]": (ex("Heart", "co2"), -1e6, 1e6),
            "Liver_EX_co2(e)_[bc]": (ex("Liver", "co2"), -1e6, 1e6),
            "Heart_DM_atp_c_": ({"Heart_atp[c]": -1}, 0, 1e6),
            "EX_o2[a]": ({"o2[a]": -1}, -1e6, 1e6)}
    lb, ub, log = apply(spec, inp)
    assert ub["Muscle_EX_glc_D(e)_[bc]"] == -10                     # muscle must take up glucose
    glc = 130 * 60 * 24 * 70 / 65 / 1000 * 1000 / 180.16          # hepatic glucose output, mmol/day
    assert (lb["Liver_EX_glc_D(e)_[bc]"], ub["Liver_EX_glc_D(e)_[bc]"]) == (pytest.approx(0.8 * glc), pytest.approx(1.2 * glc))
    assert glc * 0.8 == pytest.approx(895.204, abs=1e-3)           # value stored in Harvey 1.03d
    # hepatic alanine uptake is already tighter than the literature range from the blood step: unchanged
    assert lb["Liver_EX_ala_L(e)_[bc]"] == pytest.approx(-0.5 * 0.1 * 5360 * 0.6 * 1.44)
    assert lb["Heart_EX_co2(e)_[bc]"] == 0 and lb["Liver_EX_co2(e)_[bc]"] == -10000
    assert lb["Heart_DM_atp_c_"] == 6000
    assert (lb["EX_o2[a]"], ub["EX_o2[a]"]) == (-25000, -15000)
    assert "Reaction EX_co2[a] not in model" in log.warnings
    brain_o2 = 156 * 60 * 24 * 1400 / 100 / 1000
    assert lb["Brain_EX_o2(e)_[csf]"] == pytest.approx(-1.2 * brain_o2)
    assert brain_o2 * 1.2 == pytest.approx(3773.952)


def test_missing_literature_reaction_is_an_error_as_in_matlab():
    rxns, mets, S, lb, ub = model({})
    keep = [i for i, r in enumerate(rxns) if r != "Liver_EX_glc_D(e)_[bc]"]
    inp = inputs()
    with pytest.raises(KeyError):
        C.physiological_constraints([rxns[i] for i in keep], mets, S[:, keep], lb[keep], ub[keep],
                                    C.default_parameters(inp, "male"), inp)


def test_diet_constraints_follow_set_diet_constraints():
    rxns = ["Diet_EX_glc_D[d]", "Diet_EX_fol[d]", "Diet_EX_thm[d]", "Diet_EX_na1[d]", "Diet_EX_so4[d]",
            "Diet_EX_ac[d]", "Diet_EX_asn_L[d]", "Diet_EX_chol[d]", "Diet_EX_other[d]", "EX_glc_D[u]"]
    lb = np.full(len(rxns), -5.0); ub = np.full(len(rxns), 5.0)
    diet = [["Diet_EX_glc_D[d]", "100"], ["Diet_EX_fol[d]", "0.0001"], ["Diet_EX_thm[d]", "20"],
            ["Diet_EX_na1[d]", "170"], ["Diet_EX_so4[d]", "3"], ["Diet_EX_absent[d]", "1"]]
    agora = ["Diet_EX_ac[d]", "Diet_EX_glc_D[d]"]
    lb2, ub2, log = C.diet_constraints(rxns, lb, ub, diet, agora)
    got = {r: (l, u) for r, l, u in zip(rxns, lb2, ub2)}
    assert got["Diet_EX_glc_D[d]"] == (pytest.approx(-120), pytest.approx(-80))
    assert got["Diet_EX_fol[d]"] == (pytest.approx(-10), pytest.approx(-0.00008))      # micronutrient <= 10
    assert got["Diet_EX_thm[d]"] == (pytest.approx(-2400), pytest.approx(-16))         # micronutrient > 10
    assert got["Diet_EX_na1[d]"] == (pytest.approx(-20400), pytest.approx(-136))       # ion
    assert got["Diet_EX_so4[d]"] == (pytest.approx(-1000), pytest.approx(-2.4))
    assert got["Diet_EX_ac[d]"] == (-0.1, 0)                                           # AGORA essential, not in diet
    assert got["Diet_EX_asn_L[d]"] == (-50, 0)
    assert got["Diet_EX_chol[d]"] == (-41.251, 0)
    assert got["Diet_EX_other[d]"] == (0, 0)
    assert got["EX_glc_D[u]"] == (-5, 5)
    assert any("Diet_EX_absent[d]" in w for w in log.warnings)


def test_bile_duct_global_constraint_uses_the_toolbox_list():
    assert len(I.RUNIEM_BILE_DUCT_UB100) == len(set(I.RUNIEM_BILE_DUCT_UB100)) == 28
    src = os.path.join(ROOT, "external", "cobratoolbox", "src", "analysis", "wholeBody", "PSCMToolbox", "runIEM_HH.m")
    if os.path.exists(src):
        import re
        listed = re.findall(r"'([^']+)'", re.search(r"Rnew\s*=\s*\{(.*?)\};", open(src).read(), flags=re.S).group(1))
        assert listed == I.RUNIEM_BILE_DUCT_UB100


def test_extracted_inputs_match_the_reference_man():
    path = C.DEFAULT_INPUTS
    data = json.load(open(path))
    assert data["provenance"]["toolbox_commit"].startswith("67c790d")
    assert (len(data["hmdb"]["blood"]), len(data["hmdb"]["csf"]), len(data["hmdb"]["urine"])) == (475, 232, 253)
    params = C.default_parameters(data, "male")
    assert params.body_weight_kg == 70 and params.organ_weights_g["Brain"] == 1400
    assert params.blood_flow_fraction["Kidney"] == pytest.approx(0.201724138)
    assert C.legacy_gfr(params) == pytest.approx(129.749, abs=1e-3)
    assert len(data["diet_eu_average"]) == 95 and data["diet_eu_average"][0] == ["Diet_EX_10fthf[d]", "0.000125752"]
    assert len(data["agora_essential"]) == 137 and all(a.startswith("Diet_EX_") and a.endswith("[d]") for a in data["agora_essential"])


def test_matlab_cell_parser_drops_comments():
    import importlib.util
    spec = importlib.util.spec_from_file_location("x", os.path.join(ROOT, "scripts", "extract_wbm_constraint_inputs.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    text = "Diet = {%'EX_no(e)' '1'\n'EX_a(e)'\t'0.5' % trailing 'EX_b(e)'\n% 'EX_c(e)' '2'\n'EX_d(e)' '3.6E-06'\n};"
    assert mod.matlab_cell_strings(text, "Diet") == ["EX_a(e)", "0.5", "EX_d(e)", "3.6E-06"]
    assert mod.matlab_quoted_number("'0.05'") == 0.05 and mod.matlab_quoted_number("''") is None


def test_duplicate_concentration_rows_use_the_first_as_matlab_does():
    inp = inputs(blood=[("glc_D", 100, 200), ("glc_D", 300, 400)])
    lb, _, _ = apply({"Muscle_EX_x(e)_[bc]": (ex("Muscle", "glc_D"), -1e6, 1e6)}, inp)
    assert lb["Muscle_EX_x(e)_[bc]"] == pytest.approx(-0.2 * 0.15 * 5360 * 0.6 * 1.44)
