"""Tests for transfer study v1 components (gembench.transfer, gembench.feba_media, scripts/transfer_*.py)."""
import importlib.util
import os
import sys

import cobra
import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from gembench import feba_media as FM  # noqa: E402
from gembench import transfer as T  # noqa: E402
from gembench.gene_mapping import GeneMap  # noqa: E402


def _load_script(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, "scripts", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --- media -----------------------------------------------------------------------------------------------------------

MEDIA = "\n".join([
    "Media\tTestMM_noC\t", "Description\ttest medium\t", "Minimal\tTRUE\t", "Controlled vocabulary\tConcentration\tUnits",
    "Ammonium chloride\t1\tg/L", "PIPES sesquisodium salt\t30\tmM", "Thiamine HCl\t0.01\tmM", "Test mix\t1\tX", "\t\t",
    "Media\tTestRich\t", "Description\trich\t", "Minimal\tFALSE\t", "Controlled vocabulary\tConcentration\tUnits",
    "Yeast Extract\t5\tg/L", "Sodium Chloride\t5\tg/L", "\t\t"]) + "\n"
MIXES = "\n".join([
    "Media\tTest mix\t", "Description\t100X\t", "X\t100\t", "Controlled vocabulary\tConcentration\tUnits",
    "Calcium Chloride Dihydrate\t0.1\tg/L", "calcium pantothenate\t5\tmg/L", "\t\t"]) + "\n"
COMPONENTS = "\n".join([
    "# test dictionary", "component\tclass\tbigg_ids\tcounter_ions\tnote",
    "Ammonium chloride\tinorganic\tnh4;cl\t\t", "PIPES sesquisodium salt\tignore\t\tna1\tbuffer",
    "Thiamine HCl\ttrace_organic\tthm\tcl\tvitamin", "Calcium Chloride Dihydrate\tinorganic\tca2;cl\t\t",
    "calcium pantothenate\ttrace_organic\tpnto__R\tca2\tvitamin", "Yeast Extract\tundefined\t\t\t",
    "Sodium Chloride\tinorganic\tna1;cl\t\t"]) + "\n"


@pytest.fixture
def feba(tmp_path):
    (tmp_path / "media").write_text(MEDIA)
    (tmp_path / "mixes").write_text(MIXES)
    (tmp_path / "comp.tsv").write_text(COMPONENTS)
    return FM.parse_blocks(str(tmp_path / "media")), FM.parse_blocks(str(tmp_path / "mixes")), \
        FM.load_component_table(str(tmp_path / "comp.tsv"))


def test_medium_row_expands_mixes_drops_buffers_and_limits_only_trace_moieties(feba):
    media, mixes, comps = feba
    assert media["TestMM_noC"]["minimal"] is True and mixes["Test mix"]["x"] == "100"
    row = FM.medium_row("TestMM_noC", media, mixes, comps, aerobic=True)
    ids = row["bigg_ids"].split(";")
    assert row["status"] == "mapped" and row["aerobic"] == "yes"
    assert {"h2o", "h", "co2", "nh4", "cl", "na1", "thm", "ca2", "pnto__R"} <= set(ids)
    assert set(FM.TRACE_METALS) <= set(ids)
    # trace organics are limited; their counter-ions (cl, ca2) are not
    assert set(row["trace_components"].split(";")) == {"thm", "pnto__R"}
    assert row["ignored"] == ["PIPES sesquisodium salt"]


def test_medium_with_complex_ingredient_is_undefined_and_unknown_components_are_reported(feba):
    media, mixes, comps = feba
    assert FM.medium_row("TestRich", media, mixes, comps, aerobic=True)["status"] == "undefined"
    assert FM.medium_row("absent", media, mixes, comps, aerobic=True)["status"] == "not_found"
    media["TestMM_noC"]["components"].append(("Mystery salt", "1", "g/L"))
    row = FM.medium_row("TestMM_noC", media, mixes, comps, aerobic=False)
    assert row["status"] == "unmapped_components" and row["unmapped"] == ["Mystery salt"] and row["aerobic"] == "no"


def test_component_dictionaries_cannot_redefine_a_component(tmp_path):
    (tmp_path / "a.tsv").write_text(COMPONENTS)
    (tmp_path / "b.tsv").write_text("component\tclass\tbigg_ids\tcounter_ions\tnote\nammonium chloride\tinorganic\tnh4\t\t\n")
    with pytest.raises(ValueError, match="more than once"):
        FM.load_component_table([str(tmp_path / "a.tsv"), str(tmp_path / "b.tsv")])


def test_aerobic_flag_is_the_majority_over_carbon_source_experiments():
    e = pd.DataFrame({"media": ["M", "M", "M", "N"], "expGroup": ["carbon source"] * 3 + ["stress"],
                      "aerobic": ["Aerobic", "Aerobic", "Anaerobic", "Anaerobic"]})
    assert FM.aerobic_flag(e, "M") is True
    assert FM.aerobic_flag(e, "N") is None


def test_rule_reproduces_development_media_components():
    media = FM.parse_blocks(os.path.join(ROOT, "data/studies/transfer_v1/feba/media"))
    mixes = FM.parse_blocks(os.path.join(ROOT, "data/studies/transfer_v1/feba/mixes"))
    comps = FM.load_component_table(os.path.join(ROOT, "data/studies/transfer_v1/feba_component_bigg.tsv"))
    ref = pd.read_table(os.path.join(ROOT, "data/reference/fitness_browser_media_bigg.tsv"), dtype=str,
                        keep_default_na=False).set_index("media")
    for name in ref.index:
        aerobic = ref.at[name, "aerobic"] == "yes"
        row = FM.medium_row(name, media, mixes, comps, aerobic)
        assert row["status"] == "mapped", name
        assert set(row["bigg_ids"].split(";")) == set(ref.at[name, "bigg_ids"].split(";")), name


# --- reference condition -----------------------------------------------------------------------------------------------

def _c(name, media, ids, n):
    return {"key": f"{name} | {media}", "name": name, "media": media, "bigg_ids": ids, "n_experiments": n}


def test_reference_condition_uses_main_medium_then_central_substrates_then_most_experiments():
    ux = {"EX_glc__D_e", "EX_lac__L_e", "EX_pyr_e", "EX_adn_e", "EX_cytd_e"}
    conds = [_c("D-Glucose", "Minor", ["glc__D"], 5), _c("Pyruvate", "Main", ["pyr"], 2), _c("L-Lactate", "Main", ["lac__L"], 1),
             _c("Adenosine", "Main", ["adn"], 3), _c("Cytidine", "Main", ["cytd"], 3), _c("Mix", "Main", ["glc__D", "pyr"], 1)]
    ref = T.choose_reference_condition(conds, ux)
    assert ref["name"] == "L-Lactate"          # glucose only on the minor medium; L-lactate precedes pyruvate
    no_central = [c for c in conds if c["name"] not in ("Pyruvate", "L-Lactate")]
    assert T.choose_reference_condition(no_central, ux)["name"] == "Adenosine"   # most experiments, then name
    assert T.choose_reference_condition([], ux) is None


# --- gene-rule transforms --------------------------------------------------------------------------------------------

def _model(rules):
    m = cobra.Model("t")
    a = cobra.Metabolite("a_c", compartment="c")
    b = cobra.Metabolite("b_c", compartment="c")
    for i, rule in enumerate(rules):
        r = cobra.Reaction(f"R{i}")
        r.add_metabolites({a: -1, b: 1})
        m.add_reactions([r])
        r.gene_reaction_rule = rule
    return m


def test_conjunction_normalisation_drops_standalone_members_of_a_complex():
    m = _model(["A or B or (A and B)", "A or C or (A and B)", "A or C", "(A and B) or (C and D)"])
    T.normalize_conjunction_rules(m)
    rules = [m.reactions.get_by_id(f"R{i}").gene_reaction_rule for i in range(4)]
    assert rules[0] == "A and B"
    assert set(rules[1].replace("(", "").replace(")", "").split(" or ")) == {"C", "A and B"}
    assert rules[2] == "A or C" and rules[3] == "(A and B) or (C and D)"


def test_apply_decisions_adds_genome_genes_and_rejects_unknown_genes():
    m = _model(["", "g1 or g2"])
    gm = GeneMap(org_id="t", model_to_browser={"g1": "L1", "g2": "L2"}, unmapped_model_genes=[])
    dec = {"decisions": [{"reaction": "R0", "decision": "R6", "new_rule": "L7 and L8"},
                         {"reaction": "R1", "decision": "R1", "new_rule": "L1"},
                         {"reaction": "R1", "decision": "R1", "new_rule": "L99"},
                         {"reaction": "RX", "decision": "R6", "new_rule": "L7"},
                         {"reaction": "R0", "decision": "abstain", "new_rule": None}]}
    rec = T.apply_decisions(m, dec, gm, {"L1", "L2", "L7", "L8"})
    assert m.reactions.R0.gene_reaction_rule == "L7 and L8"
    assert m.reactions.R1.gene_reaction_rule == "g1"
    assert [r.get("not_applied", "") for r in rec][2:] == ["unknown genes ['L99']", "reaction absent from model"]


# --- analysis --------------------------------------------------------------------------------------------------------

def test_sign_test_and_organism_aggregate():
    agg = _load_script("transfer_aggregate")
    assert agg.sign_test_p(6, 0) == pytest.approx(2 / 64)
    assert agg.sign_test_p(3, 3) == 1.0
    rows = [{"org": o, "A": "B0", "B": "U", "union_delta": d, "union_delta_ci95": [d - 0.01, d + 0.01]}
            for o, d in [("a", 0.02), ("b", 0.01), ("c", -0.0005)]]
    rows.append({"org": "d", "A": "B0", "B": "U", "missing": True})
    res = agg.aggregate(rows, n_boot=2000)[0]
    assert res["n_evaluable"] == 3 and res["n_organisms"] == 4
    assert res["n_improved"] == 2 and res["n_unchanged"] == 1 and res["n_worsened"] == 0
    assert res["mean_delta"] == pytest.approx((0.02 + 0.01 - 0.0005) / 3)
    assert [o["status"] for o in res["per_organism"]][-1] == "missing run"


def test_union_pairs_predict_no_effect_for_genes_absent_from_a_model():
    cmp_ = _load_script("transfer_compare")
    params = {"growth_threshold": 1e-3, "fitness_threshold": -2.0}
    a = {"browser_genes": np.array(["g1"]), "conditions": np.array(["c1"]), "wt_growth": np.array([0.5]),
         "sim_growth": np.array([[0.5]]), "fitness": np.array([[0.1]]), "params": params, "benchmark": "x"}
    b = {"browser_genes": np.array(["g1", "g2"]), "conditions": np.array(["c1"]), "wt_growth": np.array([0.4]),
         "sim_growth": np.array([[0.4], [0.0]]), "fitness": np.array([[0.1], [-3.0]]), "params": params, "benchmark": "x"}
    sa, sb, fit, genes, grows = cmp_.union_pairs(a, b)
    assert genes == ["g1", "g2"] and grows == ["c1"]
    assert sa[1, 0] == 0.5 and sb[1, 0] == 0.0 and fit[1, 0] == -3.0
