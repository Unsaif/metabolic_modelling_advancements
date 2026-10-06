"""Check 3a: why iAH991 cannot grow under the protocol, and what the exploratory supplement does.

Claims tested (Paper 1, Results): the biomass needs adenosylcobalamin; the model cannot make B12; B12 is absent
from the screening medium; there is no route from sulfate to cysteine; cysteine (and methionine) are capped at
0.001. Also: whether the two open sinks (sink_chols, sink_s) or the three demands matter, whether the exploratory
medium supports growth without a carbon source, and what limits growth in the exploratory run.

All growth values use this check's own protocol implementation (protocol.py) on the BiGG view after medium
completion, D-glucose at -10 unless stated.

Usage: python -I check_iah991_growth.py   (writes check_iah991_growth.json)
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_model  # noqa: E402
from protocol import apply_condition, complete_medium, growth, mapped_conditions  # noqa: E402

import cobra  # noqa: E402

OUT = os.path.dirname(os.path.abspath(__file__))
MED = "Varel_Bryant_medium"
SUP = {"cbl1": -0.001, "h2s": -1000.0}


def grow_with(m, carbon=("glc__D",), supplement=None, extra=None, close=()):
    with m:
        apply_condition(m, MED, list(carbon), supplement)
        for rid, lb in (extra or {}).items():
            m.reactions.get_by_id(rid).lower_bound = lb
        for rid in close:
            m.reactions.get_by_id(rid).lower_bound = 0.0
        return growth(m)


def max_production(m, met, open_all=True, close=()):
    with m:
        if open_all:
            for r in m.reactions:
                if r.boundary:
                    r.lower_bound = -1000.0
        for rid in close:
            m.reactions.get_by_id(rid).bounds = (0.0, 1000.0) if rid.startswith(("EX_", "sink_")) else (0.0, 0.0)
        dm = cobra.Reaction("CHECK_DM", lower_bound=0.0, upper_bound=1000.0)
        dm.add_metabolites({m.metabolites.get_by_id(met): -1.0})
        m.add_reactions([dm])
        m.objective = dm
        return m.slim_optimize(error_value=float("nan"))


def main():
    m = read_model("models/curated/iAH991/iAH991_bigg_view.xml.gz")
    added = complete_medium(m, [MED], SUP)
    rep = {"medium_completion_added_with_supplement": added}
    bm = m.reactions.get_by_id("Biomass_BT_v2")
    rep["biomass_adocbl_coefficient"] = bm.metabolites.get(m.metabolites.get_by_id("adocbl_c"))
    rep["growth_cap_from_B12_at_0.001"] = 0.001 / abs(rep["biomass_adocbl_coefficient"])
    # B12 de novo: every boundary open except cobalamin uptakes
    b12_ex = [r.id for r in m.exchanges if any(x.id in ("cbl1_e", "adocbl_e", "cbl2_e", "cbi_e") for x in r.metabolites)]
    rep["cobalamin_exchanges"] = b12_ex
    rep["max_adocbl_without_cobalamin_uptake"] = max_production(m, "adocbl_c", close=b12_ex)
    rep["max_adocbl_with_cbl1"] = max_production(m, "adocbl_c")
    # sulfur: cysteine from sulfate only (every boundary open except other sulfur sources)
    s_sources = [r.id for r in m.reactions if r.boundary and any("S" in (x.formula or "") for x in r.metabolites)
                 and not any(x.id.startswith("so4") for x in r.metabolites)]
    rep["sulfur_boundaries_closed_for_sulfate_test"] = s_sources
    rep["max_cysteine_from_sulfate_only"] = max_production(m, "cys__L_c", close=s_sources)
    rep["max_cysteine_with_h2s"] = max_production(m, "cys__L_c", close=[x for x in s_sources if x != "EX_h2s_e"])
    # transsulfuration (methionine -> cysteine)?
    rep["max_cysteine_from_methionine_only"] = max_production(m, "cys__L_c", close=[x for x in s_sources if x != "EX_met__L_e"])
    # growth dissection on D-glucose
    g = {}
    g["strict protocol"] = grow_with(m)
    g["+ B12 only (cbl1 -0.001)"] = grow_with(m, supplement={"cbl1": -0.001})
    g["+ sulfide only (h2s -1000)"] = grow_with(m, supplement={"h2s": -1000.0})
    g["+ B12 and sulfide (exploratory)"] = grow_with(m, supplement=SUP)
    g["+ B12 (-0.001), cysteine uncapped (-10)"] = grow_with(m, supplement={"cbl1": -0.001}, extra={"EX_cys__L_e": -10.0})
    g["+ B12 (-0.001), methionine uncapped (-10)"] = grow_with(m, supplement={"cbl1": -0.001}, extra={"EX_met__L_e": -10.0})
    g["+ B12 uncapped (-10), sulfide"] = grow_with(m, supplement={"cbl1": -10.0, "h2s": -1000.0})
    g["exploratory, no carbon source"] = grow_with(m, carbon=(), supplement=SUP)
    g["exploratory, no carbon source, sinks closed"] = grow_with(m, carbon=(), supplement=SUP, close=("sink_chols", "sink_s"))
    g["exploratory + glucose, sinks closed"] = grow_with(m, supplement=SUP, close=("sink_chols", "sink_s"))
    g["exploratory + glucose, sink_s closed"] = grow_with(m, supplement=SUP, close=("sink_s",))
    g["exploratory + glucose, sink_chols closed"] = grow_with(m, supplement=SUP, close=("sink_chols",))
    g["exploratory, no carbon, B12 and sulfide uncapped"] = grow_with(m, carbon=(), supplement={"cbl1": -10.0, "h2s": -1000.0})
    rep["growth_on_glucose"] = g
    # sinks: what flows through them at the exploratory glucose optimum
    with m:
        apply_condition(m, MED, ["glc__D"], SUP)
        sol = m.optimize()
        rep["exploratory_glucose_fluxes"] = {r: float(sol.fluxes[r]) for r in ("sink_chols", "sink_s", "DM_4HBA", "DM_5DRIB", "DM_AMOB",
                                                                               "EX_cbl1_e", "EX_h2s_e", "EX_cys__L_e", "EX_met__L_e",
                                                                               "EX_glc__D_e", "Biomass_BT_v2")}
        rep["chol_c_reactions"] = [f"{r.id}: {r.reaction} {r.bounds}" for r in m.metabolites.get_by_id("chol_c").reactions]
    # every mapped condition: exploratory growth, and growth with B12 uncapped (is growth B12-limited?)
    rows = []
    for c in mapped_conditions("Btheta"):
        a = grow_with(m, carbon=c["bigg_ids"], supplement=SUP)
        b = grow_with(m, carbon=c["bigg_ids"], supplement={"cbl1": -10.0, "h2s": -1000.0})
        rows.append({"condition": c["name"], "bigg_ids": c["bigg_ids"], "exploratory": a, "B12_uncapped": b,
                     "B12_limited": a > 1e-3 and b > a * 1.001})
    rep["conditions"] = rows
    rep["n_conditions_B12_limited"] = sum(r["B12_limited"] for r in rows)
    with open(os.path.join(OUT, "check_iah991_growth.json"), "w") as fh:
        json.dump(rep, fh, indent=1, default=float)
    for k, v in rep.items():
        if k != "conditions":
            print(k, json.dumps(v, default=float)[:700])
    for r in rows:
        print(f"  {r['condition'][:34]:34s} exploratory={r['exploratory']:.4f}  B12 uncapped={r['B12_uncapped']:.4f}  B12-limited={r['B12_limited']}")


if __name__ == "__main__":
    main()
