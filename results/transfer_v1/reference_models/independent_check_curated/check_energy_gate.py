"""Check 4: the energy-generating-cycle gate.

(a) The published iGD1575 SBML as read by cobrapy, before any translation (identifiers cpd..._c0): close every
    boundary reaction, relax every other reaction to include zero flux (keeping its direction), add a dissipation
    reaction for each currency and maximise it. Currencies: ATP, NADH, NADPH, ubiquinol-8 (the gate's four that a
    two-compartment model has) and, as the analogue of the gate's h_p check for models without a periplasm,
    extracellular-to-cytosolic proton flow.
(b) The cycles: minimal-total-flux solutions at ATP-hydrolysis flux 1, enumerated by blocking each reaction of a
    found cycle in turn (depth 2), and the reactions whose removal alone ends ATP production from nothing.
(c) The same gate on the BiGG view after the protocol's medium completion (re-implemented here), on iSO783 as
    scored, and on the Smeli and MR1 draft arms (B0, UNQ, M; their saved models plus medium completion), to see
    whether the gate was applied the same way.

Usage: python -I check_energy_gate.py   (writes check_energy_gate.json)
"""
from __future__ import annotations

import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import MODELS, dump, medium_components, p, read_model  # noqa: E402

import cobra  # noqa: E402
from cobra.flux_analysis import pfba  # noqa: E402

CURRENCIES_BIGG = {
    "atp": {"atp_c": -1, "h2o_c": -1, "adp_c": 1, "pi_c": 1, "h_c": 1},
    "nadh": {"nadh_c": -1, "nad_c": 1, "h_c": 1},
    "nadph": {"nadph_c": -1, "nadp_c": 1, "h_c": 1},
    "q8h2": {"q8h2_c": -1, "q8_c": 1, "h_c": 2},
    "h_p": {"h_p": -1, "h_c": 1},
    "h_e (no-periplasm analogue)": {"h_e": -1, "h_c": 1},
}
# ModelSEED: cpd00002 ATP, cpd00001 H2O, cpd00008 ADP, cpd00009 Pi, cpd00067 H+, cpd00004 NADH, cpd00003 NAD,
# cpd00005 NADPH, cpd00006 NADP, cpd15561 ubiquinol-8, cpd15560 ubiquinone-8
CURRENCIES_MS = {
    "atp": {"cpd00002_c0": -1, "cpd00001_c0": -1, "cpd00008_c0": 1, "cpd00009_c0": 1, "cpd00067_c0": 1},
    "nadh": {"cpd00004_c0": -1, "cpd00003_c0": 1, "cpd00067_c0": 1},
    "nadph": {"cpd00005_c0": -1, "cpd00006_c0": 1, "cpd00067_c0": 1},
    "q8h2": {"cpd15561_c0": -1, "cpd15560_c0": 1, "cpd00067_c0": 2},
    "h_e (no-periplasm analogue)": {"cpd00067_e0": -1, "cpd00067_c0": 1},
}
GATE_KEYS = ("atp", "nadh", "nadph", "q8h2", "h_p")     # what gembench.checks.energy_from_nothing tests


def close_and_relax(model):
    for r in model.reactions:
        if r.boundary:
            r.bounds = (0.0, 0.0)
        else:
            r.bounds = (min(0.0, r.lower_bound), max(0.0, r.upper_bound))


def add_dissipation(model, stoich, rid="EGC_test"):
    dm = cobra.Reaction(rid, lower_bound=0.0, upper_bound=1000.0)
    dm.add_metabolites({model.metabolites.get_by_id(m): c for m, c in stoich.items()})
    model.add_reactions([dm])
    model.objective = dm
    model.objective_direction = "max"
    return dm


def egc_values(model, currencies):
    out = {}
    with model:
        close_and_relax(model)
        for key, st in currencies.items():
            if any(m not in model.metabolites for m in st):
                out[key] = None          # not testable in this model
                continue
            with model:
                add_dissipation(model, st)
                v = model.slim_optimize(error_value=float("nan"))
                out[key] = float(v)
    return out


def describe(model, rid, flux):
    r = model.reactions.get_by_id(rid)
    eq = r.build_reaction_string(use_metabolite_names=True)
    return {"id": rid, "name": r.name, "flux": round(float(flux), 6), "bounds": list(r.bounds),
            "equation_ids": r.reaction, "equation_names": eq, "gene_rule": r.gene_reaction_rule}


def enumerate_cycles(model, stoich, max_lps=400, depth=2):
    cycles, seen, lps = [], set(), [0]
    with model:
        close_and_relax(model)
        dm = add_dissipation(model, stoich)

        def solve(blocked):
            lps[0] += 1
            with model:
                for b in blocked:
                    model.reactions.get_by_id(b).bounds = (0.0, 0.0)
                v = model.slim_optimize(error_value=0.0)
                if v is None or v <= 1e-6:
                    return None
                dm.bounds = (1.0, 1.0)
                sol = pfba(model, fraction_of_optimum=1.0)
                return {k: x for k, x in sol.fluxes.items() if abs(x) > 1e-7 and k != dm.id}

        def search(blocked, d):
            if lps[0] >= max_lps:
                return
            sup = solve(blocked)
            if sup is None:
                return
            key = frozenset(sup)
            if key not in seen:
                seen.add(key)
                cycles.append({"blocked_to_find": sorted(blocked), "reactions": [describe(model, k, x) for k, x in sorted(sup.items())]})
            if d < depth:
                for k in sorted(sup):
                    search(blocked | {k}, d + 1)

        search(frozenset(), 0)
        union = sorted(set().union(*[{x["id"] for x in c["reactions"]} for c in cycles])) if cycles else []
        cuts = []
        for k in union:
            with model:
                model.reactions.get_by_id(k).bounds = (0.0, 0.0)
                v = model.slim_optimize(error_value=0.0)
            if v <= 1e-6:
                cuts.append(k)
    return {"n_lps": lps[0], "n_distinct_minimal_flux_supports": len(cycles), "cycles": cycles,
            "reactions_in_any_found_cycle": len(union),
            "single_reaction_cuts": [describe(model, k, 0.0) for k in cuts]}


def complete_medium(model, media, exclude=("pnto__R", "fol", "hco3")):
    """Re-implementation of the protocol's medium completion: an exchange and a gene-less uptake for each medium
    component the model has in the cytosol but cannot exchange."""
    added = []
    for medium in media:
        comps, _ = medium_components(medium)
        for c in comps:
            ex_id = f"EX_{c}_e"
            if ex_id in model.reactions or f"{c}_c" not in model.metabolites or c in exclude:
                continue
            cyt = model.metabolites.get_by_id(f"{c}_c")
            if f"{c}_e" not in model.metabolites:
                model.add_metabolites([cobra.Metabolite(f"{c}_e", name=cyt.name, compartment="e")])
            ext = model.metabolites.get_by_id(f"{c}_e")
            ex = cobra.Reaction(ex_id, lower_bound=0.0, upper_bound=1000.0)
            ex.add_metabolites({ext: -1.0})
            tr = cobra.Reaction(f"MEDt_{c}", lower_bound=0.0, upper_bound=1000.0)
            tr.add_metabolites({ext: -1.0, cyt: 1.0})
            model.add_reactions([ex, tr])
            added.append(ex_id)
    return sorted(set(added))


def read_gz(path):
    with gzip.open(path, "rt") as fh:
        m = cobra.io.read_sbml_model(fh)
    m.solver = "glpk"
    return m


def main():
    rep = {}
    # (a) + (b): published iGD1575, as read, no translation
    pub = read_model(MODELS["iGD1575"]["orig"])
    rep["iGD1575_published_as_read"] = {"values": egc_values(pub, CURRENCIES_MS)}
    print("published iGD1575:", rep["iGD1575_published_as_read"]["values"], flush=True)
    rep["iGD1575_published_as_read"]["atp_cycles"] = enumerate_cycles(pub, CURRENCIES_MS["atp"])
    c = rep["iGD1575_published_as_read"]["atp_cycles"]
    print("  ATP cycles:", c["n_distinct_minimal_flux_supports"], "single cuts:", [x["id"] for x in c["single_reaction_cuts"]], flush=True)
    for key in ("nadh", "nadph", "q8h2", "h_e (no-periplasm analogue)"):
        if (rep["iGD1575_published_as_read"]["values"].get(key) or 0) > 1e-6:
            rep["iGD1575_published_as_read"][f"{key}_cycles"] = enumerate_cycles(pub, CURRENCIES_MS[key], max_lps=150, depth=1)

    # (c) the BiGG views after medium completion, as the protocol prepares them
    for label in MODELS:
        v = read_model(MODELS[label]["view"])
        added = complete_medium(v, [MODELS[label]["medium"]])
        vals = egc_values(v, CURRENCIES_BIGG)
        rep[f"{label}_view_after_medium_completion"] = {"medium_completion_added": added, "values": vals,
                                                         "gate_keys_testable": [k for k in GATE_KEYS if vals.get(k) is not None],
                                                         "gate_fails": any((vals.get(k) or 0) > 1e-6 for k in GATE_KEYS)}
        print(label, "view:", added, vals, flush=True)
        if (vals.get("h_e (no-periplasm analogue)") or 0) > 1e-6:
            rep[f"{label}_view_after_medium_completion"]["h_e_cycles"] = enumerate_cycles(v, CURRENCIES_BIGG["h_e (no-periplasm analogue)"], max_lps=150, depth=1)

    # draft arms
    for org, medium in (("Smeli", "RCH2_defined_noCarbon"), ("MR1", "ShewMM_noCarbon")):
        for arm in ("B0", "UNQ", "M"):
            m = read_gz(p(f"results/transfer_v1/development/{org}/{arm}/model.xml.gz"))
            added = complete_medium(m, [medium])
            vals = egc_values(m, CURRENCIES_BIGG)
            rep[f"{org}_{arm}_after_medium_completion"] = {"medium_completion_added": added, "values": vals,
                                                            "gate_keys_testable": [k for k in GATE_KEYS if vals.get(k) is not None],
                                                            "gate_fails": any((vals.get(k) or 0) > 1e-6 for k in GATE_KEYS)}
            print(org, arm, added, vals, flush=True)
    dump(rep, "check_energy_gate.json")


if __name__ == "__main__":
    main()
