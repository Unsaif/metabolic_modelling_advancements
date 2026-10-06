"""Check 2: iAH991 translation (rebuilt SBML -> BiGG view) and the identifiers the protocol uses for B. theta.

(a) Relabelling only: reaction by reaction, stoichiometry (mapped back with the renames), bounds, gene rules (as
    Boolean functions) and objective; every reaction, metabolite and gene accounted for; optimal growth with both
    versions in the protocol medium for each mapped condition (strict and with the exploratory supplement), with
    every boundary open, and in 30 random media; 60 single-gene deletions.
(b) Every protocol identifier (Varel_Bryant_medium components, the 25 mapped carbon sources, the energy-check
    currencies) traced to the model metabolite behind it, with the S10b name and formula (from the rebuilt SBML,
    cross-checked against this check's own reading of S10b) next to the BiGG name and formula (CarveMe universe,
    iJN1463, iML1515). Absent identifiers are searched for under other names.
(c) The arabinan override (arabinan101 -> araban__L): what each identifier denotes, in iAH991, in the BiGG universe
    and in the B. theta drafts.
(d) Energy gate on the BiGG view after medium completion (strict and supplemented), with the h_e analogue.
(e) Gene mapping BT_0554 -> BT0554: rates and ten spot checks.

Usage: python -I check_iah991_translation.py <reparse.json>   (writes check_iah991_translation.json)
"""
from __future__ import annotations

import gzip
import json
import math
import os
import random
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
from common import ROOT, norm_name, p, universe_metabolites  # noqa: E402
from check_translation import rule_dnf  # noqa: E402
from check_energy_gate import CURRENCIES_BIGG, egc_values  # noqa: E402
from protocol import apply_condition, complete_medium, growth, mapped_conditions  # noqa: E402
from common import medium_components  # noqa: E402

import cobra  # noqa: E402

ORIG = "models/curated/iAH991/iAH991_rebuilt.xml"
VIEW = "models/curated/iAH991/iAH991_bigg_view.xml.gz"
TRANS = "models/curated/iAH991/iAH991_translation.json"
SUP = {"cbl1": -0.001, "h2s": -1000.0}
TOL = 1e-9


def read(rel):
    path = p(rel)
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            m = cobra.io.read_sbml_model(fh)
    else:
        m = cobra.io.read_sbml_model(path)
    m.solver = "glpk"
    return m


def formula_no_h(f):
    if not f:
        return None
    out = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", f):
        if el != "H":
            out[el] = out.get(el, 0) + (int(n) if n else 1)
    return dict(sorted(out.items()))


def main():
    rp = json.load(open(sys.argv[1]))
    s10b = {x["cells"][1]["despaced"]: {"name": " ".join(x["cells"][2]["lines"]), "formula": x["cells"][3]["despaced"],
                                        "seed": x["cells"][0]["despaced"]}
            for x in rp["263-306"] if "cells" in x and not x["cells"][0]["despaced"].startswith("SEED")}
    orig, view = read(ORIG), read(VIEW)
    tr = json.load(open(p(TRANS)))
    met_back = {new: old for old, new in tr["metabolite_renames"].items()}
    rxn_back = {new: old for old, new in tr["exchange_renames"].items()}
    rep = {"renames": {"metabolites": tr["metabolite_renames"], "exchanges": len(tr["exchange_renames"])}}
    # (a) structure
    mism = {"missing": [], "stoichiometry": [], "bounds": [], "gene_rule": [], "objective": []}
    used = set()
    for vr in view.reactions:
        oid = rxn_back.get(vr.id, vr.id)
        if oid not in orig.reactions:
            mism["missing"].append(vr.id)
            continue
        used.add(oid)
        orr = orig.reactions.get_by_id(oid)
        vst = {met_back.get(m.id, m.id): c for m, c in vr.metabolites.items()}
        ost = {m.id: c for m, c in orr.metabolites.items()}
        if set(vst) != set(ost) or any(abs(vst[k] - ost[k]) > TOL for k in ost):
            mism["stoichiometry"].append(vr.id)
        if any(abs(a - b) > TOL for a, b in zip(vr.bounds, orr.bounds)):
            mism["bounds"].append(vr.id)
        if rule_dnf(vr) != rule_dnf(orr):
            mism["gene_rule"].append(vr.id)
        if abs(vr.objective_coefficient - orr.objective_coefficient) > TOL:
            mism["objective"].append(vr.id)
    rep["structure"] = {"n_reactions": [len(orig.reactions), len(view.reactions)], "n_metabolites": [len(orig.metabolites), len(view.metabolites)],
                        "n_genes": [len(orig.genes), len(view.genes)],
                        "genes_identical": sorted(g.id for g in orig.genes) == sorted(g.id for g in view.genes),
                        "unaccounted_published_reactions": sorted(set(r.id for r in orig.reactions) - used),
                        "mismatch_counts": {k: len(v) for k, v in mism.items()}, "mismatches": mism}
    # growth tests (ids are equal up to the renames; map constraints back)
    tests = []
    vb = [r for r in view.reactions if r.boundary]

    def run(name, cons, ko=None):
        with view, orig:
            for vid, b in cons.items():
                view.reactions.get_by_id(vid).bounds = b
                orig.reactions.get_by_id(rxn_back.get(vid, vid)).bounds = b
            if ko:
                view.genes.get_by_id(ko).knock_out()
                orig.genes.get_by_id(ko).knock_out()
            a = orig.slim_optimize(error_value=float("nan"))
            b = view.slim_optimize(error_value=float("nan"))
        eq = (math.isnan(a) and math.isnan(b)) or abs(a - b) <= 1e-6 * max(1.0, abs(a))
        tests.append({"test": name, "published": a, "view": b, "equal": eq})

    run("published bounds", {})
    run("every boundary open", {r.id: (-1000.0, 1000.0) for r in vb})
    exv = {r.id for r in view.exchanges}
    for sup_name, sup in (("strict", {}), ("supplemented", SUP)):
        for c in mapped_conditions("Btheta"):
            with view:
                apply_condition(view, "Varel_Bryant_medium", c["bigg_ids"], sup)
                cons = {r.id: r.bounds for r in view.reactions if r.boundary}
            run(f"{sup_name}: {c['name']}", cons)
    rng = random.Random(1)
    for k in range(30):
        run(f"random medium {k}", {r.id: ((0.0, 1000.0) if rng.random() < 0.3 else (-rng.choice([1.0, 10.0, 1000.0]), 1000.0)) for r in vb})
    with view:
        apply_condition(view, "Varel_Bryant_medium", ["glc__D"], SUP)
        cons = {r.id: r.bounds for r in view.reactions if r.boundary}
    genes = sorted(g.id for g in view.genes)
    rng.shuffle(genes)
    for gid in genes[:60]:
        run(f"supplemented glucose, knockout {gid}", cons, ko=gid)
    rep["growth_tests"] = {"n": len(tests), "n_equal": sum(t["equal"] for t in tests),
                           "n_positive": sum(1 for t in tests if (t["view"] or 0) > 1e-6), "not_equal": [t for t in tests if not t["equal"]]}
    # (b) identifiers
    uni = universe_metabolites()
    orig_names = {m.id: (m.name, m.formula) for m in orig.metabolites}
    comps, trace = medium_components("Varel_Bryant_medium")
    conds = mapped_conditions("Btheta")
    carbon = sorted({b for c in conds for b in c["bigg_ids"]})
    currencies = ["atp", "adp", "pi", "h", "h2o", "nadh", "nad", "nadph", "nadp", "q8h2", "q8"]
    rows = []
    for role, ids, comp in (("medium", comps, "e"), ("medium (cytosol, for completion)", comps, "c"), ("carbon", carbon, "e"),
                            ("currency", currencies, "c"), ("supplement", list(SUP), "e")):
        for b in ids:
            vid = f"{b}_{comp}"
            row = {"role": role, "bigg_id": b, "view_metabolite": vid if vid in view.metabolites else "",
                   "exchange_in_view": (f"EX_{b}_e" in view.reactions) if comp == "e" else None}
            if vid in view.metabolites:
                pid = met_back.get(vid, vid)
                base = pid.rsplit("_", 1)[0]
                printed = base.replace("__", "-")
                name, formula = orig_names.get(pid, ("", ""))
                u = uni.get(b)
                s = s10b.get(printed, {})
                row.update({"published_id": pid, "s10b_id": printed, "s10b_name_own_reading": s.get("name"),
                            "s10b_formula_own_reading": s.get("formula"), "model_name": name, "model_formula": formula,
                            "bigg_name": u[0] if u else None, "bigg_formula": u[1] if u else None,
                            "name_match": bool(u) and norm_name(u[0]) == norm_name(name),
                            "formula_match_ignoring_H": bool(u) and bool(formula) and formula_no_h(u[1]) == formula_no_h(formula),
                            "renamed": pid != vid})
            rows.append(row)
    rep["protocol_identifiers"] = rows
    # absent identifiers: search names
    absent = [r for r in rows if not r["view_metabolite"] and r["role"] in ("medium", "carbon", "currency", "supplement")]
    search = {}
    for r in absent:
        b = r["bigg_id"]
        u = uni.get(b)
        targets = {norm_name(u[0])} if u else set()
        hits = [{"id": m.id, "name": m.name, "formula": m.formula, "has_exchange": any(x.boundary for x in m.reactions)}
                for m in view.metabolites if norm_name(m.name) in targets or m.id.split("_")[0] == b.split("__")[0]]
        search[b] = {"bigg_name": u[0] if u else None, "hits": hits}
    rep["absent_search"] = search
    # (c) arabinan
    arab = {}
    for mid in ("araban__L_e", "arabinan101_e"):
        if mid in view.metabolites:
            m = view.metabolites.get_by_id(mid)
            arab[f"view:{mid}"] = {"name": m.name, "formula": m.formula, "reactions": [f"{r.id}: {r.reaction}" for r in m.reactions]}
    s = s10b.get("arabinan101")
    arab["s10b:arabinan101"] = s
    u = uni.get("araban__L")
    arab["bigg_universe:araban__L"] = {"name": u[0], "formula": u[1]} if u else None
    for arm in ("B0", "UNQ", "M"):
        with gzip.open(p(f"results/transfer_v1/development/Btheta/{arm}/model.xml.gz"), "rt") as fh:
            dm = cobra.io.read_sbml_model(fh)
        if "araban__L_e" in dm.metabolites:
            m = dm.metabolites.get_by_id("araban__L_e")
            arab[f"draft {arm}:araban__L_e"] = {"name": m.name, "formula": m.formula, "reactions": [f"{r.id}: {r.reaction}" for r in m.reactions][:6]}
        for extra in ("amylose300_e", "starch1200_e", "lmn2_e"):
            if extra in dm.metabolites:
                arab.setdefault(f"draft {arm}: other non-BiGG carbon ids present", []).append(extra)
    ct = pd.read_table(p("data/reference/fitness_browser_carbon_sources_bigg.tsv"), dtype=str, keep_default_na=False)
    arab["carbon_table_row"] = ct[ct["bigg_ids"].str.contains("araban")].to_dict("records")
    rep["arabinan"] = arab
    # (d) gate
    gate = {}
    for name, sup in (("strict", {}), ("supplemented", SUP)):
        v = read(VIEW)
        added = complete_medium(v, ["Varel_Bryant_medium"], sup)
        gate[name] = {"completion_added": added, "values": egc_values(v, CURRENCIES_BIGG)}
    rep["energy_gate"] = gate
    # (e) genes
    g = pd.read_table(p("data/fitness_browser/Btheta/genes.tsv"), dtype=str, keep_default_na=False)
    sysn = set(g["sysName"])
    desc = dict(zip(g["sysName"], g["desc"]))
    fit = pd.read_table(p("data/fitness_browser/Btheta/fit_logratios.tsv"), dtype=str, keep_default_na=False, usecols=["locusId", "sysName"])
    fitset = set(fit["sysName"].where(fit["sysName"] != "", fit["locusId"]))
    mp = {x.id: x.id.replace("_", "") for x in view.genes if x.id.replace("_", "") in sysn}
    rep["gene_mapping"] = {"model_genes": len(view.genes), "mapped": len(mp), "with_fitness": sum(1 for v in mp.values() if v in fitset),
                           "unmapped": sorted(x.id for x in view.genes if x.id not in mp)}
    rng = random.Random(20261006)
    pick = sorted(k for k, v in mp.items() if v in fitset)
    rng.shuffle(pick)
    rep["gene_spot_checks"] = [{"gene": k, "sysName": mp[k], "desc": desc.get(mp[k]),
                                "reactions": [f"{r.id} ({r.name})" for r in view.genes.get_by_id(k).reactions][:6]} for k in pick[:10]]
    json.dump(rep, open(os.path.join(HERE, "check_iah991_translation.json"), "w"), indent=1, default=float)
    print("structure", json.dumps(rep["structure"]["mismatch_counts"]), rep["structure"]["n_reactions"], rep["structure"]["genes_identical"],
          rep["structure"]["unaccounted_published_reactions"])
    print("growth tests", {k: v for k, v in rep["growth_tests"].items() if k != "not_equal"})
    for r in rows:
        if r["view_metabolite"]:
            flag = "" if (r.get("name_match") or r.get("formula_match_ignoring_H")) else "  <-- REVIEW"
            print(f"  {r['role'][:10]:10s} {r['bigg_id']:10s} {r['published_id']:16s} | {str(r['model_name'])[:38]:38s} {str(r['model_formula'])[:16]:16s} | "
                  f"{str(r['bigg_name'])[:30]:30s} {str(r['bigg_formula'])[:16]:16s}{flag}")
        else:
            print(f"  {r['role'][:10]:10s} {r['bigg_id']:10s} ABSENT")
    print("absent search", json.dumps(search)[:1500])
    print("arabinan", json.dumps(arab, default=str)[:3000])
    print("gate", json.dumps(gate))
    print("genes", {k: v for k, v in rep["gene_mapping"].items()})
    for s in rep["gene_spot_checks"]:
        print("  ", s["gene"], "|", s["desc"], "|", s["reactions"][:3])


if __name__ == "__main__":
    main()
