"""Check 2: is every identifier the protocol uses mapped to the same compound?

For iSO783 (MR1, ShewMM_noCarbon) and iGD1575 (Smeli, RCH2_defined_noCarbon), every BiGG identifier the protocol
looks up -- medium components, the carbon sources of the mapped conditions, the energy-check currencies -- is
traced to the published metabolite behind it, and compared on the evidence available:

  * the published name (iGD1575; the BioModels iSO783 file has no names or formulas, only identifiers);
  * the ModelSEED compound behind the published identifier: iGD1575 identifiers are ModelSEED compounds; iSO783
    identifiers are placed through the alias table's own "iSO783" source (ModelSEED's mapping of this very model),
    which the translation script did not use;
  * the ModelSEED compound's name, formula and BiGG aliases (ModelSEED database 194ac8a, compounds.tsv);
  * the BiGG name and formula (CarveMe universe, iJN1463, iML1515).

For every identifier reported absent, the published model is searched for the compound under any identifier
(ModelSEED compound with that BiGG alias; same name; same compound in another compartment).

The same comparison is run over every renamed metabolite of both models as a broader check of the translation.

Usage: python -I check_identifier_mapping.py <ModelSEED compounds.tsv at 194ac8a>
Writes check_identifier_mapping.json, identifier_mapping_protocol.tsv and modelseed_compounds_subset.tsv.
"""
from __future__ import annotations

import collections
import hashlib
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (MODELS, OUT, back_maps, carbon_table, dump, fb_conditions, medium_components, modelseed_aliases,  # noqa: E402
                    norm_name, read_model, translation, universe_metabolites)

CURRENCY = ["atp", "adp", "pi", "h", "h2o", "nadh", "nad", "nadph", "nadp", "q8h2", "q8"]


def formula_no_h(f):
    """Element counts without hydrogen (protonation states differ between namespaces)."""
    if not f or f in ("null", "None"):
        return None
    counts = collections.Counter()
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", f):
        if el != "H":
            counts[el] += int(n) if n else 1
    return dict(sorted(counts.items()))


def load_compounds(path):
    cp = pd.read_table(path, dtype=str, keep_default_na=False)
    cp = cp[["id", "abbreviation", "name", "formula", "charge", "aliases"]]
    return cp.set_index("id")


def ms_names(row):
    names = {row["name"], row["abbreviation"]}
    for part in row["aliases"].split("|"):
        if part.startswith("Name:"):
            names |= {x.strip() for x in part[5:].split(";")}
    return {n for n in names if n}


def ms_bigg(row):
    for part in row["aliases"].split("|"):
        if part.startswith("BiGG:"):
            return {x.strip() for x in part[5:].split(";")}
    return set()


class Evidence:
    def __init__(self, compounds_path):
        self.cp = load_compounds(compounds_path)
        al = modelseed_aliases()
        self.al = al
        self.iso_to_cpd = collections.defaultdict(set)
        for c, e in zip(al.loc[al.source == "iSO783", "cpd"], al.loc[al.source == "iSO783", "ext"]):
            self.iso_to_cpd[e].add(c)
        self.bigg_to_cpd = collections.defaultdict(set)
        for c, e in zip(al.loc[al.source == "BiGG", "cpd"], al.loc[al.source == "BiGG", "ext"]):
            self.bigg_to_cpd[e].add(c)
        self.uni = universe_metabolites()

    def cpds_for_published(self, label, base):
        if label == "iGD1575":
            return {base} if base.startswith("cpd") else set()
        # iSO783: SBML ids encode '-' as '__' and '(', ')' as _LPAREN_/_RPAREN_; the alias table uses the raw ids
        raw = base.replace("_LPAREN_", "(").replace("_RPAREN_", ")").replace("_COMMA_", ",").replace("__", "-")
        cands = {base, raw, raw.replace("-", "_")}
        if base == "o2__":       # superoxide, 'o2-' in the original
            cands.add("o2-")
        out = set()
        for c in cands:
            out |= self.iso_to_cpd.get(c, set())
        return out

    def cpd_info(self, cpd):
        if cpd not in self.cp.index:
            return {"cpd": cpd, "missing_from_compounds_tsv": True}
        r = self.cp.loc[cpd]
        return {"cpd": cpd, "name": r["name"], "formula": r["formula"], "charge": r["charge"], "bigg_aliases": sorted(ms_bigg(r))}

    def verdict(self, label, pub_base, pub_name, bigg_id):
        """Is the published metabolite pub_base (name pub_name) the BiGG compound bigg_id?"""
        cpds = self.cpds_for_published(label, pub_base)
        infos = [self.cpd_info(c) for c in sorted(cpds)]
        uni = self.uni.get(bigg_id)
        alias_ok = any(bigg_id in i.get("bigg_aliases", []) for i in infos)
        name_ok = False
        if uni:
            un = norm_name(uni[0])
            names = {norm_name(pub_name)} if pub_name else set()
            for c in cpds:
                if c in self.cp.index:
                    names |= {norm_name(n) for n in ms_names(self.cp.loc[c])}
            name_ok = un in names
        formula_ok = None
        if uni and uni[1]:
            fu = formula_no_h(uni[1])
            fms = [formula_no_h(i.get("formula")) for i in infos if i.get("formula")]
            formula_ok = any(f == fu for f in fms) if fms else None
        if alias_ok or name_ok:
            v = "same compound"
        elif formula_ok:
            v = "formula agrees, alias and name do not: review"
        elif not cpds and pub_base == bigg_id:
            v = "identical identifier, no ModelSEED placement of the published id"
        else:
            v = "REVIEW"
        return {"published_id_base": pub_base, "published_name": pub_name, "modelseed": infos,
                "bigg": {"id": bigg_id, "name": uni[0] if uni else None, "formula": uni[1] if uni else None},
                "bigg_alias_of_modelseed_compound": alias_ok, "name_match": name_ok, "formula_match_ignoring_H": formula_ok,
                "verdict": v}


def split(mid):
    m = re.match(r"^(.*)_([a-z]\d?)$", mid)
    return (m.group(1), m.group(2)) if m else (mid, "")


def protocol_ids(label):
    cfg = MODELS[label]
    comps, trace = medium_components(cfg["medium"])
    conds = [c for c in fb_conditions(cfg["org"]) if c["bigg_ids"]]
    carbon = sorted({b for c in conds for b in c["bigg_ids"]})
    return comps, carbon, conds


def search_absent(label, orig, ev, bigg_id, extra_names=()):
    """Look for the compound behind bigg_id anywhere in the published model (also by the Fitness Browser condition
    name, for identifiers that are not BiGG identifiers at all)."""
    hits = []
    target_cpds = ev.bigg_to_cpd.get(bigg_id, set())
    uni = ev.uni.get(bigg_id)
    target_names = set()
    for c in target_cpds:
        if c in ev.cp.index:
            target_names |= {norm_name(n) for n in ms_names(ev.cp.loc[c])}
    if uni:
        target_names.add(norm_name(uni[0]))
    for n in extra_names:
        target_names.add(norm_name(re.sub(r"\(\s*[-+]\s*\)", "", n)))
    target_names.discard("")
    for m in orig.metabolites:
        base, comp = split(m.id)
        cpds = ev.cpds_for_published(label, base)
        why = []
        if cpds & target_cpds:
            why.append("ModelSEED compound has this BiGG alias")
        if norm_name(m.name) in target_names:
            why.append("name")
        for c in cpds:
            if c in ev.cp.index and norm_name(ev.cp.loc[c]["name"]) in target_names:
                why.append("ModelSEED name of the published id")
        if label == "iSO783" and base.replace("__", "-") == bigg_id.replace("__", "-"):
            why.append("same identifier stem")
        if why:
            hits.append({"published_id": m.id, "name": m.name, "compartment": m.compartment, "why": sorted(set(why)),
                         "reactions": [r.id for r in m.reactions][:12], "n_reactions": len(m.reactions),
                         "has_boundary_reaction": any(r.boundary for r in m.reactions)})
    return {"target_modelseed_compounds": sorted(target_cpds), "bigg_name": uni[0] if uni else None,
            "bigg_id_in_universe": bigg_id in ev.uni, "hits_in_published_model": hits}


def main():
    compounds_path = sys.argv[1]
    sha = hashlib.sha256(open(compounds_path, "rb").read()).hexdigest()
    ev = Evidence(compounds_path)
    report = {"modelseed_compounds_tsv": {"url": "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/194ac8afe48f8a606c0dd07ba3c7af10c02ba2fd/Biochemistry/compounds.tsv",
                                          "sha256": sha}}
    rows = []
    used_cpds = set()
    for label, cfg in MODELS.items():
        orig, view = read_model(cfg["orig"]), read_model(cfg["view"])
        met_back, rxn_back = back_maps(label)
        names = {m.id: m.name for m in orig.metabolites}
        comps, carbon, conds = protocol_ids(label)
        rep = {"medium": cfg["medium"], "n_conditions_mapped": len(conds), "items": [], "absent": {}}
        items = [("medium", c) for c in comps] + [("carbon", c) for c in carbon] + [("currency", c) for c in CURRENCY]
        ex_ids = {r.id for r in view.reactions}
        for role, b in items:
            for comp in (("e",) if role != "currency" else ("c",)) + (("c",) if role == "medium" else ()):
                vid = f"{b}_{comp}"
                present = vid in view.metabolites
                has_ex = f"EX_{b}_e" in ex_ids if comp == "e" else None
                row = {"model": label, "role": role, "bigg_id": b, "compartment": comp, "view_metabolite": vid if present else "",
                       "exchange_in_view": has_ex}
                if present:
                    pid = met_back.get(vid, vid)
                    base, _ = split(pid)
                    v = ev.verdict(label, base, names.get(pid, ""), b)
                    used_cpds |= {i["cpd"] for i in v["modelseed"]}
                    row.update({"published_id": pid, "published_name": names.get(pid, ""), "verdict": v["verdict"],
                                "modelseed": ";".join(f"{i['cpd']}={i.get('name')}" for i in v["modelseed"]),
                                "bigg_name": v["bigg"]["name"]})
                    rep["items"].append({"role": role, "compartment": comp, **v})
                else:
                    row.update({"published_id": "", "published_name": "", "verdict": "absent", "modelseed": "", "bigg_name": ""})
                    if comp == "e" or role == "currency":
                        key = f"{b}_{comp}"
                        cnames = [c["name"] for c in conds if b in c["bigg_ids"] and len(c["bigg_ids"]) == 1]
                        rep["absent"][key] = search_absent(label, orig, ev, b, cnames)
                        rep["absent"][key]["condition_names_searched"] = cnames
                        used_cpds |= set(rep["absent"][key]["target_modelseed_compounds"])
                rows.append(row)
        # conditions whose carbon source(s) have no exchange
        rep["conditions_with_absent_exchange"] = [{"condition": c["name"], "missing": [x for x in c["bigg_ids"] if f"EX_{x}_e" not in ex_ids]}
                                                  for c in conds if any(f"EX_{x}_e" not in ex_ids for x in c["bigg_ids"])]
        # broader check: every renamed metabolite
        broad = collections.Counter()
        review = []
        for old, new in translation(label)["metabolite_renames"].items():
            nb, _ = split(new)
            ob, _ = split(old)
            v = ev.verdict(label, ob, names.get(old, ""), nb)
            used_cpds |= {i["cpd"] for i in v["modelseed"]}
            broad[v["verdict"]] += 1
            if v["verdict"] != "same compound":
                review.append({"published": old, "published_name": names.get(old, ""), "view": new, "verdict": v["verdict"],
                               "modelseed": [(i["cpd"], i.get("name"), i.get("formula")) for i in v["modelseed"]],
                               "bigg": (v["bigg"]["name"], v["bigg"]["formula"])})
        rep["all_renames"] = {"counts": dict(broad), "not_confirmed": review}
        # iSO783: identifiers kept as they are, checked against ModelSEED's own iSO783 placement
        if label == "iSO783":
            kept = collections.Counter()
            disagree = []
            for m in view.metabolites:
                if m.id in met_back:
                    continue
                base, _ = split(m.id)
                cpds = ev.cpds_for_published(label, base)
                if not cpds:
                    kept["no ModelSEED placement"] += 1
                    continue
                aliases = set().union(*[set(ev.cpd_info(c).get("bigg_aliases", [])) for c in cpds])
                if base in aliases:
                    kept["kept id is a BiGG alias of ModelSEED's compound for it"] += 1
                else:
                    kept["kept id is NOT a BiGG alias of ModelSEED's compound"] += 1
                    disagree.append({"id": m.id, "modelseed": [ev.cpd_info(c) for c in sorted(cpds)]})
                used_cpds |= cpds
            rep["identifiers_kept"] = {"counts": dict(kept), "disagreements": disagree}
        report[label] = rep
        print(label, "absent:", sorted(rep["absent"]), "| verdicts:",
              dict(collections.Counter(r["verdict"] for r in rows if r["model"] == label)), "| renames:", dict(broad), flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "identifier_mapping_protocol.tsv"), sep="\t", index=False)
    sub = ev.cp.loc[sorted(c for c in used_cpds if c in ev.cp.index)].reset_index()
    sub["aliases"] = sub["aliases"].map(lambda a: "|".join(x for x in a.split("|") if x.startswith(("Name:", "BiGG:", "KEGG:"))))
    sub.to_csv(os.path.join(OUT, "modelseed_compounds_subset.tsv"), sep="\t", index=False)
    dump(report, "check_identifier_mapping.json")


if __name__ == "__main__":
    main()
