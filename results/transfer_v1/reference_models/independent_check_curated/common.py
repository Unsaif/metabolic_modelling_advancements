"""Shared helpers for the independent check of the curated reference models (iSO783, iGD1575).

Written for this check only. Nothing here imports gembench or the repository's scripts: the Fitness Browser
condition selection, the media and the energy check are re-implemented from their descriptions so that a fault in
the repository code cannot hide itself.
"""
from __future__ import annotations

import gzip
import json
import logging
import os
import re
import warnings

import pandas as pd

ROOT = "/home/claude/mma"
OUT = os.path.join(ROOT, "results", "transfer_v1", "reference_models", "independent_check_curated")
logging.getLogger("cobra").setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

import cobra  # noqa: E402

MODELS = {
    "iSO783": {"orig": "models/curated/iSO783/MODEL1507180036_url.xml",
               "view": "models/curated/iSO783/iSO783_bigg_view.xml.gz",
               "trans": "models/curated/iSO783/iSO783_translation.json",
               "org": "MR1", "medium": "ShewMM_noCarbon"},
    "iGD1575": {"orig": "models/curated/iGD1575/ncomms12219-s7.xml",
                "view": "models/curated/iGD1575/iGD1575_bigg_view.xml.gz",
                "trans": "models/curated/iGD1575/iGD1575_translation.json",
                "org": "Smeli", "medium": "RCH2_defined_noCarbon"},
}


def p(rel):
    return os.path.join(ROOT, rel)


def read_model(rel):
    path = p(rel)
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            m = cobra.io.read_sbml_model(fh)
    else:
        m = cobra.io.read_sbml_model(path)
    m.solver = "glpk"
    return m


def translation(label):
    return json.load(open(p(MODELS[label]["trans"])))


def back_maps(label):
    """view id -> published id, for metabolites and reactions."""
    t = translation(label)
    met_back = {new: old for old, new in t["metabolite_renames"].items()}
    rxn_back = {new: old for old, new in t["exchange_renames"].items()}
    if len(met_back) != len(t["metabolite_renames"]) or len(rxn_back) != len(t["exchange_renames"]):
        raise SystemExit(f"{label}: renames are not one-to-one")
    return met_back, rxn_back


# ---------------------------------------------------------------------------------------------------------------
# Study tables and Fitness Browser conditions (own implementation of the selection described in
# gembench/fitness_browser.py: carbon-source group, condition_2 empty or DMSO, medium in the media table,
# experiment present in the fitness table; grouped by condition_1 x media; mapped when the carbon table lists it).

def media_table():
    return pd.read_table(p("data/reference/fitness_browser_media_bigg.tsv"), dtype=str, keep_default_na=False)


def carbon_table():
    return pd.read_table(p("data/reference/fitness_browser_carbon_sources_bigg.tsv"), dtype=str, keep_default_na=False)


def medium_components(medium):
    mm = media_table()
    row = mm[mm["media"] == medium].iloc[0]
    comps = [c for c in row["bigg_ids"].split(";") if c]
    if row["aerobic"].strip().lower() in ("yes", "true", "1") and "o2" not in comps:
        comps.append("o2")
    trace = set(row["trace_components"].split(";")) - {""}
    return comps, trace


TRACE_METALS = ("fe2", "fe3", "mn2", "zn2", "cu2", "cobalt2", "mobd", "ni2", "sel", "slnt", "tungs")


def medium_bounds(medium):
    """{EX_<id>_e: lower bound} as the study defines them (bulk -1000, trace organics -0.001, trace metals -0.1)."""
    comps, trace = medium_components(medium)
    out = {}
    for c in comps:
        out[f"EX_{c}_e"] = -0.001 if c in trace else (-0.1 if c in TRACE_METALS else -1000.0)
    return out


def fb_fitness_columns(org):
    with open(p(f"data/fitness_browser/{org}/fit_logratios.tsv")) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    meta = {"orgId", "locusId", "sysName", "geneName", "desc"}
    return [c.split(" ")[0] for c in header if c not in meta]


def fb_conditions(org):
    exp = pd.read_table(p(f"data/fitness_browser/{org}/experiments.tsv"), dtype=str, keep_default_na=False)
    cols = set(fb_fitness_columns(org))
    known_media = set(media_table()["media"])
    ct = carbon_table()
    if ct["condition"].duplicated().any():
        raise SystemExit("duplicate condition in the carbon table")
    cmap = dict(zip(ct["condition"], ct["bigg_ids"]))
    sel = exp[(exp["expGroup"] == "carbon source") & (exp["condition_2"].isin(["", "Dimethyl Sulfoxide"]))]
    out = {}
    for _, r in sel.iterrows():
        if r["media"] not in known_media or r["expName"] not in cols:
            continue
        key = f"{r['condition_1']} | {r['media']}"
        if key not in out:
            ids = [x for x in cmap.get(r["condition_1"], "").split(";") if x]
            out[key] = {"name": r["condition_1"], "media": r["media"], "bigg_ids": ids, "experiments": []}
        out[key]["experiments"].append(r["expName"])
    return list(out.values())


# ---------------------------------------------------------------------------------------------------------------
# Reference data for compound identity

def modelseed_aliases():
    al = pd.read_table(p("data/reference/namespace/Unique_ModelSEED_Compound_Aliases.txt"), dtype=str, keep_default_na=False)
    al.columns = ["cpd", "ext", "source"]
    return al


def universe_metabolites():
    """BiGG base id -> (name, formula, charge) from the CarveMe universe and the two local BiGG models."""
    out = {}
    for rel in ("results/ppnp_repair_2026_09_06/evidence/source/universe_bacteria.xml.gz", "models/bigg/iJN1463.xml",
                "models/iML1515/iML1515.xml"):
        path = p(rel)
        if not os.path.exists(path):
            continue
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fh:
            text = fh.read()
        for m in re.finditer(r"<species\b([^>]*)>", text):
            attrs = dict(re.findall(r'([\w:]+)="([^"]*)"', m.group(1)))
            sid = attrs.get("id", "")
            mm = re.match(r"^M_(.+)_([a-z]{1,2})$", sid)
            if not mm:
                continue
            base = mm.group(1)
            if base not in out:
                out[base] = (attrs.get("name", ""), attrs.get("fbc:chemicalFormula", ""), attrs.get("fbc:charge", ""))
    return out


def norm_name(s):
    s = (s or "").lower()
    s = re.sub(r"_(c0|e0|b|c|e)$", "", s)
    return re.sub(r"[^a-z0-9]", "", s)


def dump(obj, name):
    path = os.path.join(OUT, name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=float)
        fh.write("\n")
    return path
