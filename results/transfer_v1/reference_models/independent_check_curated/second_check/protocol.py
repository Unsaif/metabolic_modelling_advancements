"""This check's own implementation of the transfer protocol's simulation (no gembench import), with the
exploratory medium supplement used for iAH991 as an option.

  * every exchange closed (0, 1000); medium completion: an exchange and a gene-less uptake for each component of
    the medium (and of the supplement, as the repository's patched base_medium does) that the model has in the
    cytosol but cannot exchange (pantothenate, folate and bicarbonate excluded);
  * per condition: the medium at the study's bounds (bulk -1000, trace organics -0.001, trace metals -0.1), the
    supplement, each carbon source at -10; wild-type growth; if >= 1e-3, every listed gene knocked out in turn;
  * infeasible = 0; any other non-optimal status raises.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import fb_conditions, medium_bounds, medium_components, p  # noqa: E402

import cobra  # noqa: E402

ST = 1e-3
EXCLUDE = ("pnto__R", "fol", "hco3")


def fitness_table(org):
    fit = pd.read_table(p(f"data/fitness_browser/{org}/fit_logratios.tsv"), dtype=str, keep_default_na=False)
    meta = [c for c in ["orgId", "locusId", "sysName", "geneName", "desc"] if c in fit.columns]
    cols = [c for c in fit.columns if c not in meta]
    f = fit[cols].apply(pd.to_numeric, errors="coerce")
    f.columns = [c.split(" ")[0] for c in cols]
    f.index = fit["sysName"].where(fit["sysName"] != "", fit["locusId"]).values
    return f


def complete_medium(model, media, supplement=None):
    added = []
    for medium in media:
        comps, _ = medium_components(medium)
        comps = list(comps) + [c for c in (supplement or {}) if c not in comps]
        for c in comps:
            ex_id = f"EX_{c}_e"
            if ex_id in model.reactions or f"{c}_c" not in model.metabolites or c in EXCLUDE:
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
    for ex in model.exchanges:
        ex.bounds = (0.0, 1000.0)
    return sorted(set(added))


def growth(model):
    v = model.slim_optimize(error_value=float("nan"))
    st = model.solver.status
    if st == "infeasible":
        return 0.0
    if st != "optimal" or not np.isfinite(v):
        raise RuntimeError(st)
    return float(v)


def apply_condition(model, medium, carbon, supplement=None, carbon_uptake=-10.0):
    for ex in model.exchanges:
        ex.bounds = (0.0, 1000.0)
    for rid, lb in medium_bounds(medium).items():
        if rid in model.reactions:
            model.reactions.get_by_id(rid).lower_bound = lb
    for c, lb in (supplement or {}).items():
        if f"EX_{c}_e" in model.reactions:
            model.reactions.get_by_id(f"EX_{c}_e").lower_bound = lb
    for b in carbon:
        if f"EX_{b}_e" in model.reactions:
            model.reactions.get_by_id(f"EX_{b}_e").lower_bound = carbon_uptake


def simulate(model, genes, conds, supplement=None, knockouts=True):
    sim = np.zeros((len(genes), len(conds)))
    wt = np.zeros(len(conds))
    for j, c in enumerate(conds):
        with model:
            apply_condition(model, c["media"], c["bigg_ids"], supplement)
            wt[j] = growth(model)
            if knockouts and wt[j] >= ST:
                for i, g in enumerate(genes):
                    with model:
                        model.genes.get_by_id(g).knock_out()
                        sim[i, j] = growth(model)
    return sim, wt


def mapped_conditions(org):
    return [c for c in fb_conditions(org) if c["bigg_ids"]]


def fitness_matrix(org, conds, browser_genes):
    fit = fitness_table(org)
    fm = pd.DataFrame({f"{c['name']} | {c['media']}": fit[c["experiments"]].mean(axis=1, skipna=True) for c in conds})
    return fm.loc[browser_genes].to_numpy()
