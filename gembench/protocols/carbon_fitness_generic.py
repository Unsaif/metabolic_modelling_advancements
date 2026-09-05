"""Protocol: gene-knockout growth across carbon sources vs RB-TnSeq fitness, for any organism
with Fitness Browser data and a BiGG-namespace model.

Generalises the Bernstein et al. 2023 E. coli protocol: for every carbon-source condition the base
medium is the Fitness Browser medium of that experiment (mapped to BiGG components), the carbon
source(s) are opened at `carbon_uptake`, wild-type growth is computed, and every model gene with
fitness data is deleted in turn.  Conditions on which the wild-type model does not grow are kept
in the result (they are a model prediction in their own right — the organism did grow, or there
would be no fitness data) but excluded from gene-level scoring by default.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import cobra
import numpy as np
import pandas as pd
from cobra.flux_analysis import single_gene_deletion

from ..fitness_browser import Condition, FitnessBrowserOrganism, base_medium, replicate_averaged_fitness
from ..gene_mapping import GeneMap
from ..media import apply_medium


@dataclass
class GenericParams:
    carbon_uptake: float = -10.0
    growth_threshold: float = 1e-3          # absolute objective value below this = no growth
    fitness_threshold: float = -2.0
    drop_rich_medium_essentials: bool = True
    rich_medium_uptake: float = -1000.0
    knockout_genes: List[str] = field(default_factory=list)   # model genes absent from the assayed strain
    processes: int = 2
    solver: str = "glpk"
    max_conditions: Optional[int] = None    # for smoke tests


@dataclass
class GenericResult:
    model_id: str
    org_id: str
    model_genes: List[str]              # model gene ids simulated
    browser_genes: List[str]            # matching Fitness Browser sysNames (same order)
    conditions: List[Condition]
    sim_growth: np.ndarray              # genes x conditions
    wt_growth: np.ndarray               # conditions
    fitness: np.ndarray                 # genes x conditions (replicate-averaged)
    missing_medium_components: Dict[str, List[str]]
    missing_carbon_exchanges: Dict[str, List[str]]
    dropped_rich_essential_genes: List[str]
    counts: Dict[str, int]
    timings_s: Dict[str, float]
    params: GenericParams

    def condition_table(self) -> pd.DataFrame:
        rows = []
        for j, c in enumerate(self.conditions):
            rows.append({"condition": c.name, "media": c.media, "bigg_ids": ";".join(c.bigg_ids),
                         "mapping_confidence": c.mapping_confidence, "n_experiments": len(c.experiments),
                         "wt_growth": float(self.wt_growth[j]),
                         "wt_grows": bool(self.wt_growth[j] >= self.params.growth_threshold),
                         "missing_exchanges": ";".join(self.missing_carbon_exchanges.get(c.key, [])),
                         "n_ko_no_growth": int((self.sim_growth[:, j] < self.params.growth_threshold).sum()),
                         "n_genes_fitness_below_threshold": int(np.nansum(self.fitness[:, j] < self.params.fitness_threshold))})
        return pd.DataFrame(rows)


def run(model: cobra.Model, org: FitnessBrowserOrganism, conditions: List[Condition], gene_map: GeneMap,
        params: Optional[GenericParams] = None, verbose: bool = True) -> GenericResult:
    p = params or GenericParams()
    t0 = time.time()
    model = model.copy()
    model.solver = p.solver
    for ex in model.exchanges:
        ex.lower_bound = 0.0
        ex.upper_bound = 1000.0

    # 1. genes: model genes that map to a Browser gene with fitness data
    fit_index = set(org.fitness.index)
    pairs = [(g.id, gene_map.model_to_browser[g.id]) for g in model.genes
             if g.id in gene_map.model_to_browser and gene_map.model_to_browser[g.id] in fit_index]
    # several model genes can map to one Browser gene (rare); keep the first
    seen, uniq = set(), []
    for mg, bg in pairs:
        if bg not in seen:
            uniq.append((mg, bg)); seen.add(bg)
    model_genes = [mg for mg, _ in uniq]
    browser_genes = [bg for _, bg in uniq]
    counts = {"model_genes": len(model.genes), "model_genes_mapped": len(gene_map.model_to_browser),
              "genes_with_fitness": len(model_genes)}

    # 2. strain adjustment (genes absent from the assayed strain)
    for gid in p.knockout_genes:
        if gid in model.genes:
            model.genes.get_by_id(gid).knock_out()

    # 3. rich-medium essentials
    dropped_rich: List[str] = []
    t1 = time.time()
    if p.drop_rich_medium_essentials and model_genes:
        with model:
            for ex in model.exchanges:
                ex.lower_bound = p.rich_medium_uptake
            res = single_gene_deletion(model, model_genes, processes=p.processes)
            growth = _growth_by_gene(res, model_genes)
        dropped_rich = [g for g in model_genes if growth[g] < p.growth_threshold]
        keep = [i for i, g in enumerate(model_genes) if g not in set(dropped_rich)]
        model_genes = [model_genes[i] for i in keep]; browser_genes = [browser_genes[i] for i in keep]
    timings = {"rich_medium_essentials_s": time.time() - t1}
    counts["genes_after_adjustment"] = len(model_genes)

    # 4. conditions
    conds = [c for c in conditions if c.bigg_ids]
    if p.max_conditions:
        conds = conds[: p.max_conditions]
    counts["conditions_total"] = len(conditions)
    counts["conditions_mapped"] = len(conds)
    fit_df = replicate_averaged_fitness(org, conds)
    fitness = fit_df.loc[browser_genes, [c.key for c in conds]].to_numpy(dtype=float)

    # 5. simulate
    sim = np.zeros((len(model_genes), len(conds)))
    wt = np.zeros(len(conds))
    missing_medium: Dict[str, List[str]] = {}
    missing_carbon: Dict[str, List[str]] = {}
    media_cache: Dict[str, object] = {}
    t2 = time.time()
    for j, c in enumerate(conds):
        if c.media not in media_cache:
            media_cache[c.media] = base_medium(c.media)
        with model:
            miss = apply_medium(model, media_cache[c.media], close_all=True)
            missing_medium[c.media] = miss
            absent = []
            for ex_id in c.exchanges:
                if ex_id in model.reactions:
                    model.reactions.get_by_id(ex_id).lower_bound = p.carbon_uptake
                else:
                    absent.append(ex_id)
            missing_carbon[c.key] = absent
            wt[j] = _nan0(model.slim_optimize())
            if wt[j] >= p.growth_threshold and model_genes:
                res = single_gene_deletion(model, model_genes, processes=p.processes)
                growth = _growth_by_gene(res, model_genes)
                sim[:, j] = [growth[g] for g in model_genes]
            else:
                sim[:, j] = 0.0
        if verbose:
            print(f"  [{org.org_id}/{model.id}] {j+1:2d}/{len(conds)} {c.name[:32]:32s} on {c.media[:22]:22s} "
                  f"wt={wt[j]:.4f} KO-no-growth={(sim[:, j] < p.growth_threshold).sum():4d} "
                  f"absent={absent} ({time.time()-t2:.0f}s)", flush=True)
    timings["knockout_simulation_s"] = time.time() - t2
    timings["total_s"] = time.time() - t0
    counts["conditions_wt_grows"] = int((wt >= p.growth_threshold).sum())

    return GenericResult(model_id=model.id, org_id=org.org_id, model_genes=model_genes, browser_genes=browser_genes,
                         conditions=conds, sim_growth=sim, wt_growth=wt, fitness=fitness,
                         missing_medium_components=missing_medium, missing_carbon_exchanges=missing_carbon,
                         dropped_rich_essential_genes=dropped_rich, counts=counts, timings_s=timings, params=p)


def _nan0(x) -> float:
    return 0.0 if (x is None or (isinstance(x, float) and np.isnan(x))) else float(x)


def _growth_by_gene(res, genes: List[str]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for ids, g in zip(res["ids"], res["growth"]):
        out[next(iter(ids))] = _nan0(g)
    for g in genes:
        out.setdefault(g, 0.0)
    return out
