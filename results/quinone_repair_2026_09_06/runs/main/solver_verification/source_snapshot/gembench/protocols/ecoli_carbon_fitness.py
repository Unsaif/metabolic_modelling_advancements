"""Protocol: E. coli gene-knockout growth across carbon sources vs RB-TnSeq fitness.

Reproduces Bernstein et al. 2023 (Mol Syst Biol 19:e11566; code at
github.com/dbernste/E_coli_GEM_validation) with every step made explicit and
switchable.  Defaults reproduce the published analysis *as coded* in the
notebook; see `carbon_sources_to_remove` for a documented code/text discrepancy.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import cobra
import numpy as np
from cobra.flux_analysis import single_gene_deletion

from ..datasets import CarbonFitnessDataset
from ..media import M9_MINIMAL_NO_CARBON, Medium, apply_medium, bigg_exchange

# Genes deleted in E. coli BW25113 relative to K-12 MG1655 (araBAD, rhaBAD, lacZ);
# ids as used by Bernstein et al. (looked up in iML1515).
BW25113_DELETED_GENES = ["b0062", "b0063", "b3903", "b3904", "b0061", "b0344", "b3902"]


@dataclass
class ProtocolParams:
    base_medium: Medium = field(default_factory=lambda: M9_MINIMAL_NO_CARBON)
    carbon_uptake: float = -10.0
    growth_threshold: float = 1e-3          # FBA objective below this = no growth (Bernstein sim_thresh)
    fitness_threshold: float = -2.0         # experimental fitness below this = important (Bernstein fit_thresh)
    strain_adjustment: bool = True          # knock out BW25113-deleted genes and drop them from the data
    drop_rich_medium_essentials: bool = True  # drop genes essential with all exchanges open (Bernstein adj_essential)
    rich_medium_uptake: float = -1000.0
    # Bernstein's notebook comment says 'sucrose and mannitol' but the code removes
    # BiGG ids ['man', 'sucr'] where 'man' is D-mannose (mannitol is 'mnl').
    carbon_sources_to_remove: List[str] = field(default_factory=lambda: ["man", "sucr"])
    keep_atpm: bool = True
    processes: int = 2
    solver: str = "glpk"


@dataclass
class ProtocolResult:
    model_id: str
    genes: List[str]
    carbon_sources: List[str]
    sim_growth: np.ndarray            # genes x carbon sources (absolute objective after knockout)
    wt_growth: np.ndarray             # carbon sources
    fitness: np.ndarray               # genes x carbon sources (replicate-averaged experimental fitness)
    missing_medium_components: List[str]
    missing_carbon_exchanges: List[str]
    dropped_strain_genes: List[str]
    dropped_rich_essential_genes: List[str]
    dropped_carbon_sources: List[str]
    counts: Dict[str, int]
    timings_s: Dict[str, float]
    params: ProtocolParams


def run(model: cobra.Model, data: CarbonFitnessDataset, params: Optional[ProtocolParams] = None,
        verbose: bool = True) -> ProtocolResult:
    p = params or ProtocolParams()
    t0 = time.time()
    model = model.copy()
    model.solver = p.solver
    for ex in model.exchanges:            # close everything; media are applied explicitly
        ex.lower_bound = 0.0
        ex.upper_bound = 1000.0

    # 1. gene matching (model gene ids are b-numbers in all four E. coli BiGG models)
    model_gene_ids = {g.id for g in model.genes}
    gene_idx = [i for i, g in enumerate(data.genes) if g in model_gene_ids]
    genes = [data.genes[i] for i in gene_idx]
    fitness = data.fitness[gene_idx, :].copy()
    counts = {"model_genes": len(model_gene_ids), "genes_matched": len(genes),
              "carbon_sources_in_data": len(data.carbon_sources)}

    # 2. strain adjustment
    dropped_strain: List[str] = []
    if p.strain_adjustment:
        for gid in BW25113_DELETED_GENES:
            if gid in model_gene_ids:
                model.genes.get_by_id(gid).knock_out()
        keep = [i for i, g in enumerate(genes) if g not in BW25113_DELETED_GENES]
        dropped_strain = [g for g in genes if g in BW25113_DELETED_GENES]
        genes = [genes[i] for i in keep]; fitness = fitness[keep, :]

    # 3. drop genes essential in an unrestricted (all-uptake) medium
    dropped_rich: List[str] = []
    t1 = time.time()
    if p.drop_rich_medium_essentials:
        with model:
            for ex in model.exchanges:
                ex.lower_bound = p.rich_medium_uptake
                ex.upper_bound = 1000.0
            res = single_gene_deletion(model, genes, processes=p.processes)
            growth = _growth_by_gene(res, genes)
        dropped_rich = [g for g in genes if growth[g] < p.growth_threshold]
        keep = [i for i, g in enumerate(genes) if g not in set(dropped_rich)]
        genes = [genes[i] for i in keep]; fitness = fitness[keep, :]
    timings = {"rich_medium_essentials_s": time.time() - t1}

    # 4. carbon-source filtering
    carbon = list(data.carbon_sources)
    dropped_carbon = [c for c in carbon if c in set(p.carbon_sources_to_remove)]
    keep_c = [j for j, c in enumerate(carbon) if c not in set(p.carbon_sources_to_remove)]
    carbon = [carbon[j] for j in keep_c]; fitness = fitness[:, keep_c]
    counts.update({"genes_after_adjustment": len(genes), "carbon_sources_after_adjustment": len(carbon)})

    # 5. media check
    missing_medium = apply_medium(model, p.base_medium, close_all=True)
    missing_carbon = [c for c in carbon if bigg_exchange(c) not in model.reactions]

    # 6. simulate knockouts per carbon source
    sim = np.zeros((len(genes), len(carbon)))
    wt = np.zeros(len(carbon))
    t2 = time.time()
    for j, c in enumerate(carbon):
        ex_id = bigg_exchange(c)
        with model:
            if ex_id in model.reactions:
                model.reactions.get_by_id(ex_id).lower_bound = p.carbon_uptake
            wt[j] = _nan0(model.slim_optimize())
            res = single_gene_deletion(model, genes, processes=p.processes)
            growth = _growth_by_gene(res, genes)
            sim[:, j] = [growth[g] for g in genes]
        if verbose:
            print(f"  [{model.id}] {j+1:2d}/{len(carbon)} {c:10s} wt={wt[j]:.4f} "
                  f"no-growth KOs={(sim[:, j] < p.growth_threshold).sum():4d}  ({time.time()-t2:.0f}s)", flush=True)
    timings["knockout_simulation_s"] = time.time() - t2
    timings["total_s"] = time.time() - t0

    return ProtocolResult(model_id=model.id, genes=genes, carbon_sources=carbon, sim_growth=sim, wt_growth=wt,
                          fitness=fitness, missing_medium_components=missing_medium,
                          missing_carbon_exchanges=missing_carbon, dropped_strain_genes=dropped_strain,
                          dropped_rich_essential_genes=dropped_rich, dropped_carbon_sources=dropped_carbon,
                          counts=counts, timings_s=timings, params=p)


def _nan0(x: float) -> float:
    return 0.0 if (x is None or np.isnan(x)) else float(x)


def _growth_by_gene(res, genes: List[str]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for ids, g in zip(res["ids"], res["growth"]):
        gid = next(iter(ids))
        out[gid] = _nan0(g)
    for g in genes:
        out.setdefault(g, 0.0)
    return out
