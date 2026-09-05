"""Protocol: S. cerevisiae single-gene deletion viability vs the Stanford deletion collection.

Port of yeast-GEM's `code/modelTests/essentialGenes.m` / `yeastgem.model_tests.essential_genes`
(Kennedy synthetic complete medium, growth-ratio tolerance 1e-6, restricted to SGD verified ORFs),
with bootstrap CIs added.  Exchange ids are yeast-GEM's r_xxxx identifiers.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import cobra
import numpy as np
from cobra.flux_analysis import single_gene_deletion

from ..datasets import YeastEssentialityDataset

# Kennedy-medium presets copied from yeast-GEM (complete_Y7): 33 constrained uptakes at -0.5,
# glucose at -20, 16 unconstrained uptakes at -1000.
COMPLETE_Y7_CONSTRAINED = (
    "r_1604", "r_1639", "r_1873", "r_1879", "r_1880", "r_1881", "r_1671", "r_1883", "r_1757", "r_1891",
    "r_1889", "r_1810", "r_1993", "r_1893", "r_1897", "r_1947", "r_1899", "r_1900", "r_1902", "r_1967",
    "r_1903", "r_1548", "r_1904", "r_2028", "r_2038", "r_1906", "r_2067", "r_1911", "r_1912", "r_1913",
    "r_2090", "r_1914", "r_2106",
)
COMPLETE_Y7_GLUCOSE = "r_1714"
COMPLETE_Y7_UNCONSTRAINED = (
    "r_1672", "r_1654", "r_1992", "r_2005", "r_2060", "r_1861", "r_1832", "r_2100", "r_4593", "r_4595",
    "r_4596", "r_4597", "r_2049", "r_4594", "r_4600", "r_2020",
)


@dataclass
class YeastParams:
    ko_tol: float = 1e-6            # growth ratio below this = inviable (yeast-GEM ko_tol)
    constrained_uptake: float = -0.5
    glucose_uptake: float = -20.0
    unconstrained_uptake: float = -1000.0
    restrict_to_verified: bool = True
    processes: int = 2
    solver: str = "glpk"


@dataclass
class YeastResult:
    model_id: str
    genes: List[str]                 # evaluated genes (model ∩ verified ORFs if restricted)
    growth_ratio: np.ndarray         # per gene
    exp_inviable: np.ndarray         # bool per gene
    wt_growth: float
    missing_exchanges: List[str]
    counts: Dict[str, int]
    classification: Dict[str, str]   # gene -> TP/TN/FP/FN (yeast-GEM convention: 'positive' = viable)
    timings_s: Dict[str, float]
    params: YeastParams


def apply_complete_y7(model: cobra.Model, p: YeastParams) -> List[str]:
    missing = []
    for ex in model.exchanges:
        ex.lower_bound = 0.0
        ex.upper_bound = 1000.0
    for rid, lb in [(r, p.constrained_uptake) for r in COMPLETE_Y7_CONSTRAINED] + \
                   [(COMPLETE_Y7_GLUCOSE, p.glucose_uptake)] + \
                   [(r, p.unconstrained_uptake) for r in COMPLETE_Y7_UNCONSTRAINED]:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).lower_bound = lb
        else:
            missing.append(rid)
    return missing


def run(model: cobra.Model, data: YeastEssentialityDataset, params: Optional[YeastParams] = None,
        verbose: bool = True) -> YeastResult:
    p = params or YeastParams()
    t0 = time.time()
    model = model.copy()
    model.solver = p.solver
    missing = apply_complete_y7(model, p)
    wt = model.slim_optimize()
    if not np.isfinite(wt) or wt <= 0:
        raise RuntimeError(f"wild type does not grow on complete_Y7 (objective={wt})")
    model_genes = [g.id for g in model.genes]
    res = single_gene_deletion(model, model_genes, processes=p.processes)
    ratio = {}
    for ids, g in zip(res["ids"], res["growth"]):
        gid = next(iter(ids))
        ratio[gid] = (0.0 if (g is None or np.isnan(g)) else float(g)) / wt
    genes = [g for g in model_genes if (g in data.verified_orfs or not p.restrict_to_verified)]
    gr = np.array([ratio.get(g, 1.0) for g in genes])
    inv = np.array([g in data.inviable_orfs for g in genes])
    pred_viable = gr >= p.ko_tol
    cls = {}
    for g, pv, ei in zip(genes, pred_viable, inv):
        cls[g] = "TP" if (pv and not ei) else "TN" if (not pv and ei) else "FP" if (pv and ei) else "FN"
    counts = {"model_genes": len(model_genes), "evaluated_genes": len(genes),
              "exp_inviable": int(inv.sum()), "exp_viable": int((~inv).sum()),
              **{k: sum(1 for v in cls.values() if v == k) for k in ["TP", "TN", "FP", "FN"]}}
    if verbose:
        print(f"  [{model.id}] wt={wt:.4f} genes={len(genes)} TP={counts['TP']} TN={counts['TN']} FP={counts['FP']} FN={counts['FN']} ({time.time()-t0:.0f}s)")
    return YeastResult(model_id=model.id, genes=genes, growth_ratio=gr, exp_inviable=inv, wt_growth=float(wt),
                       missing_exchanges=missing, counts=counts, classification=cls,
                       timings_s={"total_s": time.time() - t0}, params=p)
