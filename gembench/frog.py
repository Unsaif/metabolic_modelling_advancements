"""FROG-style reproducibility reports and cross-solver / cross-tool comparison.

A FROG report (BioModels FBC curation standard) records, for a model as shipped:
  F  — FBA objective value
  R  — reaction deletion objective values
  O  — objective (already covered by F)
  G  — gene deletion objective values
plus flux variability (FVA) at a stated fraction of the optimum.

This module computes the same quantities with COBRApy under a chosen solver and
writes them as plain TSV/JSON so that results from other tools (the MATLAB COBRA
Toolbox, COBREXA.jl, cobrar) can be dropped into the same fixture format and
compared.  Only the *as-shipped* bounds are used — no medium is applied — because
the question here is tool/solver agreement, not biology.
"""
from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import cobra
import numpy as np
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis, single_gene_deletion, single_reaction_deletion


@dataclass
class FrogReport:
    model_id: str
    solver: str
    objective_value: float
    status: str
    fva: pd.DataFrame                      # index reaction id; columns minimum, maximum
    gene_deletion: pd.Series                # index gene id; objective after deletion (NaN = infeasible)
    reaction_deletion: pd.Series            # index reaction id; objective after deletion
    fva_fraction: float
    fva_subset: Optional[List[str]]         # None = all reactions
    timings_s: Dict[str, float] = field(default_factory=dict)
    tool: str = "cobrapy"
    tool_version: str = cobra.__version__

    def write(self, outdir: str) -> None:
        os.makedirs(outdir, exist_ok=True)
        meta = {"model_id": self.model_id, "tool": self.tool, "tool_version": self.tool_version, "solver": self.solver,
                "objective_value": self.objective_value, "status": self.status, "fva_fraction": self.fva_fraction,
                "fva_subset_size": None if self.fva_subset is None else len(self.fva_subset), "timings_s": self.timings_s}
        with open(os.path.join(outdir, "objective.json"), "w") as fh:
            json.dump(meta, fh, indent=2)
        self.fva.to_csv(os.path.join(outdir, "fva.tsv"), sep="\t")
        self.gene_deletion.rename("objective").to_csv(os.path.join(outdir, "gene_deletion.tsv"), sep="\t")
        self.reaction_deletion.rename("objective").to_csv(os.path.join(outdir, "reaction_deletion.tsv"), sep="\t")


def frog(model: cobra.Model, solver: str = "glpk", fva_fraction: float = 1.0, max_fva_reactions: Optional[int] = None,
         max_deletion_reactions: Optional[int] = None, processes: int = 1, seed: int = 0) -> FrogReport:
    model = model.copy()
    model.solver = solver
    timings = {}
    t = time.time()
    sol = model.optimize()
    timings["fba_s"] = time.time() - t
    obj = float(sol.objective_value) if sol.status == "optimal" else float("nan")

    rng = np.random.default_rng(seed)
    rxn_ids = [r.id for r in model.reactions]
    fva_subset = None
    if max_fva_reactions is not None and len(rxn_ids) > max_fva_reactions:
        fva_subset = sorted(rng.choice(rxn_ids, max_fva_reactions, replace=False).tolist())
    t = time.time()
    if sol.status == "optimal":
        fva = flux_variability_analysis(model, reaction_list=fva_subset, fraction_of_optimum=fva_fraction, processes=processes)
    else:
        fva = pd.DataFrame(columns=["minimum", "maximum"])
    timings["fva_s"] = time.time() - t

    t = time.time()
    gd = single_gene_deletion(model, processes=processes) if sol.status == "optimal" else None
    gene_del = pd.Series({next(iter(i)): g for i, g in zip(gd["ids"], gd["growth"])}) if gd is not None else pd.Series(dtype=float)
    timings["gene_deletion_s"] = time.time() - t

    del_subset = rxn_ids
    if max_deletion_reactions is not None and len(rxn_ids) > max_deletion_reactions:
        del_subset = sorted(rng.choice(rxn_ids, max_deletion_reactions, replace=False).tolist())
    t = time.time()
    rd = single_reaction_deletion(model, reaction_list=del_subset, processes=processes) if sol.status == "optimal" else None
    rxn_del = pd.Series({next(iter(i)): g for i, g in zip(rd["ids"], rd["growth"])}) if rd is not None else pd.Series(dtype=float)
    timings["reaction_deletion_s"] = time.time() - t

    return FrogReport(model_id=model.id, solver=solver, objective_value=obj, status=sol.status, fva=fva,
                      gene_deletion=gene_del, reaction_deletion=rxn_del, fva_fraction=fva_fraction,
                      fva_subset=fva_subset, timings_s=timings)


def compare(a: FrogReport, b: FrogReport, abs_tol: float = 1e-6, rel_tol: float = 1e-6, growth_tol: float = 1e-6) -> Dict[str, object]:
    """Compare two FROG reports. Reports max absolute/relative differences and the number of
    growth/no-growth disagreements (objective above vs below growth_tol) for deletions."""
    out: Dict[str, object] = {"model_id": a.model_id, "a": f"{a.tool}/{a.solver}", "b": f"{b.tool}/{b.solver}"}
    oa, ob = a.objective_value, b.objective_value
    out["objective_abs_diff"] = abs(oa - ob) if np.isfinite(oa) and np.isfinite(ob) else float("nan")
    out["objective_rel_diff"] = out["objective_abs_diff"] / max(abs(oa), abs(ob), 1e-12)
    out["objective_agree"] = bool(out["objective_abs_diff"] <= max(abs_tol, rel_tol * max(abs(oa), abs(ob))))

    def _series_cmp(sa: pd.Series, sb: pd.Series, name: str) -> None:
        common = sa.index.intersection(sb.index)
        va = sa.loc[common].astype(float).fillna(0.0).to_numpy()
        vb = sb.loc[common].astype(float).fillna(0.0).to_numpy()
        diff = np.abs(va - vb)
        out[f"{name}_n"] = int(len(common))
        out[f"{name}_max_abs_diff"] = float(diff.max()) if len(diff) else float("nan")
        out[f"{name}_n_exceeding_tol"] = int(np.sum(diff > np.maximum(abs_tol, rel_tol * np.maximum(np.abs(va), np.abs(vb)))))
        out[f"{name}_growth_call_disagreements"] = int(np.sum((va > growth_tol) != (vb > growth_tol)))

    _series_cmp(a.gene_deletion, b.gene_deletion, "gene_deletion")
    _series_cmp(a.reaction_deletion, b.reaction_deletion, "reaction_deletion")
    common = a.fva.index.intersection(b.fva.index)
    if len(common):
        dmin = np.abs(a.fva.loc[common, "minimum"].to_numpy() - b.fva.loc[common, "minimum"].to_numpy())
        dmax = np.abs(a.fva.loc[common, "maximum"].to_numpy() - b.fva.loc[common, "maximum"].to_numpy())
        out["fva_n"] = int(len(common))
        out["fva_max_abs_diff"] = float(max(dmin.max(), dmax.max()))
        out["fva_n_exceeding_1e-6"] = int(np.sum((dmin > 1e-6) | (dmax > 1e-6)))
        out["fva_n_exceeding_1e-3"] = int(np.sum((dmin > 1e-3) | (dmax > 1e-3)))
    else:
        out["fva_n"] = 0
    return out
