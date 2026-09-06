"""For every carbon-source condition on which a draft's wild type fails to grow: the minimal reaction set from the
CarveMe universe that restores growth, next to the genes the fitness data say the organism actually uses there.

The model is built as in the best benchmark arm (gap-filled draft + universe patches + model patches + accepted
gene-rule patches + medium completion). For each failing condition the carbon source's exchange is added if the
draft lacks it (exchanges are otherwise never gap-filled), the medium is applied, the carbon source opened at
-10 mmol/gDW/h, and gembench.gapfill finds the cheapest additions (transport 1.5x, no oxygen release). The
condition-specific genes are those with fitness below the threshold on this condition and a median fitness
above -1 on the other conditions where the wild-type model grows: the catabolic route the organism used.

The two lists are the evidence for a model patch with genes (data/reference/model_patches_v*.json): a reaction
the gap-filler proposes is adopted only when the genome annotation and the condition-specific genes support it.

Usage: python scripts/condition_gapfill.py [--orgs Btheta,Putida,MR1,Smeli] [--time-limit 240] [--min-growth 0.05]
Output: results/carbon_fitness_multi/condition_gapfill_<org>.tsv and condition_gapfill_<org>.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import os
import sys
import time

import cobra
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
logging.getLogger("cobra").setLevel(logging.ERROR)

from gembench.fitness_browser import base_medium, carbon_source_conditions, load_organism, replicate_averaged_fitness  # noqa: E402
from gembench.gapfill import _is_transport, gapfill  # noqa: E402
from gembench.media import apply_medium  # noqa: E402
from gembench.patches import apply_gpr_patches, apply_model_patches, apply_universe_patches, load_patch_files  # noqa: E402
from gembench.protocols.carbon_fitness_generic import GenericParams, complete_medium_transport  # noqa: E402
from run_carbon_fitness_generic import gene_map_for, load_model  # noqa: E402

OUT = os.path.join(ROOT, "results", "carbon_fitness_multi")
UNIVERSE = os.path.join(ROOT, "external", "carveme", "universe_bacteria.xml.gz")


def ensure_exchange(model: cobra.Model, universe: cobra.Model, met_id: str) -> bool:
    """Add EX_<met>_e (and the extracellular metabolite) from the universe if the model lacks it."""
    ex_id = f"EX_{met_id}_e"
    if ex_id in model.reactions:
        return False
    e_id = f"{met_id}_e"
    if e_id not in model.metabolites:
        if e_id not in universe.metabolites:
            return False
        u = universe.metabolites.get_by_id(e_id)
        model.add_metabolites([cobra.Metabolite(e_id, name=u.name, formula=u.formula, charge=u.charge, compartment="e")])
    ex = cobra.Reaction(ex_id, name=f"{met_id} exchange (added for the condition)", lower_bound=0.0, upper_bound=1000.0)
    ex.add_metabolites({model.metabolites.get_by_id(e_id): -1.0})
    model.add_reactions([ex])
    return True


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--orgs", default="Btheta,Putida,MR1,Smeli")
    ap.add_argument("--variant", default="gapfilled")
    ap.add_argument("--patch", default=os.path.join(ROOT, "data/reference/universe_patches_v0.1.json"))
    ap.add_argument("--model-patch", default=os.path.join(ROOT, "data/reference/model_patches_v0.2.json"))
    ap.add_argument("--gpr-patch", default=",".join(os.path.join(ROOT, f) for f in ["data/reference/gpr_patches_v0.2.json", "data/reference/gpr_patches_v0.3.json"]))
    ap.add_argument("--time-limit", type=float, default=240.0)
    ap.add_argument("--min-growth", type=float, default=0.05)
    ap.add_argument("--conditions", nargs="*", default=None, help="restrict to these condition names")
    args = ap.parse_args()
    patches = json.load(open(args.patch))
    model_patches = json.load(open(args.model_patch)) if args.model_patch else None
    gpr = load_patch_files(args.gpr_patch.split(","))
    with gzip.open(UNIVERSE, "rt") as fh:
        universe = cobra.io.read_sbml_model(fh)
    # universe-level reaction patches without a taxon scope apply to the universe itself (e.g. directionality), so
    # that the gap-filler cannot re-introduce a corrected error
    u_applied = apply_universe_patches(universe, "_universe_", {"patches": [pt for pt in patches["patches"] if not pt.get("scope_orgs")]})
    print(f"universe: {len(universe.reactions)} reactions; patches applied to the universe: {[a['reaction'] for a in u_applied]}", flush=True)
    params = GenericParams()

    for org in args.orgs.split(","):
        fb = load_organism(org)
        model, _, _ = load_model(org, args.variant)
        model.solver = "glpk"
        apply_universe_patches(model, org, patches)
        gm = gene_map_for(org, model, set(fb.genes["sysName"]), args.variant)
        if model_patches and apply_model_patches(model, org, model_patches, gm, verbose=False):
            gm = gene_map_for(org, model, set(fb.genes["sysName"]), args.variant)
        if any(a.get("genes_added") for a in apply_gpr_patches(model, org, gm, gpr, verbose=False)):
            gm = gene_map_for(org, model, set(fb.genes["sysName"]), args.variant)
        conds = [c for c in carbon_source_conditions(fb) if c.bigg_ids]
        if args.conditions:
            conds = [c for c in conds if c.name in set(args.conditions)]
        complete_medium_transport(model, sorted({c.media for c in conds}), params.medium_completion_exclude)
        fit = replicate_averaged_fitness(fb, conds)
        desc = dict(zip(fb.genes["sysName"], fb.genes["desc"])) if "desc" in fb.genes.columns else {}

        # wild-type growth per condition on the current model
        wt = {}
        for c in conds:
            with model:
                apply_medium(model, base_medium(c.media), close_all=True)
                for ex_id in c.exchanges:
                    if ex_id in model.reactions:
                        model.reactions.get_by_id(ex_id).lower_bound = params.carbon_uptake
                v = model.slim_optimize()
                wt[c.key] = 0.0 if v is None or np.isnan(v) else float(v)
        growing = [c for c in conds if wt[c.key] >= params.growth_threshold]
        failing = [c for c in conds if wt[c.key] < params.growth_threshold]
        print(f"== {org}: {len(failing)} failing of {len(conds)} conditions; {len(growing)} growing", flush=True)
        rows, details = [], []
        for c in failing:
            t0 = time.time()
            with model:
                added_ex = [m for m in c.bigg_ids if ensure_exchange(model, universe, m)]
                med = base_medium(c.media)
                res = gapfill(model, universe, med, extra_uptakes={ex: params.carbon_uptake for ex in c.exchanges},
                              min_growth=args.min_growth, time_limit_s=args.time_limit, verbose=False)
            added = []
            for rid in res.added_reactions:
                r = universe.reactions.get_by_id(rid)
                added.append({"id": rid, "name": r.name, "reaction": r.reaction, "transport": _is_transport(r)})
            # condition-specific genes
            f_here = fit[c.key]
            others = [g.key for g in growing]
            f_other = fit[others].median(axis=1) if others else pd.Series(np.nan, index=fit.index)
            spec = fit.index[(f_here < params.fitness_threshold) & (f_other > -1.0)]
            spec = sorted(spec, key=lambda g: f_here[g])
            genes = [{"gene": g, "fitness": round(float(f_here[g]), 2), "median_other": round(float(f_other[g]), 2) if np.isfinite(f_other[g]) else None,
                      "in_model": g in set(gm.model_to_browser.values()), "desc": desc.get(g, "")} for g in spec[:25]]
            rows.append({"org": org, "condition": c.name, "media": c.media, "bigg_ids": ";".join(c.bigg_ids),
                         "exchange_added": ";".join(added_ex), "gapfill_status": res.status, "growth_after": round(res.growth_after, 4),
                         "rejected_energy_cycle_sets": " | ".join(";".join(x) for x in res.rejected_energy_cycles),
                         "n_added": len(added), "added_reactions": ";".join(a["id"] for a in added),
                         "added_transport": ";".join(a["id"] for a in added if a["transport"]),
                         "n_condition_specific_genes": len(spec),
                         "condition_specific_genes": ";".join(f"{g['gene']}({g['fitness']})" for g in genes[:12]),
                         "seconds": round(time.time() - t0)})
            details.append({"condition": c.name, "media": c.media, "bigg_ids": c.bigg_ids, "exchange_added": added_ex,
                            "gapfill": {"status": res.status, "growth_after": res.growth_after, "added": added, "seconds": res.seconds,
                                        "rejected_energy_cycle_sets": res.rejected_energy_cycles},
                            "condition_specific_genes": genes})
            print(f"   {c.name[:34]:34s} {res.status:12s} added={[a['id'] for a in added]} growth={res.growth_after:.3f} "
                  f"specific genes={len(spec)}: {', '.join(g['gene'] for g in genes[:6])} ({time.time()-t0:.0f}s)", flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(OUT, f"condition_gapfill_{org}.tsv"), sep="\t", index=False)
        json.dump(details, open(os.path.join(OUT, f"condition_gapfill_{org}.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
