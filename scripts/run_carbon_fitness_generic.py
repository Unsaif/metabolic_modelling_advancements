"""Run the multi-organism carbon-source fitness benchmark.

Usage: python scripts/run_carbon_fitness_generic.py [--orgs Btheta,Putida,MR1,Smeli,Keio] [--max-conditions N]
                                                    [--no-drop-rich] [--processes 2]

Models: EMBL GEMs (CarveMe) for Btheta/Putida/MR1/Smeli; iML1515 for Keio (control: same pipeline,
fresh Fitness Browser download and this repository's own condition/media mapping, instead of the
Bernstein et al. files used in Sprint 1).
Results: results/carbon_fitness_multi/<org>/<model>/ (card JSON+MD, npz matrices, condition table, gene map).
"""
from __future__ import annotations

import argparse
import glob
import gzip
import json
import logging
import os
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
logging.getLogger("cobra").setLevel(logging.ERROR)

import cobra  # noqa: E402

from gembench import metrics as M  # noqa: E402
from gembench.cards import BenchmarkCard, LeakageCard, ModelProvenance, now, sha256_of  # noqa: E402
from gembench.fitness_browser import carbon_source_conditions, load_organism  # noqa: E402
from gembench.gene_mapping import GeneMap, build_gene_map  # noqa: E402
from gembench.protocols import carbon_fitness_generic as P  # noqa: E402
from gembench.protocols.ecoli_carbon_fitness import BW25113_DELETED_GENES  # noqa: E402

EMBL = {"Btheta": "Bacteroides_thetaiotaomicron_VPI_5482", "Putida": "Pseudomonas_putida_KT2440",
        "MR1": "Shewanella_oneidensis_MR_1", "Smeli": "Sinorhizobium_meliloti_1021"}
OUT = os.path.join(ROOT, "results", "carbon_fitness_multi")


def load_model(org: str, variant: str = "shipped"):
    if org != "Keio" and variant == "gapfilled":
        path = os.path.join(ROOT, "models", "gapfilled", f"{org}.xml.gz")
        with gzip.open(path, "rt") as fh:
            m = cobra.io.read_sbml_model(fh)
        info = json.load(open(os.path.join(ROOT, "models", "gapfilled", f"{org}_gapfill.json")))
        return m, path, ("EMBL GEMs (CarveMe draft, Machado et al. 2018) minimally gap-filled by gembench.gapfill for growth on "
                         f"{info['medium']} + {info['reference_carbon_exchange']}: added {[a['id'] for a in info['added_reactions']]}")
    if org == "Keio":
        path = os.path.join(ROOT, "models", "iML1515", "iML1515.xml")
        if not os.path.exists(path):
            path = glob.glob(os.path.join(ROOT, "models", "iML1515*"))[0]
            if os.path.isdir(path):
                path = glob.glob(os.path.join(path, "*.xml"))[0]
        m = cobra.io.read_sbml_model(path)
        return m, path, "BiGG iML1515 (Monk et al. 2017) via github.com/dbernste/E_coli_GEM_validation Models/"
    path = glob.glob(os.path.join(ROOT, "models", "embl", f"{EMBL[org]}*.xml.gz"))[0]
    with gzip.open(path, "rt") as fh:
        m = cobra.io.read_sbml_model(fh)
    m.id = m.id or EMBL[org]
    return m, path, "EMBL GEMs (CarveMe draft, Machado et al. 2018; github.com/cdanielmachado/embl_gems)"


def gene_map_for(org: str, model, browser_sysnames) -> GeneMap:
    if org == "Keio":   # iML1515 gene ids are b-numbers = Fitness Browser sysName
        ids = [g.id for g in model.genes]
        mp = {g: g for g in ids if g in browser_sysnames}
        return GeneMap(org_id=org, model_to_browser=mp, unmapped_model_genes=[g for g in ids if g not in mp],
                       stats={"model_genes": len(ids), "mapped": len(mp)}, provenance={"method": "identity (b-numbers)"})
    return build_gene_map(org, [g.id for g in model.genes], browser_sysnames)


def score(res: P.GenericResult, use_conditions: np.ndarray) -> dict:
    sim = res.sim_growth[:, use_conditions]
    fit = res.fitness[:, use_conditions]
    ok_rows = np.isfinite(fit).any(axis=1)
    sim, fit = sim[ok_rows], fit[ok_rows]
    out = {"n_genes": int(sim.shape[0]), "n_conditions": int(sim.shape[1]),
           "n_gene_condition_pairs": int(np.isfinite(fit).sum())}
    if sim.size == 0:
        return out
    st, ft = res.params.growth_threshold, res.params.fitness_threshold
    for name, fn in [("aucpr_bernstein", lambda s, f: M.aucpr_bernstein(s, f, st)),
                     ("aucpr_standard", lambda s, f: M.aucpr_standard(s, f, ft)),
                     ("auroc_standard", lambda s, f: M.auroc_standard(s, f, ft)),
                     ("mcc", lambda s, f: M.mcc(s, f, st, ft)),
                     ("balanced_accuracy", lambda s, f: M.balanced_accuracy(s, f, st, ft)),
                     ("accuracy", lambda s, f: M.accuracy(s, f, st, ft))]:
        try:
            pt, lo, hi = M.bootstrap_ci(fn, sim, fit, n_boot=500)
        except Exception as e:  # noqa: BLE001
            pt, lo, hi = float("nan"), float("nan"), float("nan")
        out[name] = {"point": pt, "ci95": [lo, hi]}
    out["confusion"] = M.confusion(sim, fit, st, ft)
    return out


def per_condition(res: P.GenericResult) -> pd.DataFrame:
    rows = []
    st, ft = res.params.growth_threshold, res.params.fitness_threshold
    for j, c in enumerate(res.conditions):
        s, f = res.sim_growth[:, j], res.fitness[:, j]
        ok = np.isfinite(f)
        d = {"condition": c.name, "media": c.media, "wt_growth": float(res.wt_growth[j]),
             "n_genes": int(ok.sum())}
        if res.wt_growth[j] >= st and ok.sum() > 0:
            d.update(M.confusion(s[ok], f[ok], st, ft))
            d["mcc"] = M.mcc(s[ok], f[ok], st, ft)
            d["aucpr_bernstein"] = M.aucpr_bernstein(s[ok], f[ok], st)
            d["auroc_standard"] = M.auroc_standard(s[ok], f[ok], ft)
        rows.append(d)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--orgs", default="Btheta,Putida,MR1,Smeli,Keio")
    ap.add_argument("--max-conditions", type=int, default=None)
    ap.add_argument("--no-drop-rich", action="store_true")
    ap.add_argument("--processes", type=int, default=2)
    ap.add_argument("--solver", default="glpk")
    ap.add_argument("--variant", default="shipped", choices=["shipped", "gapfilled"])
    args = ap.parse_args()

    for org in args.orgs.split(","):
        t0 = time.time()
        fb = load_organism(org)
        model, mpath, msource = load_model(org, args.variant)
        gm = gene_map_for(org, model, set(fb.genes["sysName"]))
        gm.stats["mapped_with_fitness_data"] = sum(1 for v in gm.model_to_browser.values() if v in set(fb.fitness.index))
        conds = carbon_source_conditions(fb)
        print(f"== {org}: model {model.id} ({len(model.genes)} genes); gene map {gm.stats}; "
              f"{len(conds)} carbon-source conditions, {sum(1 for c in conds if c.bigg_ids)} mapped", flush=True)
        params = P.GenericParams(drop_rich_medium_essentials=not args.no_drop_rich, processes=args.processes,
                                 solver=args.solver, max_conditions=args.max_conditions,
                                 knockout_genes=BW25113_DELETED_GENES if org == "Keio" else [])
        res = P.run(model, fb, conds, gm, params)

        outdir = os.path.join(OUT, org, f"{model.id}__{args.variant}" if org != "Keio" else model.id)
        os.makedirs(outdir, exist_ok=True)
        grows = res.wt_growth >= params.growth_threshold
        results = {
            "condition_level": {"n_conditions_mapped": int(len(res.conditions)),
                                "n_conditions_wt_grows": int(grows.sum()),
                                "wt_growth_recall": float(grows.mean()) if len(grows) else float("nan"),
                                "conditions_with_absent_exchange": int(sum(1 for c in res.conditions if res.missing_carbon_exchanges.get(c.key)))},
            "gene_level_conditions_where_wt_grows": score(res, grows),
            "gene_level_all_mapped_conditions": score(res, np.ones(len(res.conditions), dtype=bool)),
            "gene_map": gm.stats, "counts": res.counts, "timings_s": res.timings_s,
            "dropped_rich_medium_essentials": len(res.dropped_rich_essential_genes),
        }
        unmapped = [c for c in conds if not c.bigg_ids]
        card = BenchmarkCard(
            benchmark=f"carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism {org}",
            created=now(),
            dataset_provenance={**fb.provenance, "conditions_unmapped": "; ".join(sorted(set(c.name for c in unmapped)))},
            model=ModelProvenance(model_id=model.id, file=os.path.relpath(mpath, ROOT), source=msource,
                                  n_reactions=len(model.reactions), n_metabolites=len(model.metabolites),
                                  n_genes=len(model.genes), sha256=sha256_of(mpath)),
            protocol={"variant": args.variant, "params": res.params.__dict__, "media_mapping": "data/reference/fitness_browser_media_bigg.tsv",
                      "carbon_source_mapping": "data/reference/fitness_browser_carbon_sources_bigg.tsv",
                      "gene_mapping": gm.provenance,
                      "condition_selection": "expGroup == 'carbon source'; condition_2 empty or DMSO; media with a BiGG mapping; replicates averaged per condition x medium",
                      "scoring": "gene-level metrics pooled over conditions where the wild-type model grows (>= growth_threshold); "
                                 "condition-level recall = fraction of experimentally growing carbon sources on which the wild-type model grows"},
            leakage=LeakageCard(
                ground_truth_used_in_model_curation="no for EMBL draft models (automated reconstruction from genome annotation); "
                                                    "for iML1515: partly (E. coli curation used phenotype data)",
                ground_truth_public_since="Fitness Browser releases 2015-2018 (Price et al. 2018)",
                frontier_model_training_exposure="Fitness Browser tables are public and partly in training corpora; the mapping tables here are new",
                held_out_recommendation="unpublished RB-TnSeq experiments, or organisms added to the Browser after the model's training cut-off",
                notes=["Draft models are untouched by any phenotype data, so this is a true prospective test of automated reconstruction."]),
            results=results,
            warnings=[f"{len(unmapped)} of {len(conds)} conditions have no BiGG mapping"] +
                     [f"medium '{m}' components absent from the model: {v}" for m, v in res.missing_medium_components.items() if v],
        )
        card.write(os.path.join(outdir, "card.json"), os.path.join(outdir, "card.md"))
        np.savez_compressed(os.path.join(outdir, "matrices.npz"), sim_growth=res.sim_growth, wt_growth=res.wt_growth,
                            fitness=res.fitness, model_genes=np.array(res.model_genes), browser_genes=np.array(res.browser_genes),
                            conditions=np.array([c.key for c in res.conditions]))
        res.condition_table().to_csv(os.path.join(outdir, "conditions.tsv"), sep="\t", index=False)
        per_condition(res).to_csv(os.path.join(outdir, "per_condition_metrics.tsv"), sep="\t", index=False)
        pd.DataFrame({"model_gene": list(gm.model_to_browser), "browser_gene": list(gm.model_to_browser.values())}).to_csv(
            os.path.join(outdir, "gene_map.tsv"), sep="\t", index=False)
        with open(os.path.join(outdir, "gene_map_unmapped.txt"), "w") as fh:
            fh.write("\n".join(gm.unmapped_model_genes))
        with open(os.path.join(outdir, "dropped_rich_medium_essentials.txt"), "w") as fh:
            fh.write("\n".join(res.dropped_rich_essential_genes))
        g = results["gene_level_conditions_where_wt_grows"]
        print(f"== {org} done in {time.time()-t0:.0f}s: WT grows on {grows.sum()}/{len(grows)} mapped conditions; "
              f"gene-level AUC-PR(Bernstein)={g.get('aucpr_bernstein', {}).get('point', float('nan')):.3f} "
              f"MCC={g.get('mcc', {}).get('point', float('nan')):.3f} over {g.get('n_gene_condition_pairs')} pairs", flush=True)


if __name__ == "__main__":
    main()
