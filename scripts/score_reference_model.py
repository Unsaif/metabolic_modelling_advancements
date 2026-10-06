"""Score a hand-curated reference model with the transfer-study protocol (post hoc, descriptive).

Paper 1 asks how close corrected drafts get to hand-curated models. This script scores a published curated model
exactly as the transfer study scored its arms: the organism's media and carbon-source tables, the fixed parameters
(carbon uptake -10, growth threshold 1e-3, fitness threshold -2, medium completion, GLPK) and the same scoring code.
It does not gap-fill the curated model and changes nothing in it. Gene identifiers of BiGG models are locus tags,
which match the Fitness Browser sysNames, so genes map by identity.

This is context for the paper, not a test of the corrections: the organisms scored here (P. putida, E. coli) were
used in development, and their fitness data were seen long before.

Usage:
  python scripts/score_reference_model.py --org Putida --model models/bigg/iJN1463.xml --label iJN1463
  python scripts/score_reference_model.py --org Keio --model models/iML1515/iML1515.xml --label iML1515
  python -I scripts/score_reference_model.py --org MR1 --model models/curated/iSO783/iSO783_bigg_view.xml.gz --label iSO783
  python -I scripts/score_reference_model.py --org Smeli --model models/curated/iGD1575/iGD1575_bigg_view.xml.gz \
      --label iGD1575 --gene-normalize sm_prefix
Models outside BiGG are first given BiGG identifiers by scripts/translate_curated_model.py (relabelling only).
Writes results/transfer_v1/reference_models/<org>/<label>/ (card.json, matrices.npz, conditions.tsv,
per_condition_metrics.tsv), in the same format as the transfer arms.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import run_transfer_study as RT  # noqa: E402
from gembench import fitness_browser as FB  # noqa: E402
from gembench.cards import BenchmarkCard, LeakageCard, ModelProvenance, now, sha256_of  # noqa: E402
from gembench.gene_mapping import GeneMap  # noqa: E402
from gembench.protocols import carbon_fitness_generic as P  # noqa: E402
from run_carbon_fitness_generic import per_condition, score  # noqa: E402

OUT = os.path.join(ROOT, "results", "transfer_v1", "reference_models")
# Organisms outside the transfer study's configuration (never written to the locked study files).
EXTRA_CONFIG = {"Keio": {"role": "reference", "fb_dir": "data/fitness_browser/Keio"}}


def config(org):
    if org in EXTRA_CONFIG:
        c = dict(EXTRA_CONFIG[org])
        c["org"] = org
        return c
    return RT.load_config(org)


GENE_NORMALIZE = {
    "identity": (lambda g: g, "identity (BiGG gene ids are locus tags = Fitness Browser sysName)"),
    "remove_underscore": (lambda g: g.replace("_", ""), "locus tag with the underscore removed (BT_0554 -> BT0554 = Fitness Browser sysName)"),
    "sm_prefix": (lambda g: ("SM_b" + g[3:]) if g[:3].lower() == "smb" else ("SM" + g[2:]) if g[:2].lower() == "sm" else g,
                  "S. meliloti locus tag in the Fitness Browser's form (smc04029 -> SMc04029, sma2091 -> SMa2091, smb21184 -> SM_b21184)"),
}


def identity_gene_map(org, model, sysnames, normalize="identity"):
    fn, how = GENE_NORMALIZE[normalize]
    ids = [g.id for g in model.genes]
    mp = {g: fn(g) for g in ids if fn(g) in sysnames}
    return GeneMap(org_id=org, model_to_browser=mp, unmapped_model_genes=[g for g in ids if g not in mp],
                   stats={"model_genes": len(ids), "mapped": len(mp)},
                   provenance={"method": how})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--org", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--source", default="")
    ap.add_argument("--processes", type=int, default=2)
    ap.add_argument("--gene-normalize", choices=sorted(GENE_NORMALIZE), default="identity")
    args = ap.parse_args()
    cfg = config(args.org)
    outdir = os.path.join(OUT, args.org, args.label)
    if os.path.exists(os.path.join(outdir, "card.json")):
        raise SystemExit(f"{outdir} exists; results are never overwritten")
    os.makedirs(outdir, exist_ok=True)
    fb = FB.load_organism(args.org, data_dir=os.path.dirname(os.path.join(ROOT, cfg["fb_dir"])))
    sysnames = set(fb.genes["sysName"])
    model = RT.read_sbml(os.path.join(ROOT, args.model))
    model.solver = "glpk"
    t0 = time.time()
    with RT.reference_tables(cfg) as (media_table, carbon_table):
        gm = identity_gene_map(args.org, model, sysnames, args.gene_normalize)
        gm.stats["mapped_with_fitness_data"] = sum(1 for v in gm.model_to_browser.values() if v in set(fb.fitness.index))
        conds = FB.carbon_source_conditions(fb)
        params = P.GenericParams(processes=args.processes, **RT.FIXED)
        res = P.run(model, fb, conds, gm, params, verbose=False)
    grows = res.wt_growth >= params.growth_threshold
    results = {"condition_level": {"n_conditions_mapped": int(len(res.conditions)), "n_conditions_wt_grows": int(grows.sum()),
                                   "conditions_with_absent_exchange": int(sum(1 for c in res.conditions if res.missing_carbon_exchanges.get(c.key)))},
               "gene_level_conditions_where_wt_grows": score(res, grows),
               "gene_map": gm.stats, "counts": res.counts, "timings_s": res.timings_s}
    leakage = LeakageCard(
        ground_truth_used_in_model_curation="unknown for the published curated model; its authors may have used these or related phenotypes",
        ground_truth_public_since="Fitness Browser releases include Price et al. 2018",
        frontier_model_training_exposure="not applicable (no AI step)",
        held_out_recommendation="reference point only; not a held-out test",
        notes=["Scored post hoc with the transfer study's fixed protocol, for context in Paper 1."])
    card = BenchmarkCard(
        benchmark=f"carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism {args.org}",
        created=now(), dataset_provenance=fb.provenance,
        model=ModelProvenance(model_id=model.id, file=args.model, source=args.source or args.label,
                              n_reactions=len(model.reactions), n_metabolites=len(model.metabolites), n_genes=len(model.genes),
                              sha256=sha256_of(os.path.join(ROOT, args.model))),
        protocol={"study": "transfer_v1 reference models (post hoc)", "arm": f"REF_{args.label}", "applied": [],
                  "params": res.params.__dict__, "media_mapping": os.path.relpath(media_table, ROOT),
                  "carbon_source_mapping": os.path.relpath(carbon_table, ROOT), "gene_mapping": gm.provenance,
                  "role": cfg.get("role", "reference")},
        leakage=leakage, results=results, warnings=[])
    card.write(os.path.join(outdir, "card.json"), os.path.join(outdir, "card.md"))
    np.savez_compressed(os.path.join(outdir, "matrices.npz"), sim_growth=res.sim_growth, wt_growth=res.wt_growth,
                        fitness=res.fitness, model_genes=np.array(res.model_genes), browser_genes=np.array(res.browser_genes),
                        conditions=np.array([c.key for c in res.conditions]))
    res.condition_table().to_csv(os.path.join(outdir, "conditions.tsv"), sep="\t", index=False)
    per_condition(res).to_csv(os.path.join(outdir, "per_condition_metrics.tsv"), sep="\t", index=False)
    g = results["gene_level_conditions_where_wt_grows"]
    print(f"{args.org}/{args.label}: WT grows {int(grows.sum())}/{len(grows)}; MCC={g.get('mcc', {}).get('point', float('nan')):.4f} "
          f"pairs={g.get('n_gene_condition_pairs')} genes mapped {gm.stats['mapped']}/{gm.stats['model_genes']} ({time.time()-t0:.0f}s)",
          flush=True)


if __name__ == "__main__":
    main()
