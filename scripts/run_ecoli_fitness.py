"""Run the E. coli carbon-source fitness protocol on one or more models and write benchmark cards.

Usage: python scripts/run_ecoli_fitness.py iJR904 iAF1260 iJO1366 iML1515 [--variant mannitol]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import cobra
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gembench import metrics as M  # noqa: E402
from gembench.cards import BenchmarkCard, LeakageCard, ModelProvenance, now, sha256_of  # noqa: E402
from gembench.datasets import load_ecoli_bw25113_carbon_fitness  # noqa: E402
from gembench.protocols import ecoli_carbon_fitness as P  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODELS = os.path.join(ROOT, "models")
RESULTS = os.path.join(ROOT, "results", "ecoli_carbon_fitness")

MODEL_SOURCES = {
    "iJR904": "BiGG (2003, Reed et al.); copy from github.com/dbernste/E_coli_GEM_validation/Models",
    "iAF1260": "BiGG (2007, Feist et al.); copy from github.com/dbernste/E_coli_GEM_validation/Models",
    "iJO1366": "BiGG (2011, Orth et al.); copy from github.com/dbernste/E_coli_GEM_validation/Models",
    "iML1515": "BiGG (2017, Monk et al.); copy from github.com/dbernste/E_coli_GEM_validation/Models",
    "iML1515_corrected": "Bernstein et al. 2023 iML1515 with all corrections (Analysis/iML1515_model_adjusted_all_corrections.xml)",
}
MODEL_FILES = {
    "iML1515_corrected": os.path.join(ROOT, "external", "E_coli_GEM_validation", "Analysis",
                                      "iML1515_model_adjusted_all_corrections.xml"),
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("models", nargs="+")
    ap.add_argument("--variant", default="bernstein_as_coded",
                    choices=["bernstein_as_coded", "mannitol_excluded", "no_exclusion"])
    ap.add_argument("--processes", type=int, default=2)
    ap.add_argument("--n-boot", type=int, default=1000)
    args = ap.parse_args()

    data = load_ecoli_bw25113_carbon_fitness()
    print(f"dataset: {data.provenance['n_genes']} genes, {len(data.carbon_sources)} carbon sources: {data.carbon_sources}")

    remove = {"bernstein_as_coded": ["man", "sucr"], "mannitol_excluded": ["mnl", "sucr"], "no_exclusion": []}[args.variant]
    params = P.ProtocolParams(carbon_sources_to_remove=remove, processes=args.processes)

    summary = {}
    for name in args.models:
        path = MODEL_FILES.get(name, os.path.join(MODELS, f"{name}.xml"))
        t = time.time()
        model = cobra.io.read_sbml_model(path)
        model.id = name
        print(f"\n== {name}: {len(model.reactions)} rxns, {len(model.metabolites)} mets, {len(model.genes)} genes (loaded {time.time()-t:.1f}s)")
        res = P.run(model, data, params)
        sim, fit = res.sim_growth, res.fitness
        gt, ft = params.growth_threshold, params.fitness_threshold

        r = {}
        r["aucpr_bernstein"] = M.bootstrap_ci(lambda s, f: M.aucpr_bernstein(s, f, gt), sim, fit, args.n_boot)
        r["aucpr_standard"] = M.bootstrap_ci(lambda s, f: M.aucpr_standard(s, f, ft), sim, fit, args.n_boot)
        r["auroc_standard"] = M.bootstrap_ci(lambda s, f: M.auroc_standard(s, f, ft), sim, fit, args.n_boot)
        r["mcc"] = M.bootstrap_ci(lambda s, f: M.mcc(s, f, gt, ft), sim, fit, args.n_boot)
        r["balanced_accuracy"] = M.bootstrap_ci(lambda s, f: M.balanced_accuracy(s, f, gt, ft), sim, fit, args.n_boot)
        r["accuracy"] = M.bootstrap_ci(lambda s, f: M.accuracy(s, f, gt, ft), sim, fit, args.n_boot)
        conf = M.confusion(sim, fit, gt, ft)
        per_carbon = {c: {"aucpr_bernstein": M.aucpr_bernstein(sim[:, j], fit[:, j], gt),
                          "mcc": M.mcc(sim[:, j], fit[:, j], gt, ft),
                          "wt_growth": float(res.wt_growth[j]),
                          "n_pred_no_growth": int((sim[:, j] < gt).sum()),
                          "n_exp_important": int((fit[:, j] < ft).sum())}
                     for j, c in enumerate(res.carbon_sources)}

        results = {k: {"point": v[0], "ci95": [v[1], v[2]]} for k, v in r.items()}
        results["confusion_growth_vs_important"] = conf
        results["counts"] = res.counts
        results["dropped_strain_genes"] = res.dropped_strain_genes
        results["n_dropped_rich_medium_essentials"] = len(res.dropped_rich_essential_genes)
        results["dropped_carbon_sources"] = res.dropped_carbon_sources
        results["missing_medium_components"] = res.missing_medium_components
        results["missing_carbon_exchanges"] = res.missing_carbon_exchanges
        results["per_carbon_source"] = per_carbon
        results["timings_s"] = res.timings_s

        card = BenchmarkCard(
            benchmark="ecoli_carbon_fitness_v0 (Bernstein et al. 2023 protocol)",
            created=now(),
            dataset_provenance=data.provenance,
            model=ModelProvenance(model_id=name, file=os.path.relpath(path, ROOT), source=MODEL_SOURCES.get(name, "unknown"),
                                  n_reactions=len(model.reactions), n_metabolites=len(model.metabolites),
                                  n_genes=len(model.genes), sha256=sha256_of(path)),
            protocol={"variant": args.variant, **{k: (v if not hasattr(v, "uptakes") else {"medium": v.name, "n_components": len(v.uptakes)})
                                                    for k, v in res.params.__dict__.items()}},
            leakage=LeakageCard(
                ground_truth_used_in_model_curation=("yes — this model was corrected against this dataset (in-sample)" if name == "iML1515_corrected"
                                                     else "no — BiGG E. coli models predate the 2018 Fitness Browser release; not used in their curation"),
                ground_truth_public_since="2018 (Price et al., Nature; Fitness Browser)",
                frontier_model_training_exposure="almost certainly present in pre-2026 training corpora (public TSVs, GitHub copies)",
                held_out_recommendation="for any claim about model-driven curation, hold out organisms/conditions absent from the "
                                        "Fitness Browser at training cut-off, or use newly generated phenotypes",
                notes=["Sucrose/'man' exclusion: notebook comment says mannitol, code removes D-mannose ('man'); "
                       "variant flag records which was used."]),
            results=results,
            warnings=[w for w in [
                f"missing medium components: {res.missing_medium_components}" if res.missing_medium_components else "",
                f"missing carbon exchanges: {res.missing_carbon_exchanges}" if res.missing_carbon_exchanges else ""] if w],
        )
        out = os.path.join(RESULTS, args.variant)
        card.write(os.path.join(out, f"{name}.card.json"), os.path.join(out, f"{name}.card.md"))
        np.savez_compressed(os.path.join(out, f"{name}.arrays.npz"), sim_growth=sim, fitness=fit, wt_growth=res.wt_growth,
                            genes=np.array(res.genes), carbon_sources=np.array(res.carbon_sources))
        summary[name] = {k: results[k] for k in ["aucpr_bernstein", "aucpr_standard", "mcc", "balanced_accuracy", "accuracy"]}
        summary[name]["counts"] = res.counts
        print(f"   AUC-PR (Bernstein) {r['aucpr_bernstein'][0]:.3f} [{r['aucpr_bernstein'][1]:.3f}, {r['aucpr_bernstein'][2]:.3f}] | "
              f"AUC-PR (std) {r['aucpr_standard'][0]:.3f} | MCC {r['mcc'][0]:.3f} | bal.acc {r['balanced_accuracy'][0]:.3f} | "
              f"acc {r['accuracy'][0]:.3f} | genes {res.counts['genes_after_adjustment']} | carbon {len(res.carbon_sources)} | {res.timings_s['total_s']:.0f}s")
    with open(os.path.join(RESULTS, args.variant, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
