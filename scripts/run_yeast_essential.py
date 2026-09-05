"""Run the yeast essential-gene protocol on yeast-GEM and write a benchmark card."""
from __future__ import annotations

import argparse
import json
import os
import sys

import cobra
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gembench import metrics as M  # noqa: E402
from gembench.cards import BenchmarkCard, LeakageCard, ModelProvenance, now, sha256_of  # noqa: E402
from gembench.datasets import load_yeast_essentiality  # noqa: E402
from gembench.protocols import yeast_essential as P  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=os.path.join(ROOT, "models", "yeast-GEM.xml"))
    ap.add_argument("--name", default="yeast-GEM")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--processes", type=int, default=2)
    args = ap.parse_args()

    data = load_yeast_essentiality()
    model = cobra.io.read_sbml_model(args.model)
    model.id = args.name
    print(f"== {args.name}: {len(model.reactions)} rxns, {len(model.metabolites)} mets, {len(model.genes)} genes")
    res = P.run(model, data, P.YeastParams(processes=args.processes))

    # metrics: 'sim' = growth ratio, 'fit' = -1 for inviable, +1 for viable, so that the shared
    # helpers (which call fitness < -2 'important') can be reused with fit_thresh = 0.
    sim = res.growth_ratio[:, None]
    fit = np.where(res.exp_inviable, -1.0, 1.0)[:, None]
    tol = res.params.ko_tol
    r = {"mcc": M.bootstrap_ci(lambda s, f: M.mcc(s, f, tol, 0.0), sim, fit, args.n_boot),
         "accuracy": M.bootstrap_ci(lambda s, f: M.accuracy(s, f, tol, 0.0), sim, fit, args.n_boot),
         "balanced_accuracy": M.bootstrap_ci(lambda s, f: M.balanced_accuracy(s, f, tol, 0.0), sim, fit, args.n_boot),
         "aucpr_inviable_as_positive": M.bootstrap_ci(lambda s, f: M.aucpr_standard(s, f, 0.0), sim, fit, args.n_boot)}
    results = {k: {"point": v[0], "ci95": [v[1], v[2]]} for k, v in r.items()}
    results["counts"] = res.counts
    results["wt_growth"] = res.wt_growth
    results["missing_exchanges"] = res.missing_exchanges
    results["reference_yeastGEM_9.1.1_testResults"] = {"TP": 933, "TN": 65, "FP": 94, "FN": 15,
                                                       "source": "external/yeast-GEM/data/testResults/essentialGenes.tsv"}
    card = BenchmarkCard(
        benchmark="yeast_deletion_viability_v0 (yeast-GEM essentialGenes protocol)",
        created=now(),
        dataset_provenance=data.provenance,
        model=ModelProvenance(model_id=args.name, file=os.path.relpath(args.model, ROOT),
                              source="github.com/SysBioChalmers/yeast-GEM (main branch, model/yeast-GEM.xml)",
                              n_reactions=len(model.reactions), n_metabolites=len(model.metabolites),
                              n_genes=len(model.genes), sha256=sha256_of(os.path.realpath(args.model))),
        protocol=res.params.__dict__,
        leakage=LeakageCard(
            ground_truth_used_in_model_curation="yes — yeast-GEM has used this exact essential-gene test as a curation "
                                                "regression check since yeast7/8; treat as in-sample",
            ground_truth_public_since="2002 (Giaever et al.) / Stanford deletion project downloads",
            frontier_model_training_exposure="almost certainly in training corpora",
            held_out_recommendation="use condition-specific deletion phenotypes not used by the yeast-GEM test suite "
                                    "(e.g. Nichols-style chemical genomics for yeast, or newer Tn-seq/CRISPRi screens)"),
        results=results,
    )
    out = os.path.join(ROOT, "results", "yeast_essential")
    card.write(os.path.join(out, f"{args.name}.card.json"), os.path.join(out, f"{args.name}.card.md"))
    with open(os.path.join(out, f"{args.name}.classification.tsv"), "w") as fh:
        fh.write("gene\tclassification\tgrowth_ratio\n")
        for g, gr in zip(res.genes, res.growth_ratio):
            fh.write(f"{g}\t{res.classification[g]}\t{gr:.6g}\n")
    print(f"   accuracy {r['accuracy'][0]:.4f} [{r['accuracy'][1]:.4f}, {r['accuracy'][2]:.4f}] | MCC {r['mcc'][0]:.3f} "
          f"[{r['mcc'][1]:.3f}, {r['mcc'][2]:.3f}] | counts {res.counts}")


if __name__ == "__main__":
    main()
