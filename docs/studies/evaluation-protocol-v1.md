# Evaluation protocol v1: development now, independent evaluation next

6 September 2026. This document governs interpretation of new experiments; it does not retrospectively preregister previous decisions. See the [exposure registry](../../data/studies/exposure_registry_v1.json) and [scientific audit](../reviews/2026-09-06-scientific-audit.md).

## What the present experiments can establish

The four draft-model organisms (Btheta, Putida, MR1 and Smeli) and the Keio E. coli control are development data. Their fitness outcomes, growth coverage and errors have already influenced model corrections or protocol choices. Repeating simulations with frozen settings can establish reproducibility and isolate effects of particular changes. It cannot estimate performance on independent organisms, prove a biological mechanism, or remove earlier selection bias.

The remaining four downloaded organisms (Bvulgatus_CL09T03C04, Koxy, SynE and DvH) have unknown prior exposure. Downloading a dataset does not establish that it was inspected, but absence of documented inspection does not establish independence either. Their phenotype contents remain quarantined during this sprint. They are not designated a holdout. No independent evaluation set has yet been established.

The [eight-run completion sprint](../sprints/2026-09-06-development-sprint.md) fixes both existing patch configurations before the new simulations. The [quinone sensitivity plan](../../data/studies/quinone_biomass_sensitivity_v1.json) separately fixes four interventions before their simulations. These are prespecified development experiments with local records, not public preregistrations. The preserved failed sensitivity attempt and its restart are part of the record.

## Measurements and decision rules

1. Report every planned arm and every failure. A solver error is an unknown result, not a gene-essentiality observation. Preserve failed attempts before retrying; document any changes to code, inputs or settings.
2. Report wild-type growth coverage as both numerator and denominator for all mapped conditions. Also report unmapped conditions where available. A model that predicts no wild-type growth cannot supply an interpretable relative knockout-growth assessment for that condition.
3. Report gene predictions only for finite observations and successful simulations. For pairwise comparisons, align gene and condition identities, verify identical observed fitness values and thresholds, and use the intersection of conditions where both wild types grow. Display the number of paired observations and the conditions gained or lost. An increase in each arm's separately computed MCC can be caused by a change in which conditions are scored.
4. Use the existing absolute growth cutoff of 0.001 and experimental fitness cutoff of -2, with inclusive comparisons, for these replication experiments. These are protocol assumptions; they are not newly established biological thresholds. No threshold tuning on the reported outcomes is permitted within these studies.
5. Report paired MCC differences and gene-cluster bootstrap intervals alongside the raw confusion matrices, growth coverage and changed predictions. The intervals describe resampling of the observed genes within one organism. They do not account for model selection, shared pathways, dependence between organisms or generalization to new datasets. Undefined metrics remain undefined.
6. A prediction that improves agreement is not sufficient evidence to accept a biochemical correction. Require a traceable, organism-appropriate source for reaction chemistry, gene assignment, direction or biomass composition. Keep provisional hypotheses separate from established corrections. Do not fit biomass coefficients to the benchmark score.

## Procedure for a future independent evaluation

The next validation needs an external evaluator or custodian and an explicit exposure record. The evaluator should choose eligible organism-level datasets before seeing comparative model results, document their source versions and experimental design, and retain outcome access until the method is frozen. A genuinely new dataset is preferable if inspection history of existing downloads cannot be established. Whole organisms should be kept together: randomly splitting gene-condition rows from an organism used in model curation would not test organism-level transfer.

Before outcomes are revealed, publish an immutable study record containing the eligible set and exclusions, exact baseline and intervention code, reconstruction and annotation versions, model-generation rules, media mapping procedure, all thresholds, intended comparisons and missing-data handling. Document any literature or curated-model evidence used for those organisms: phenotype findings embedded in that evidence are another potential source of exposure. Store outcome files separately from method development and record who accessed them. A file hash verifies bytes; it does not enforce that separation or certify a declaration.

Evaluate all eligible organisms, not only those on which the model grows or the intervention helps. The primary report should show each organism's paired MCC difference on common observations and its change in wild-type growth coverage. Show the distribution of those organism-level effects and its uncertainty; choose any aggregate statistic and resampling scheme before unblinding, once the size and grouping of the independent panel are known. Do not let a large gene-condition matrix from one organism silently dominate the conclusion.

Freeze all annotation-driven correction rules before unblinding. If a rule cannot resolve an organism's case without its phenotype outcomes, record an abstention. Adjudication after looking at validation outcomes is a new development iteration and requires another independent evaluation. A failed independent test is a result to report.

## Reusable file freeze

`scripts/freeze_study.py` records explicitly named regular files, their SHA256 hashes and sizes, the normalized plan, a UTC timestamp and the Git working-tree snapshot. It verifies every named file and manifest content consistency, rejects unsafe paths and refuses overwriting manifests. It does not freeze unnamed files, sign a manifest, authenticate a timestamp or establish data independence.

```sh
.venv/bin/python scripts/freeze_study.py freeze --root . \
  --plan data/studies/development_method_v1.json \
  --output results/study_freezes/development_method_v1.json
.venv/bin/python scripts/freeze_study.py verify --root . \
  --manifest results/study_freezes/development_method_v1.json
```

The development-method manifest packages the evaluation procedure and input exposure record. Its timestamp may follow some development results. For the actual time ordering of each simulation, consult that experiment's own `prespecified_manifest.json` and the preserved attempt records. Create a new version when named files change; do not overwrite an old freeze to make it pass.
