# Putida quinone repair experiment v1

6 September 2026. This is a development experiment on already exposed data. The source-derived proposal, biological evidence, intervention recipe and execution code are recorded before the new model simulations. This ordering makes the experiment traceable; it does not establish independent validation.

The earlier [quinone investigation](../evidence/08-quinone-biomass-audit.md) found that the prepared v0.4 draft can circulate its represented quinone pool but cannot synthesize it. Removing the biomass requirement conceals this gap. This experiment asks whether five consecutive missing reactions from the local iJN1463 reference connect the pool to nutrient-derived precursors, and how requiring synthesis changes the model's predictions.

## Fixed arms

| Arm | Five missing reactions | Biomass ubiquinol demand | Gene rules |
|---|---|---|---|
| `baseline` | No | None | Prepared v0.4 |
| `demand_only` | No | Original template amount | Prepared v0.4 |
| `curated_path_only` | Yes | None | Copied source rules for added reactions |
| `curated_path_template` | Yes | Original template amount | Copied source rules for added reactions |
| `curated_path_curated` | Yes | Local curated-model amount | Copied source rules for added reactions |
| `sequence_path_template` | Yes | Original template amount | Three explicitly declared sequence-supported hypotheses |

The source copy comprises `OHPHM`, `OMPHHX`, `OMBZLM`, `OMMBLHX` and `DMQMT`. The two demand coefficients come from existing biomass reactions and are not fitted. The sequence arm changes `OPHHX`, `OMMBLHX` and `OPHBDC` as described in [evidence brief 09](../evidence/09-putida-quinone-repair-sources.md). All six are provisional representations. The historical v0.4 patch remains unchanged.

Source stoichiometry and reaction bounds are retained. Existing gene identities are reused through the recorded aliases. Cytosolic compartment labels are normalized through an explicit mapping; existing target metabolite metadata is preserved. Expected source-versus-target charge differences are enumerated, while unexpected identity or chemistry differences stop the transfer. Source elemental/charge balance and the target's inherited metadata limitations are reported separately.

## Checks and interpretation

Each arm receives a fixed glucose-medium growth test, separate net quinone/quinol production probes, the [declared-pool certificate](cofactor-audit-method.md), energy-from-nothing checks and targeted deletions of 13 pathway/support loci. Existing bounds and other model constraints remain in force during production probing. The prepared draft's ATP maintenance lower bound is zero: retaining that bound does not test a positive maintenance demand. A missing or orphan gene is unrepresented, and deletion results in a nongrowing wild type cannot establish a gene requirement.

Five arms also receive the established carbon-fitness protocol. The sequence arm can omit a duplicate benchmark only after verifying identical wild-type linear programs and that every changed gene rule involves exclusively loci without exported fitness rows. Under those conditions all covered single-gene deletion problems are also identical. This equivalence does not validate the differing predictions for unscored genes.

Every scored arm is compared with baseline on common genes and finite observations in conditions where both wild types grow. Native growth coverage and gene coverage remain separate. No-growth arms have no interpretable gene-level MCC. All adverse and null outcomes are retained; no reaction, gene rule or coefficient is selected using these scores.

The missing terminal genes have no exported fitness measurements in this dataset. Missing measurements do not establish absence from the original mutant library or biological essentiality. See the [coverage audit](../../results/quinone_repair_2026_09_06/evidence/benchmark_coverage.md).

The experiment can establish that a model representation supplies its cofactor and passes the stated checks. It cannot establish native KT2440 quinone chain length, a physiological biomass coefficient, hydroxylase donor specificity, or accessory-gene requirements. The source's half-oxygen abstractions and an inherited upstream redox shortcut remain explicit limitations. Agreement with iJN1463 is source consistency, not independent experimental evidence.

## Reproduction

Use the audit environment and a fresh output directory:

```sh
.venv/bin/python scripts/run_quinone_repair.py --out results/quinone_repair_new
```

The runner expands [the recipe](../../data/studies/quinone_repair_v1.json) into a file-hash manifest before model preparation. It verifies the frozen inputs throughout execution and refuses existing output directories. `--physical-only` runs the physical checks without loading numeric fitness values; it also requires its own fresh directory. A failed attempt is preserved and any reviewed retry uses a new directory. Exported candidate models describe the prepared reconstruction; the benchmark cards record medium completion and simulation settings.
