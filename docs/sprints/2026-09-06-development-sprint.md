# Development sprint: completion of cycles 6 and 7

All eight specified runs completed successfully. Cycle 7 preserved wild-type growth coverage and improved matched gene-level MCC in three organisms; B. thetaiotaomicron was unchanged. These are retrospective development results on phenotypes already used during curation. The freeze below records the configurations kept unchanged throughout this completion exercise; it does not turn these data into an independent test.

## Configuration fixed before simulation

The [manifest](../../results/development_sprint_2026_09_06/prespecified_manifest.json) was written at **2026-09-06 10:21:58 UTC**, before the first run began at 10:22:26 UTC. Its SHA-256 is `393352467b6dfee213d246c0196e04421854da01decd211f16128be0d9a67809`. It records eight commands, 50 input-file hashes, parameters, and dependency versions. The pipeline is based on audited commit `e4193ad`; the new execution controller is separately fingerprinted in the manifest.

Each organism used the existing gap-filled draft, universe patches v0.1, GPR patches v0.2/v0.3/v0.4, medium completion, retained rich-medium-essential genes, GLPK, two deletion workers, and all mapped conditions. Cycle 6 used model patches v0.3; cycle 7 used v0.4. Runs executed sequentially. No patch contents, statuses, or gene rules were changed or selected from these outcomes.

## Completed runs

MCC below uses each arm's own scored genes and the conditions where its wild type grows. Coverage is growing/mapped conditions; unmapped experimental conditions remain outside this denominator.

| Organism | Cycle | Scored genes | WT growth | Native MCC | Runtime, seconds |
|---|---:|---:|---:|---:|---:|
| Btheta | 6 | 510 | 16/25 | 0.57310 | 69.3 |
| Btheta | 7 | 510 | 16/25 | 0.57310 | 71.8 |
| Putida | 6 | 1,047 | 34/43 | 0.56137 | 216.7 |
| Putida | 7 | 1,047 | 34/43 | 0.58716 | 198.3 |
| MR1 | 6 | 713 | 11/12 | 0.52916 | 63.1 |
| MR1 | 7 | 719 | 11/12 | 0.53648 | 63.4 |
| Smeli | 6 | 989 | 27/33 | 0.67458 | 164.3 |
| Smeli | 7 | 989 | 27/33 | 0.68944 | 161.0 |

This completes the previously unfinished Smeli cycle 6 arm: 27 of 33 mapped conditions grow. Total subprocess runtime was 1,007.9 seconds. Every prepared model passed the enforced energy-cycle gate; no solver or run failures occurred.

## Comparison on shared observations

The prespecified comparison aligns Browser genes and conditions, retains conditions where both wild types grow, and uses the same finite observations in both arms. Confidence intervals are paired percentile bootstraps over genes, with 500 resamples and seed 0; conditions and the earlier patch-selection process are fixed.

| Organism | Shared pairs | Matched cycle 6 MCC | Matched cycle 7 MCC | Difference, 95% interval |
|---|---:|---:|---:|---:|
| Btheta | 8,160 | 0.57310 | 0.57310 | 0.00000 [0.00000, 0.00000] |
| Putida | 35,598 | 0.56137 | 0.58716 | +0.02579 [+0.00626, +0.05146] |
| MR1 | 7,843 | 0.52916 | 0.54466 | +0.01550 [0.00000, +0.03441] |
| Smeli | 26,703 | 0.67458 | 0.68944 | +0.01486 [0.00000, +0.03767] |

There were no gained or lost growing conditions in any cycle 7 comparison. MR1's matched and native cycle 7 MCC differ because six newly represented genes are excluded from the matched comparison and included in the native score.

## Descriptive inspection, without patch adjudication

- **Btheta:** v0.4 contains no organism-specific changes, and no growth classification changed.
- **Putida:** removing the menaquinol biomass requirement and the OXCDC/PHPYROX gap-fill reactions changed 135 shared predictions to growth, all agreeing with the operational fitness threshold. These involve PP_1658, PP_2458, PP_3232, and PP_0665. Agreement does not establish that removing quinone demand without replacement is biochemically correct; quinone composition and demand require separate evidence and sensitivity analysis.
- **MR1:** replacing the ornithine, lactate, and zinc gap-fill routes changed 33 shared predictions: 28 agreed and five disagreed with the fitness threshold. Three adverse changes are arginine-pathway predictions on casamino acids; the other two involve SO4245 on glutamine and D,L-lactate. Six added genes produce 66 newly scored observations, of which 27 disagree. In particular, SO0565 and SO0566 are predicted necessary in all 11 growing conditions despite fitness values above −2 throughout. These new observations are additional coverage, not comparisons against previous gene predictions. Source-based review should resolve zinc-transporter identity, substrate specificity, redundancy, and medium assumptions before any rule revision.
- **Smeli:** the sole model-patch difference is removal of menaquinol from biomass. It changes 54 shared predictions to growth; 53 agree with the threshold. The adverse change is SM_b21107 on L-fucose, where measured fitness is −3.477. This local disagreement is retained alongside the overall score increase.

Here “agreement” uses fitness < −2 as an operational low-fitness label in a pooled mutant assay; it does not establish gene lethality or a causal biochemical mechanism. The [complete descriptive attribution](../../results/development_sprint_2026_09_06/descriptive_change_attribution.json) records every newly scored observation, changed gene, adverse change count, and remaining no-growth condition. No rules were accepted, rejected, or repaired from this inspection.

## Reproduction and checks

For a fresh reproduction from the repository root, run `.venv/bin/python scripts/run_development_sprint.py --output-dir /path/to/new-output`. The controller creates a new manifest for the current runtime. Reusing an existing output directory verifies its frozen inputs and Python/dependency environment, then reuses only complete results with matching parameters and artifact hashes; it refuses an environment mismatch. Failed or partial attempts are preserved. A resume check in the original runtime validated all eight completed runs without repeating simulations or overwriting artifacts.

The independent verifier recomputed saved MCC and Bernstein AUC-PR, checked condition-table consistency, and spot-checked fitness against the raw Browser tables: **zero problems in all eight runs**. Details are in [independent_verification.json](../../results/development_sprint_2026_09_06/independent_verification.json). Full scores, timings, paired comparisons, and paths to each card and matrix are in [summary.json](../../results/development_sprint_2026_09_06/summary.json).

The next bounded step is source-driven review of cofactor biosynthesis/demand and transporter assignments under the separate study freeze, followed by an independent evaluation design. These exposed phenotype outcomes should remain development diagnostics rather than a rule-selection test presented as validation.
