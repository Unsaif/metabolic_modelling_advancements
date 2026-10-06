# Hand-curated reference models: results (Paper 1, post hoc and descriptive)

6 October 2026.

- **Plan:** `transfer-v1-curated-references-plan.md` (commit 260b92f, before scoring).
- **Deviations:** `transfer-v1-curated-references-deviations.md`.
- **Numbers:** `results/transfer_v1/reference_models/comparison.json`, produced by `scripts/compare_reference_models.py`.
- **Checks:** `results/transfer_v1/reference_models/independent_check/` and `independent_check_curated/` (including
  `second_check/`).

All the organisms below were used in development, and their fitness data had been seen. None of this is a held-out
test.

## Models

| Organism | Curated model | Source | Outcome |
|---|---|---|---|
| *E. coli* (control) | iML1515 | Bernstein et al. validation repository | scored |
| *P. putida* KT2440 | iJN1463 | BiGG | scored |
| *S. oneidensis* MR-1 | iSO783 | BioModels MODEL1507180036 | scored |
| *B. thetaiotaomicron* VPI-5482 | iAH991 | Rebuilt from Supplementary Tables S10a/S10b, because no public file was found. The parse is exact; the plan's validation checks were partly met. | 0 of 25 conditions under the protocol. The exploratory run with vitamin B12 was scored. |
| *S. meliloti* 1021 | iGD1575 | Nature Communications Supplementary Data 6 | Not scored: makes ATP from nothing (≥4 cycles) |

## Shared genes and conditions (MCC)

| Organism | Genes × conditions | B0 | U′ | M | Curated | Curated − U′ [95%] |
|---|---|---|---|---|---|---|
| *P. putida* | 764 × 28 | 0.53 | 0.57 | 0.61 | 0.56 | −0.006 [−0.108, 0.103] |
| *S. oneidensis* | 440 × 8 | 0.48 | 0.50 | 0.48 | 0.55 | +0.051 [−0.028, 0.125] |
| *B. thetaiotaomicron* (with B12) | 369 × 14 | 0.49 | 0.50 | 0.54 | 0.52 | +0.014 [−0.073, 0.102] |

**Curated − B0:** +0.035 [−0.063, 0.140], +0.070 [−0.013, 0.149] and +0.026 [−0.060, 0.117].

**Curated − M:** −0.042 [−0.139, 0.054], +0.064 [−0.010, 0.143] and −0.013 [−0.095, 0.067].

Every interval includes zero.

## Shared errors and agreement (U′ vs curated, shared genes and conditions)

| | *P. putida* | *S. oneidensis* | *B. thetaiotaomicron* |
|---|---|---|---|
| Calls in agreement | 93.9% | 89.9% | 90.1% |
| Cohen's κ | 0.57 | 0.62 | 0.66 |
| "Important" calls made by both (Jaccard) | 0.43 | 0.51 | 0.56 |
| U′'s wrong calls also wrong in curated | 710 / 1,190 (60%) | 383 / 580 (66%) | 478 / 711 (67%) |
| Expected if independent | 85 (8.3×) | 89 (4.3×) | 104 (4.6×) |
| Shared errors that are both-missed important genes | 416 | 359 | 303 |

## Genes only the curated model contains

- **iJN1463:** 383 genes. 52 of 1,287 "important" calls are confirmed.
- **iSO783:** 127 genes. 16 of 43 confirmed.
- **iAH991:** 447 genes. 30 of 568 confirmed. 266 of these calls come from one choice: a capsule polysaccharide in
  the biomass, which needs 21 genes at once. None of those 266 is confirmed.

## Union of genes (four-way MCC)

| Organism | Curated | Draft arms |
|---|---|---|
| *P. putida* | 0.45 | 0.49 to 0.58 |
| *S. oneidensis* | 0.53 | 0.47 to 0.49 |
| *B. thetaiotaomicron* | 0.45 | 0.49 to 0.56 |

Only the curated − M differences exclude zero: *P. putida* −0.130 [−0.227, −0.039] and *B. thetaiotaomicron*
−0.112 [−0.187, −0.040].

## Coverage (mappable conditions with growth)

| Organism | Curated | Drafts | Notes |
|---|---|---|---|
| *P. putida* | 36 / 43 | 28 / 43 | |
| *S. oneidensis* | 9 / 12 | 8 (B0), 9 (U′, M) / 12 | Not the same conditions. iSO783 also grows on two unmapped dipeptides. |
| *B. thetaiotaomicron* | 0 / 25 as published; 24 / 25 with B12 | 14 / 25 | Giving the drafts B12 changes none of their calls. |

## Reading

- On the genes both models contain, there is no clear difference between the corrected drafts and the curated
  models in any of the three organisms. The intervals allow the curated model to be up to about 0.1 better.
- The errors are largely shared, mostly as important genes that both models call dispensable.
- The curated models grow in at least as many conditions; for iAH991 this needs B12.
- Their extra genes add mostly unconfirmed "important" calls.
- The reasons for the shared errors are hypotheses, untested here:
  - binary growth calls versus graded fitness;
  - no regulation of backup enzymes;
  - gene functions that are wrong or unknown in both kinds of model.

## Sensitivity (*B. thetaiotaomicron*)

| Variation | Curated − U′ |
|---|---|
| Without the arabinan condition (representation differs, about 34× the carbon) | +0.013 [−0.074, 0.100] |
| Without L-fucose and L-rhamnose (growth only 0.0012/h) | +0.023 [−0.063, 0.111] |
| Sulfide added | Changes no call |
