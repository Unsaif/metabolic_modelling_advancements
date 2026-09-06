# PpnP repair sprint — 2026-09-06

The provisional nucleoside phosphorolysis intervention removes the model's ribokinase bottleneck while retaining quinone synthesis and its biomass demand. In the exposed Putida development benchmark, every scored candidate variant corrects the same **34 PP_2458/ribokinase predictions**. These are 34 of the 65 disagreements introduced by the preceding quinone-demand repair; the 31 PP_5317/chorismate-lyase disagreements remain. This is a mechanistically explained development improvement, not independent validation or a confirmed biological repair.

The [study definition](../studies/ppnp-repair-v1.md) declared seven physical arms and four benchmark runs before simulation. The hypothesis arose from previous development errors. Historical patches remain unchanged.

## Complete benchmark results

All four runs score **1,047 genes × 34 growing conditions = 35,598 finite observations**. WT grows in 34/43 mapped conditions. These matrices have no missing fitness or nonfinite simulations. Thresholds remain growth 0.001 and fitness −2, with medium completion and rich-medium observations retained. Source: [final summary](../../results/ppnp_repair_2026_09_06/runs/main/summary.json).

| Arm | MCC | Balanced accuracy | Accuracy | Bernstein protocol PR score | Standard AP | Standard AUROC |
|---|---:|---:|---:|---:|---:|---:|
| Fresh parent | 0.574384008 | 0.763592090 | 0.956963874 | 0.499391718 | 0.401587529 | 0.761209472 |
| Both forward | 0.580981027 | 0.764098963 | 0.957918984 | 0.508115277 | 0.420933930 | 0.761837317 |
| Both forward, RHCYS closed | 0.580981027 | 0.764098963 | 0.957918984 | 0.508115277 | 0.411079005 | 0.761780205 |
| Both reversible | 0.580981027 | 0.764098963 | 0.957918984 | 0.508115277 | 0.414522442 | 0.761716395 |

With growth positive, parent confusion counts are TP/TN/FP/FN = 32,944/1,122/937/595; candidates give 32,978/1,122/937/561. All changes improve the binary comparison. The [102 change records](../../results/ppnp_repair_2026_09_06/runs/main/changed_predictions.json) repeat the same 34 corrections across three arms, not 102 independent corrections. Their gene-condition identities match the preceding study's 34 ribokinase disagreements.

Each candidate's **paired MCC difference is +0.006597020**, with 95% gene-bootstrap interval **[0, 0.022116647]** from 500 resamples, seed 0. The effect concerns one gene, not 34 independent gene-level replications. Conditions are fixed; adaptive selection is unaccounted for. The interval includes zero and does not establish independent predictive progress.

Individual card MCC intervals—parent [0.473265006, 0.663850618], candidates [0.478804626, 0.675308072]—describe each score, not its change. Other metric intervals remain in the cards and summary.

Adding PP_4248 increases model genes from 1,305 to 1,306, but adds **no scored gene**. PP_4248/PpnP and PP_4976/AhcY lack exported fitness rows. The export contains 4,778 gene rows against 5,661 annotated genes. Missing exported rows do not establish absence from a mutant library, biological essentiality or dispensability. In particular, the benchmark cannot directly validate the proposed PpnP GPR or the AhcY dependence.

## Physical explanation and controls

All seven physical arms grow in glucose and permit quinone production. The parent requires RBK; adding either forward adenosine or forward inosine phosphorolysis independently permits RBK-deletion growth near 0.91926. The principal forward arm gives 0.919264794 in that deletion context. Removing PpnP together with RBK, or closing the added reactions while RBK is disabled, returns numerical zero. RBK/AhcY deletion also blocks rescue. RBK/ADA deletion blocks the inosine-only arm but not an adenosine-containing arm. These are conditional model predictions.

The [independent balance analysis](../../results/ppnp_repair_2026_09_06/analysis/mechanism_interpretation.md) shows how three obligatory quinone methylations produce SAH. With demand ε per unit growth:

`v_RBK + v_PUNP1 + v_PUNP5 ≥ 3ε × Growth`.

Nested balances explain the AhcY/ADA controls. PpnP produces ribose-1-phosphate from nucleosides, bypassing free ribose. Saved optimal witnesses are not unique or measured fluxes.

PPM deletion does not remove rescue. A saved witness uses the inherited reverse PNP → RNMK → NMNN route, whose exact net equation is:

`ribose-1-P + ATP + H2O → ribose-5-P + ADP + Pi + H+`.

RNMK maps to PP_4218; PNP and NMNN have no GPRs. This additional, incompletely assigned route prevents interpreting PPM-deletion survival as evidence for or against PP_1777's proposed phosphopentomutase activity. Closing questionable RHCYS chemistry leaves the candidate rescue intact. Permitting reverse PpnP flux produces no additional binary benchmark changes or material difference in the saved glucose growth controls.

All **35 primary energy-from-nothing tests**—five currencies in each of seven arms—report zero. They test ATP, NADH, NADPH, reduced Q8 and the periplasmic proton currency. This is a limited consistency check. ATP maintenance has lower bound zero throughout; reaction capacities and the inherited quinone demand are modelling assumptions, not measured enzyme capacities or physiological robustness guarantees.

## Source evidence and unresolved biology

Source-balanced CarveMe reactions PUNP1/PUNP5 receive a provisional PP_4248 GPR. [Evidence note 11](../evidence/11-putida-ppnp-candidate.md) verifies its exact 94-aa reference/UniProt sequence and HAMAP annotation through an E. coli template. This establishes identity, not KT2440 catalysis.

The [substrate audit](../../results/ppnp_repair_2026_09_06/evidence/substrate_evidence.md) supports E. coli forward phosphorolysis but finds contradictory reverse-direction statements in the Sévin supplement's Table 9 and Figure S8. The primary arm uses forward chemistry; the reversible sensitivity cannot settle that contradiction.

PP_1777's annotation supports phosphomannomutase, leaving ribose-phosphate activity uncertain. Upstream quinone assignments, tail representation, donor abstractions, demand, gene-less salvage and maintenance assumptions also remain provisional.

## Numerical reproducibility and verification

The fresh parent reproduces every binary growth prediction from the previous saved parent. WT growth, fitness values and matrix labels are identical; mutant growth differs by at most 6.002e-11. Nevertheless, raw standard AP changes from 0.406846334 to 0.401587529 because tiny flux differences can reorder tied or nearly tied predictions. The candidate AP differences in the table should therefore not be treated as clean evidence of biological improvement.

A separate post-outcome [rank-sensitivity diagnostic](../studies/growth-rank-numerics.md) uses a fixed precision grid without changing the primary cards. At a rounding quantum of 1e-9, the two parent AP values both become 0.414810130 and AUROC both becomes 0.761461358. Zero-only treatment leaves some rank sensitivity. None of the 24 transformed run/policy/precision combinations changes a binary prediction. These diagnostic point estimates have no recomputed bootstrap intervals and do not select a correct biological precision. A future ranking policy needs to be specified and frozen before evaluation.

The [independent artifact verification](../../results/ppnp_repair_2026_09_06/runs/main/independent_verification.json) passed: 147 input files, seven model inventories, four benchmark cards, three paired comparisons and all 102 change records were checked. The [independent solver verification](../../results/ppnp_repair_2026_09_06/runs/main/solver_verification/summary.json) also passed: 147 case records, four skips for unrepresented PpnP functions, and 286 GLPK/HiGHS solves. Maximum objective disagreement was 2.668e-11; maximum independently computed mass-balance residual was 8.067e-10 and bound violation 5.885e-11, both below 1e-8. That checker read no fitness values and did not repeat the energy probes. Numerical agreement does not establish exact optimality or biological validity.

The full software suite passed 239 tests plus six subtests, with 18 known single-label metric warnings. No frozen source, evidence, model, parameter or primary result was amended after observing outcomes.

## Freeze and next decision

The proposal was committed and pushed at **b5c540d45230cf542a4e35dff9f89fd1021df0a5** before simulation. The [manifest](../../results/ppnp_repair_2026_09_06/runs/main/manifest.json) freezes 86 inputs with content fingerprint:

`ec7dc1dfdb7eed71276dd3c84d33984340c44d68203bf18b5aeb46efecdb2ede`.

It records a dirty repository: no tracked changes, only the untracked independent checker `scripts/verify_ppnp_repair.py`. Runtime code and biological evidence were committed before simulation. The completed checker was separately snapshotted before verification; it is not among the 86 experiment inputs. Final verification confirms all 86 inputs unchanged. Runtime was 776.43 seconds.

The exact parent model SHA256 is `4ee05b6b43f87fc0baa248654013973179e0ac89153e1fa2948e443a1e4dd34d`; the frozen source universe SHA256 is `b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19`. Source records, sequence comparisons and assay-evidence provenance remain in the [evidence folder](../../results/ppnp_repair_2026_09_06/evidence/).

Next source work should investigate PNP/NMNN assignments, PP_1777 specificity and quinone assumptions before selecting another intervention. Biological validation requires separate PP_4248 and PP_1777 assays, distinguishing mutase from ATP-dependent salvage, and testing physiological contribution. Independent predictive evaluation requires independent data. This sprint provides a reproducible hypothesis and bounded development improvement, not demonstrated advancement of the field.
