# Development sprint: restoring quinone synthesis

6 September 2026. The five-reaction candidate restores net synthesis of the represented quinone pool and permits growth with either of two preselected biomass requirements. Restoring the requirement also produces **65 new disagreements with measured gene fitness**, repeated under both coefficients. The candidate fixes the demonstrated connectivity defect but is **not accepted as a biological correction**. Native KT2440 quinone chemistry and the proposed gene assignments remain incompletely validated.

## Question and provenance

The earlier [biomass investigation](../evidence/08-quinone-biomass-audit.md) showed that v0.4 permits quinone cycling while blocking net synthesis. Deleting its biomass quinone requirement allowed growth but left that defect in place. We traced the precursor route and proposed the contiguous missing tail from the local iJN1463 source: `OHPHM → OMPHHX → OMBZLM → OMMBLHX → DMQMT`.

The [biological evidence](../evidence/09-putida-quinone-repair-sources.md), [pathway map](../../results/quinone_repair_2026_09_06/evidence/quinone_pathway_mapping.md), [experiment definition](../studies/quinone-repair-v1.md), recipe, code and tests were committed and pushed as [`78e12fd`](https://github.com/Unsaif/metabolic_modelling_advancements/commit/78e12fdce19013205c3c7ca7f3c3c9a784bfd661) before these simulations. The [execution manifest](../../results/quinone_repair_2026_09_06/runs/main/manifest.json) was created at 11:29:46 UTC from that clean checkout and hashes 48 inputs. Its content fingerprint is `d4a72ef731a9fa0451e0b1c9a0463837872e07dc1007772ce7ed2d1b53820734`.

Putida was already exposed during earlier model development. Publishing the definition before this experiment documents ordering and preserves the hypotheses; it does not make these data independent validation. No reaction, gene rule or coefficient was selected using this experiment's fitness outcomes.

The candidate JSON's `provisional_not_applied` status describes its preserved extraction-time proposal. The experiment's per-arm `intervention.json` and exported `model.xml.gz` record the actual applications. These are explicit experimental variants, not a replacement of historical patch v0.4.

## Physical results

The fixed glucose condition uses carbon uptake of 10 mmol/gDW/h, the established medium-completion assumptions, and GLPK feasibility tolerance `1e-9`. ATP maintenance remains at its inherited lower bound of zero. Production tests maximize a separate outward demand at growth lower bound zero; they measure the model's feasible production capacity, not a measured physiological rate.

| Arm | Glucose growth, h⁻¹ | Maximum Q8 demand, mmol/gDW/h | Maximum Q8H2 demand, mmol/gDW/h |
|---|---:|---:|---:|
| Baseline, no quinone requirement | 0.919371 | 0 | 0 |
| Requirement only | 0 | 0 | 0 |
| Five-reaction pathway only | 0.919371 | 0.768551 | 0.765845 |
| Pathway + template requirement | 0.919264 | 0.768551 | 0.765845 |
| Pathway + curated requirement | 0.919124 | 0.768551 | 0.765845 |
| Sequence-informed rules + template requirement | 0.919264 | 0.768551 | 0.765845 |

The template and curated coefficients are respectively `−9.68418998275772e−5` and `−0.000223` mmol/gDW. They are inherited model coefficients, not fitted values or newly measured pool sizes. Production maxima are below the probe capacity of 1,000. All five tested energy currencies have zero energy-from-nothing flux in each of the six arms.

The declared `q8_c + q8h2_c` pool has no supply reaction in baseline or requirement-only. In the four pathway arms, `DMQMT` supplies the pool. This provides a structural explanation for the growth rescue. It also explains why supplying the pathway without a growth-associated demand can leave synthesis optional: in the pathway-only arm, the pool's steady-state balance forces zero net `DMQMT` flux. Cofactor recycling can still support the same growth optimum.

All five added reactions are atom/charge balanced using the source metadata. They also report zero imbalance with preserved target metadata, but the latter contains zero-charge placeholders. The transfer normalizes the source's `c` compartment label to the draft's `C_c`, preserves existing target metadata, and records the two expected charge differences for `amet_c` and `h_c`. These bookkeeping checks do not validate the source's simplified oxygenase turnover chemistry.

An [independent numerical check](../../results/quinone_repair_2026_09_06/runs/main/solver_verification/summary.json) reconstructed each physical medium and bound setting from the saved models and verified them against the reports. GLPK and HiGHS each solved 76 problems: six wild types, 12 production probes and 58 represented targeted deletions. Twenty absent/orphan gene cases were explicitly skipped. Across 152 solves, maximum solver-objective disagreement was `1.55e−13`, maximum mass-balance residual `2.27e−11`, and maximum bound violation `2.27e−13`. This verifies numerical agreement and primal feasibility for these ordinary LPs; it is not an exact optimality or biological proof. Frozen inputs remained unchanged.

## Fitness comparison: adverse results retained

| Scored arm | Growing conditions / mapped conditions | Gene-level MCC | Change from matched baseline |
|---|---:|---:|---:|
| Baseline | 34 / 43 | 0.587165 | — |
| Requirement only | 0 / 43 | Undefined | Undefined |
| Pathway only | 34 / 43 | 0.587165 | 0 |
| Pathway + template requirement | 34 / 43 | 0.574384 | −0.012781 |
| Pathway + curated requirement | 34 / 43 | 0.574384 | −0.012781 |

There are 57 eligible condition groups, of which 43 map to the protocol. Comparisons between the growing arms use the same 1,047 measured genes and 35,598 finite gene-condition observations across the same 34 growing conditions. The requirement-only comparison has zero shared growing conditions and zero scorable pairs. Pathway additions increase model genes from 1,301 to 1,305, but add no measured genes. The paired gene-bootstrap interval for either nonzero MCC change is `[−0.031905, 0]` (500 resamples, seed 0); this describes variability over these development genes with fixed conditions and does not establish independent performance.

Each requirement-plus-pathway arm changes 65 binary predictions from growth to no growth: **34 for `PP_2458` (ribokinase, `RBK`) and 31 for `PP_5317` (chorismate lyase, `CHRPL`)**. Every associated measured fitness is at least `−2`, so all 65 changes worsen agreement under the fixed threshold. The two arms duplicate these same 65 observations; the saved file has 130 records, not 130 independent biological disagreements. `PP_5317` still permits growth on the three mapped 4-hydroxybenzoate/coumarate conditions.

These results should not be explained away as missing data: these two genes have measurements. They also do not prove the new pathway is wrong in every respect. They expose a conflict between the represented synthesis/disposal dependencies, the biomass requirement and the pooled fitness observations. Possible biological or assay explanations require evidence. No reaction, transporter, gene rule or demand coefficient was changed to recover the score.

The [independent artifact verifier](../../results/quinone_repair_2026_09_06/runs/main/independent_verification_portable.json) reproduces all five cards, native metric points and intervals, paired comparisons and changed-prediction records. It also verifies raw replicate averaging, coverage, source transfer and all six model/physical inventories without using the production scoring functions. The main run completed in 771.6 seconds; all 48 frozen inputs verified unchanged afterward.

## Post-outcome dependency trace

The observed disagreements prompted a separate [structural trace and exact certificate](../../results/quinone_repair_2026_09_06/evidence/postrun_dependency_trace.md); it was not part of choosing the candidate. Both modeled routes for disposing of S-adenosylhomocysteine, produced by the three quinone methylations, converge on free ribose. Ribokinase is the only possible sink for the declared pool `ahcys_c + rhcys_c + adn_c + ins_c + rib__D_c` under the represented reaction directions. Existing ribose uptake cannot export cytosolic ribose.

This yields a checkable bound without an optimizer. Let `t` be flux through the five-reaction tail and `ε` the positive quinone biomass requirement. The four isolated new intermediates make all five tail fluxes equal; the quinone-pool balance gives `t = ε × growth`. Each of the three methylations supplies one unit to the ribosyl pool, whose only sink is `RBK`. Every other external contribution to that pool is nonnegative. Therefore:

`v_RBK ≥ 3t = 3ε × growth`.

Deleting `RBK` forces zero growth in this representation for either positive coefficient. This explains a model dependency; it does not establish that ribokinase is biologically essential. The practical next question is whether the strain has an omitted disposal route or uses a different represented salvage mechanism. Any proposed correction needs biochemical evidence. For `CHRPL`, external 4-hydroxybenzoate and the represented coumarate degradation route bypass its precursor-production step in the three conditions that retain growth. Pooled-assay cross-feeding and alternative native precursor pathways remain untested explanations for the other disagreements.

## Gene predictions and missing measurements

The exact-source and sequence-informed representations give identical wild-type linear programs and identical problems for every benchmark-covered single-gene deletion. The [equivalence record](../../results/quinone_repair_2026_09_06/runs/main/sequence_benchmark_equivalence.json) verifies this from reaction constraints and the affected gene identities. The sequence arm therefore omits a redundant full benchmark, while retaining its separate targeted predictions.

Under the restored template requirement, the following glucose predictions differ:

| Locus | Exact-source representation | Sequence-informed representation |
|---|---|---|
| `PP_0548` / UbiX candidate | Gene unrepresented | Deletion blocks growth |
| `PP_5197` / UbiI candidate | Deletion permits growth | Deletion blocks growth |
| `PP_0427` / Coq7 candidate | Deletion permits growth | Deletion blocks growth |
| `PP_5013` / UbiB | Deletion blocks growth through inherited catalytic rule | Accessory function unrepresented |

“Unrepresented” is not a prediction of biological dispensability. UbiB's orphan gene object in the sequence arm is explicitly marked as having no represented function. UbiJ (`PP_5012`) is also unrepresented. The UbiD/UbiX AND rule expresses a hypothesized cofactor dependency, not a heteromeric catalyst.

Only `PP_5317` has an exported fitness row among the 12 immediate pathway/prenyl-support genes in the source model. All five distinct genes assigned to the missing tail lack rows. The [coverage audit](../../results/quinone_repair_2026_09_06/evidence/benchmark_coverage.md) distinguishes absent measurements from absent mutant-library members; the local files cannot establish why these rows are missing. The gene predictions above are testable hypotheses, not experimentally confirmed improvements.

## Limits and next discriminating evidence

This is a provisional repair of the represented pathway. The source uses half-oxygen hydroxylation equations without explicit reducing-donor costs. An inherited upstream `MECDPDH` reaction also omits a reducing donor and has charge imbalance under curated metabolite charges; the [upstream analysis](../../results/quinone_repair_2026_09_06/evidence/upstream_redox_uncertainty.json) preserves the alternatives. Passing the tested energy-cycle checks does not resolve those omissions or prove full thermodynamic consistency.

Native KT2440 quinone tail length, quantitative demand, hydroxylase donor specificity and accessory requirements remain unresolved. Useful independent evidence would include native quinone profiling and accession-specific hydroxylase/complementation experiments, as outlined in evidence brief 09. Adjusting these assumptions to maximize the existing benchmark score would not resolve the biology.

The reusable [cofactor checks](../studies/cofactor-audit-method.md) provide explicit conservation certificates and constrained production probes with failed-solve handling. They apply established conservation and demand-testing methods; this sprint does not claim new conservation theory.

## Validation and reproduction

The final repository suite passes **187 tests plus six subtests**. The 18 warnings are the existing single-label metric-test warnings. The new tests cover cofactor conservation/production edge cases, transfer identity and metadata guards, missing/orphan genes, freeze ordering, unsupported solver features, independent metrics, relocation and tamper detection.

Use `requirements-audit.txt` for the recorded environment. The main execution used Python 3.14.3, COBRApy 0.32.1, GLPK 5.0 and HiGHS 1.15.1. The [execution log](../../results/quinone_repair_2026_09_06/runs/main/execution.txt), candidate models, matrices, gene maps, physical contexts and verification reports are preserved.

The artifact verifier initially compared absolute checkout paths, which prevented checking committed results after relocation. Its original source and successful report are preserved. The corrected verifier checks arm/filename, model identity and actual file hashes while retaining the original path as historical provenance; a [complete relocated copy also verifies](../../results/quinone_repair_2026_09_06/evidence/independent_verification_relocated.json). No experimental input or output was changed for this correction.

Verifiers require fresh report paths and do not overwrite the recorded run:

```sh
.venv/bin/python -m pytest -q
.venv/bin/python scripts/run_quinone_repair.py --out results/quinone_repair_new
.venv/bin/python scripts/verify_quinone_repair.py --run-dir results/quinone_repair_new
.venv/bin/python scripts/verify_quinone_repair_solvers.py --study results/quinone_repair_new --out results/quinone_repair_new/solver_verification
```

## Source attribution

The candidate models retain their source layers and modification records. The [BiGG model page](https://bigg.ucsd.edu/models/iJN1463) links the UC Regents educational/research/non-profit [terms](https://bigg.ucsd.edu/license), which differ from the former README's CC-BY-SA claim. The required notice is now preserved in [models/bigg/NOTICE](../../models/bigg/NOTICE) and the [candidate model attribution](../../results/quinone_repair_2026_09_06/MODEL_SOURCES.md). No model bytes or frozen inputs were changed to add this attribution.
