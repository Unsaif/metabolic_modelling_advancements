# Whole-body IEM protocol v0.3: results and MATLAB comparison

5 October 2026, Claude (Opus 5.5).

- **Plan:** [wbm-iem-v0.3-plan.md](wbm-iem-v0.3-plan.md), frozen 4 October 02:06:26Z, fingerprint `66bd36b0…`. After a workspace reset the manifest was re-created with the same fingerprint (`results/study_freezes/wbm_iem_v0.3_plan_recovery.json`).
- **MATLAB reference run:** made on Tim's Mac through the `matlab-agent-runner` queue: MATLAB R2024b, COBRA Toolbox commit 67c790d, Gurobi 12.
- **Independent verification:** `results/wbm_iem/independent_verification_v0.3/REPORT.md`.

## Answer in one paragraph

**The Python port reproduces the COBRA Toolbox protocol.**
- The model setup is identical, bit for bit, in all 81,094 lower and upper bounds.
- Wherever both MATLAB/Gurobi and Python/HiGHS return an optimum, the direction call agrees for every biomarker (237 of 237). Values agree within 0.03%, except one at 0.3%.

**Accuracy.** Run the way the Toolbox runs it today (v0.3), **217 of 251 scored biomarker directions are correct (86.5%)**. That is close to the published 85%.

**The MATLAB run scores 204 of 252 (81.0%) by `runIEM_HH`'s own accuracy formula.** The whole gap comes from 13 biomarkers for which Gurobi returned no optimum: `runIEM_HH` silently counts these as "no change". HiGHS solves all 13 problems (feasibility certified twice, independently) and calls all 13 correctly. On the 237 biomarkers where MATLAB returned values, both implementations get the same 204 right (86.1%).

**Effect of the two setup changes:**
- Re-applying the current Toolbox constraints (1,572 bound changes) changes no direction call.
- Correcting v0.2's bile-duct step changes 5 calls (net −3), all of them decided by differences under 1%.

**Most "no change" errors are capped values.** In 18 of the 21 errors where the model predicts no change, healthy and disease reach the same maximum, a cap set by the physiological constraints.

## Runs

| Run | Model bounds | Bile-duct step | Correct / scored | Accuracy | Errors: opposite / no change | IEMs fully correct |
|---|---|---|---|---|---|---|
| v0.2 (3 Oct) | as shipped | all 261 exits (deviation) | 220 / 251 | 87.6% | 13 / 18 | 37 |
| v0.2b | as shipped | Toolbox list (28) | 217 / 251 | 86.5% | 13 / 21 | 38 |
| **v0.3 (primary)** | current Toolbox constraints re-applied | Toolbox list (28) | **217 / 251** | **86.5%** | 13 / 21 | 38 |
| MATLAB reference | `runIEM_HH` unchanged (same setup as v0.3) | Toolbox list | 204 / 237 with values | 86.1% | — | — |
| | | | 204 / 252 by `runIEM_HH`'s formula | 81.0% | | |

Notes on the table:
- **"Scored"** means an expected direction and both optima available. The unscored biomarker is HPC `EX_25aics[u]`, which the model lacks.
- **"IEMs fully correct"** means every protocol biomarker of that IEM was scored and correct.
- **`runIEM_HH`'s formula** divides by all biomarkers, so missing optima and absent biomarkers count as wrong.

Files:
- `results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json` and `..._v0.2b.json`, with their summaries;
- `Harvey_1_03d_iem_comparison_v0.3.{json,tsv}` (per biomarker, all three runs);
- `matlab_reference/matlab_vs_python_v0.3.json`.

## Python against MATLAB

**Bounds.** The setup bounds and the bounds after the global constraints are identical for all 81,094 reactions (`results/wbm_iem/matlab_reference/step1_setup_bounds_comparison.json`). The verifier hashed MATLAB's bounds and obtained the LP-bound hash recorded by the v0.3 run (`693af5c9…`).

**Calls.** Python and MATLAB make the same call for 238 of 251 biomarkers. Each of the 13 differences is a biomarker whose healthy optimum MATLAB recorded as NaN:

| IEM | Biomarkers |
|---|---|
| ASNSD | asparagine, blood and CSF |
| BTD | acetoacetate, acetone, 3-hydroxybutyrate, urinary ammonium |
| EF | urinary fructose |
| GMT | urinary urate |
| PHOX1 | urinary oxalate, glycolate, glyoxylate |
| SUCLA | lactate, pyruvate |

In all 13 cases the Python call is the expected one. MATLAB, scoring NaN as "unchanged", gets all 13 wrong. A 14th NaN (BTD `EX_3hpp[u]`) gives "unchanged" in both implementations.

**The 13 problems are feasible.** Two independent checks:
- `scripts/certify_iem_solutions.py` re-solved them with HiGHS and checked every vector against the LP HiGHS held (row and bound violations ≤ 2e-7, objective recomputed).
- The verifier built the LP itself from the model file and MATLAB's bounds. It obtained the same optima within 2.6e-5 and, with MATLAB's own disease values, the same calls.

We did not diagnose why Gurobi returned no optimum. The Toolbox writes NaN whenever the solver does not report an optimal solution.

**Values.** Among values above 1e-3, one differs by more than 0.1%: FED `DM_chsterol[bc]` healthy, 337.756 against 336.750. The next largest difference is 0.03%.

## What changed the calls

**v0.2 → v0.2b (Toolbox bile-duct list instead of all exits): 5 calls change.** All five rested on differences under 0.5%:

| IEM | Biomarker | v0.2 healthy → disease | v0.2b | Expected | Effect |
|---|---|---|---|---|---|
| DPYR | `EX_ura[u]`, `EX_56dura[u]` | 120.28 → 120.85 (increased) | 120.85 = 120.85 (unchanged) | increased | right → wrong |
| HCYS | `DM_Lhcystin[bc]`, `DM_hcys_L[bc]` | 28.18 → 28.27, 56.37 → 56.53 (increased) | equal (unchanged) | increased | right → wrong |
| HYPRO1 | `EX_4hpro_LT[u]` | 28.22 = 28.22 (unchanged) | 0 → 28.22 (increased) | increased | wrong → right |

**v0.2b → v0.3 (current Toolbox constraints re-applied): no call changes.** The re-application covers:
- kidney filtration from a GFR of 90 instead of 129.75 ml/min;
- CSF export bounds from 0.52 instead of 0.35 ml/min;
- three blood–brain-barrier uptakes;
- 31 AGORA diet uptakes;
- five bile acids added to the diet.

Values move, though. For example, urinary maxima limited by filtration scale by 0.69.

## The remaining errors (v0.3)

**Opposite direction (13).** Several are known physiology the model cannot reproduce:

| IEM(s) | Biomarker | Model prediction | Expected |
|---|---|---|---|
| CPS1, NAGS, OTC | blood citrulline | rises | falls |
| AADC | urinary L-DOPA and 3-methoxytyrosine | fall to zero | rise |
| PC | blood glucose | rises | falls |
| GMT | blood creatinine | falls to zero | rises |
| FED | blood cholesterol | rises | falls |
| OXOP | blood 5-oxoproline | rises | falls |
| LNS | blood folate | rises | falls |
| NAGS | urinary orotate | rises | falls |
| HCYS | blood ornithine | falls slightly | rises |
| MMA | blood 3-hydroxypropionate | falls slightly | rises |

**No change predicted (21).** In 18 of these, the healthy and disease maxima are equal and positive. Each such value is a cap:
- the urinary filtration limit for a metabolite without a measured blood range (20 µM × GFR = 2.592 mmol/day);
- a measured filtration or excretion limit (uracil 84.11);
- most likely the carnitine supply in the diet (MMA acylcarnitines at 50 mmol/day, the carnitine uptake limit set by `setDietConstraints`).

The protocol maximises each biomarker, so once both states reach the cap it cannot see a difference. The other three are zero in both states (CYP21D aldosterone and cortisol; BTD urinary 3-hydroxypropionate).

## Post hoc: how much rests on tiny differences

This analysis was not in the plan. It recomputes the calls with a minimum relative change τ between healthy and disease maxima (`results/wbm_iem/Harvey_1_03d_iem_effect_size_v0.3.json`):

| τ | v0.2 correct | v0.3 correct |
|---|---|---|
| protocol (1e-6 absolute) | 220 (87.6%) | 217 (86.5%) |
| 0.1% | 215 | 212 |
| 1% | 209 | 210 (83.7%) |
| 5% | 204 | 207 (82.5%) |
| 10% | 204 | 207 |

Ten of v0.3's 217 correct calls rest on changes below 5%, and five below 0.1%. Those include HPC coproporphyrin III (1e-7 relative) and GACR blood glutamine. The bile-duct correction showed how easily such calls flip.

## Corrections to the record

1. **v0.2 was not "without constraints".** Harvey 1.03d ships with the diet and physiological constraints applied. This is documented in the plan and corrected in the Sprint 5 note.
2. **v0.2 had a bile-duct deviation.** It opened all 261 bile-duct exits instead of the Toolbox's 28. v0.2's 220 should be read as 217 under the Toolbox protocol.
3. **Relation to the published 85%.** The Toolbox's own scoring of the same protocol gives 81.0% here, because solver failures are counted as "no change". The figure therefore depends on solver behaviour, not only on the model.

## Verification

An independent agent recomputed every number with its own code (`results/wbm_iem/independent_verification_v0.3/REPORT.md`):
- **Confirmed:** six of seven claims, exactly.
- **Seventh claim (certified optima):** holds for 12 of 13 as first stated. A fresh re-solve of BTD `EX_nh4[u]` differs from the v0.3 value by 8.8e-5 (2e-6 relative); the verifier's own solve is within 4e-6 of v0.3, and the call is unchanged.
- **Further checks:**
  - the freeze verifies;
  - every record carries its run's fingerprint;
  - the records saved at the mid-run checkpoint are unchanged in the final files;
  - the 1,572 setup changes all fall in the documented categories.

## Limits

- **Scope:** one model (Harvey 1.03d, male) and directions only. The ground truth is the protocol's literature labels.
- **Single MATLAB run:** the reference is one run with Gurobi 12. The published 2020 analysis used other model and code versions.
- **Restarts:** both Python runs were restarted after workspace resets and resumed under unchanged fingerprints. The original 02:06Z freeze record survives only in the Claude project; the re-created manifest has the same content fingerprint.
- **Post hoc analysis:** the effect-size analysis is descriptive.

## Next

1. **Report two points to the Toolbox maintainers:**
   - `runIEM_HH` scores solver failures as "no change";
   - the shipped Harvey 1.03d carries older constraint parameters (GFR, CSF flow) than the current code re-applies.
2. **Harvetta**, with the same pipeline. It runs in Python, and in MATLAB through the queue.
3. **Make the predictions robust:**
   - report solver failures separately;
   - add a minimum effect size;
   - address the saturation at caps that causes 18 of the "no change" errors, for example by comparing flux ranges or minimum excretion as well as maxima.
4. **Extend the ground truth** to the 279-pair biomarker table, with per-disease precision and recall: the road to Recon4IMD and ChatRD.
