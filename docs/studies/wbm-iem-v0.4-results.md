# Whole-body IEM protocol v0.4: Harvetta and flux ranges (results)

6 October 2026, Claude (Opus 5.5). Corrected the same day after independent verification (see the last section).

- **Plan:** [wbm-iem-v0.4-plan.md](wbm-iem-v0.4-plan.md), frozen 5 October at 13:39:37Z (fingerprint `f3dde4e2…`, manifest `results/study_freezes/wbm_iem_v0.4_plan.json`, external record in the Claude project). No Harvetta IEM result and no biomarker minimum existed at the freeze.
- **Analyses:** `scripts/v04_report.py` runs every declared analysis and applies the reading rules (`results/wbm_iem/v0.4_report.json`).
- **Verification:** an independent agent recomputed every number from the raw files with its own code (`results/wbm_iem/independent_verification_v0.4/REPORT.md`).

## Answer in one paragraph

**Every Harvey finding replicated on Harvetta, and the flux-range fix failed.**

- **Setup.** The port matches the COBRA Toolbox exactly on Harvetta too: all 83,521 bounds are identical after the setup and after the global constraints.
- **Accuracy.** Harvetta gets 218 of 251 directions right (86.9%), against 217 (86.5%) for Harvey.
- **Same calls as MATLAB.** MATLAB makes the same call for all 221 Harvetta biomarkers where it returned values.
- **Solver failures matter more on Harvetta.** Gurobi returned no optimum for 30 biomarkers (Harvey: 14), and one more is absent from the model, so `runIEM_HH`'s own figure falls to 75.8%. HiGHS solves all 30; Python is right on 27 of them.
- **Capped values.** Most "no change" errors (16 of 18) have the same positive maximum in health and disease, so something other than the blocked pathway limits them. Only 7 of the 16 equal a bound on the same metabolite.
- **Flux ranges.** Comparing minima as well as maxima recovers none of the capped calls on either model. The full-range rule loses 10 (Harvey) and 13 (Harvetta) correct calls.

## Runs

| Run | What | Fingerprint | IEMs | Wall time |
|---|---|---|---|---|
| A | Harvetta, maxima (the v0.3 protocol) | `51d7bda7…` | 57, all complete | 12.9 h |
| B | Harvey, minima | `92385e32…` (LP bounds hash `693af5c9…`, as v0.3) | 57, all complete | 6.1 h |
| C | Harvetta, minima | `5b6c8d07…` | 57, all complete | 6.8 h |
| M | MATLAB `runIEM_HH` on Harvetta (Gurobi 12.0.0) | job `iem_ref_harvetta_20261005_step2b_runiem` | 57 | 3.05 h |

Notes:
- **One process each.** A, B and C each ran as a single process; no shards were needed. `logs/watch_v04.sh` started C when B finished and would have resumed any killed run; none was killed. B finished on 5 October (19:48Z), A and C on 6 October (02:36Z and 02:38Z).
- **M took two attempts.** The first submission (`…_step2_runiem`, 13:40Z) ended at 14:34Z with an error from the queue: its Docker status check on Tim's Mac timed out. MATLAB itself reported no error. After inspection it was submitted once more under a new name, as the plan allows. That attempt ran from 15:28Z to 18:31Z. Both records are in `results/wbm_iem/matlab_reference/harvetta_step2/`.
- **Versions.** MATLAB R2024b Update 9, COBRA Toolbox 67c790d and Gurobi 12.0.0, as recorded by the job runner's own verification (`tools/matlab/runner/RUNNER_VERIFICATION_2026-10-05.md`). The job records themselves name only "gurobi".

## Q1: the Harvey findings on Harvetta

| Finding | Harvey (v0.3) | Harvetta (v0.4) | Reading |
|---|---|---|---|
| (a) Setup bounds identical to MATLAB | 81,094 / 81,094 | 83,521 / 83,521 (setup and global) | replicated |
| (b) Same call wherever MATLAB has values | 237 / 237 | 221 / 221 | replicated |
| (c) Biomarkers without a MATLAB optimum (one more is absent from the model) | 14 | 30 | replicated |
| `runIEM_HH` figure vs Python accuracy among scored | 81.0% vs 86.5% | 75.8% vs 86.9% | |
| (d) No-change errors with equal positive maxima | 18 of 21 | 16 of 18 | replicated |
| (e) Correct calls resting on a change below 5% | 10 of 217 | 5 of 218 | (descriptive) |

### Accuracy (A1) and error anatomy (A3)

| | Harvey v0.3 | Harvetta v0.4 |
|---|---|---|
| Correct / scored | 217 / 251 (86.5%) | 218 / 251 (86.9%) |
| Opposite direction | 13 | 15 |
| No change predicted | 21 (18 capped, 3 zero in both states) | 18 (16 capped, 2 zero) |
| IEMs fully correct | 38 | 35 |

The published figures are 84.9% for Harvey (214 of 252) and 85.3% for Harvetta (215 of 252), from other model and code versions (Thiele et al. 2020).

**Opposite-direction errors.** Twelve are shared with Harvey: citrulline in CPS1, NAGS and OTC, AADC L-DOPA and 3-methoxytyrosine, PC glucose, GMT creatinine, FED cholesterol, OXOP 5-oxoproline, LNS folate, NAGS orotate and HCYS ornithine. Harvey's thirteenth, MMA blood 3-hydroxypropionate, is right in Harvetta, on a change of 0.3%. Harvetta adds three:
- **ASNSD blood asparagine.** The healthy maximum is essentially zero in Harvetta (1.3 × 10⁻⁶ in Python; MATLAB prints 0) but 130.7 mmol/day in Harvey (Python only; MATLAB returned no optimum there), so the call turns from "decreased" to "increased". This is a difference between the models, not a threshold effect.
- **MMA blood carnitine.** In Harvey both states sat at 50 mmol/day, the dietary carnitine bound (a no-change error). In Harvetta the healthy value is 32.8, below it, so the call becomes "increased".
- **BTD urinary 3-hydroxypropionate.** In Harvetta its healthy maximum is 1.85 × 10⁻⁶, just above the 10⁻⁶ threshold below which the protocol sets a value to zero, so the call becomes "decreased". In Harvey both values were zero and the call was "no change". The difference between the two models is numerical noise at the threshold, not biology. MATLAB returned no healthy optimum for this biomarker on either model.

**Capped values.** Of the 18 capped no-change errors in Harvey and the 16 in Harvetta, 6 and 7 equal a bound on the same metabolite: kidney filtration limits (for example 2.592 mmol/day, 20 µM × the filtration rate) and, in Harvey, the dietary carnitine bound of 50 mmol/day for blood carnitine. The three MMA acylcarnitines at 50 are most likely limited by the same carnitine supply. The others match no single bound and were not traced. Ten of Harvetta's 16 capped errors are also capped errors in Harvey.

### Python against MATLAB on Harvetta (A2)

- **Bounds.** Identical after the setup (checked before the freeze) and after the global constraints.
- **Calls.** Same call on all 221 biomarkers where MATLAB returned values. All 30 differences are biomarkers that MATLAB recorded as NaN, so `runIEM_HH` scores them "unchanged"; Python calls 27 of these 30 correctly. Its three errors are LNS blood folate and PC blood glucose, which are also errors on Harvey where MATLAB has values and agrees, and BTD urinary 3-hydroxypropionate (the threshold case above). On Harvey, MATLAB also had no healthy optimum for that BTD biomarker, but both implementations called it "unchanged", so it is not among Harvey's 13 differences.
- **Values.** Above 10⁻³, all 291 pairs of values agree within 0.008% (Harvey: 293 pairs, all within 0.03% except FED blood cholesterol at 0.3%).
- **Where the failures are.**

| IEM | Biomarkers without a MATLAB optimum |
|---|---|
| PC | 9 |
| BTD | 6 |
| HPII | 5 |
| LNS | 3 |
| CIT1, XAN1 | 2 each |
| ADSL, ASNSD, DGK | 1 each |

The 31st non-finite entry is HPC `EX_25aics[u]`, which is absent from the model (marked NA).

- **The figures.** By `runIEM_HH`'s formula, MATLAB scores 191 of 252 (75.8%). On the 221 biomarkers with values it scores 191 (86.4%).

### Effect sizes (A4, declared in advance for Harvetta)

| Minimum relative change τ | Harvey v0.3 | Harvetta v0.4 |
|---|---|---|
| protocol (10⁻⁶ absolute) | 217 | 218 |
| 0.1% | 212 | 215 |
| 1% | 210 | 214 |
| 5% | 207 | 213 |
| 10% | 207 | 212 |

Harvetta's five correct calls below 5%:
- GACR blood glutamine (2.5 × 10⁻⁵);
- HPC coproporphyrin III in CSF and blood (1.8 × 10⁻⁴, 2.5 × 10⁻⁴);
- MMA blood 3-hydroxypropionate (0.3%);
- GACR urinary glutamate (2.4%).

## Q2: flux ranges (A5)

| Model | Maxima only | R1, full range | R2, tie-break | R1 transitions | Capped cases where the minimum differs |
|---|---|---|---|---|---|
| Harvey | 217 | 207 (−10) | 217 (0) | 10 right → wrong (all conflicting); 1 error → another; 0 fixed | 1 of 18 |
| Harvetta | 218 | 205 (−13) | 218 (0) | 13 right → wrong (all conflicting); 3 errors → another; 0 fixed | 3 of 16 |

By biofluid (declared in the plan): R1 lowers the correct calls in blood from 102 to 97 and in urine from 107 to 103 in Harvey, and from 105 to 98 (blood) and 105 to 100 (urine) in Harvetta. CSF (7 of 7) is unchanged, and the one other biomarker (brain, in ADSL) is lost under R1 on both models. R2 changes no biofluid's count.

**Readings (pre-declared).**
- **R1 is harmful:** it worsens on both models.
- **R2 is not useful:** it fixed no error on either model.
- **Errors turned into other errors.** Both rules turn one capped error on Harvey and three on Harvetta into "decreased" errors: BTD urinary 2-methylcitrate on both models, plus MMA blood C4-dicarboxylic acylcarnitine (`DM_c4dc[bc]`) and PC urinary acetone on Harvetta. In each, the healthy minimum is above zero and the disease minimum is zero.

**What the data show, and our interpretation.**
- **The capped biomarkers.** Their minimum is zero in both states in 17 of 18 cases in Harvey and 13 of 16 in Harvetta. In those cases the model can send the metabolite elsewhere, so a blocked pathway forces nothing into urine or blood.
- **The full-range rule's losses.** All are "conflicting": the disease maximum rises while the disease minimum falls to zero, and the healthy minimum is above zero in every case.
- **Interpretation (not tested).** The protocol's healthy state forces maximal flux through the IEM's reactions, which forces some downstream metabolites to be made and excreted; the disease state forces nothing. On this reading the healthy reference is an extreme state, and Shlomi-style flux-range logic does not carry over to it.

## Post hoc (not in the plan, descriptive)

**Reporting capped values as indeterminate.** A biomarker whose healthy and disease maxima are equal and positive cannot be called by a maximum-only comparison. Reporting such cases as "indeterminate" instead of "no change" would remove only errors:
- Harvey: 18 indeterminate; 217 of the remaining 233 calls correct (93.1%).
- Harvetta: 16 indeterminate; 218 of 235 correct (92.8%).

This is a reporting change, not an improvement of the model. Coverage falls from 251 to 233 and 235 biomarkers.

## Deviations

1. **MATLAB job M ran twice.** See Runs.
2. **The plan quotes Harvey's older GFR for Harvetta.** For Harvetta the older value is 128.64 ml/min (20% of renal plasma flow with the female parameters), not 129.75. The port check used the correct value; only the plan's text is affected.
3. **Analysis scripts.** Two were written after the freeze, as implementations of declared definitions:
   - `scripts/iem_error_anatomy.py` (A3) reproduces the Harvey v0.3 counts;
   - `scripts/v04_report.py` applies the frozen reading rules.
4. **The plan quotes the published accuracies the wrong way round.** Thiele et al. 2020 report 84.9% for Harvey (214 of 252) and 85.3% for Harvetta (215 of 252); the plan gives 85.3% for Harvey and 84.9% for Harvetta. They were context only, and no reading depends on them.

## Corrections after independent verification (6 October)

The verifier confirmed the counts, values, hashes, times and readings, and found these problems in the first version of this note and in the Paper 2 draft. All are corrected here and in the paper.

1. **Published accuracies swapped** (Deviation 4). The same swap was in Paper 2 and the roadmap.
2. **MATLAB failures.** An earlier draft said Gurobi returned no optimum for 31 (Harvetta) and 15 (Harvey) biomarkers; the 31st and 15th are the biomarker absent from the model.
3. **Caps overstated.** "Capped" means equal positive maxima. Only 6 of 18 (Harvey) and 7 of 16 (Harvetta) equal a bound on the same metabolite. The v0.3 note's uracil example (84.11 mmol/day) matches no bound in the LP.
4. **R1 transitions incomplete.** R1 also turns 1 and 3 errors into other errors.
5. **ASNSD wording.** Python's Harvetta value is 1.3 × 10⁻⁶, not zero, and the Harvey value of 130.7 is Python's alone.
6. **Mechanism stated as fact.** The explanation of the flux-range losses is now marked as an interpretation.
7. **Smaller fixes:** a cross-reference ("below" → "above"), the source of the Gurobi version, and the wording of the first MATLAB attempt's failure.

The verifier could not check one cited figure, Shlomi et al.'s recall of 0.10 on a clinical database: the full text was blocked by rate limits and a captcha. It is on Tim's checklist.

## Limits

- **Not blind.** The IEM labels and the published figures were known. Harvetta shares structure, labels, diet and constraints with Harvey.
- **One completed MATLAB run per model, with one solver version.**
- **Direction only,** against the protocol's literature labels.
- **Caps not traced.** Most capped values are not matched to the constraint that limits them.

## Next

1. **Paper 2** is drafted (`docs/paper/paper2-draft.md`; the Claude Doc is the working copy).
2. **The cap problem needs another approach.** Options, each to be planned and frozen before testing:
   - report capped values as indeterminate;
   - define a typical healthy state rather than one with maximal IEM flux;
   - relax the binding constraint of the biomarker itself, after tracing which constraint binds.
