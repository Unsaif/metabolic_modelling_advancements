# Whole-body IEM protocol v0.4: Harvetta and flux ranges (results)

6 October 2026, Claude (Opus 5.5).

- **Plan:** [wbm-iem-v0.4-plan.md](wbm-iem-v0.4-plan.md), frozen 5 October at 13:39:37Z (fingerprint `f3dde4e2…`, manifest `results/study_freezes/wbm_iem_v0.4_plan.json`, external record in the Claude project). No Harvetta IEM result and no biomarker minimum existed at the freeze.
- **Analyses:** `scripts/v04_report.py` runs every declared analysis and applies the reading rules (`results/wbm_iem/v0.4_report.json`).

## Answer in one paragraph

**Every Harvey finding replicated on Harvetta, and the flux-range fix failed.**

- **Setup.** The port matches the COBRA Toolbox exactly on Harvetta too: all 83,521 bounds are identical after the setup and after the global constraints.
- **Accuracy.** Harvetta gets 218 of 251 directions right (86.9%), against 217 (86.5%) for Harvey.
- **Same calls as MATLAB.** MATLAB makes the same call for all 221 Harvetta biomarkers where it returned values.
- **Solver failures matter more on Harvetta.** Gurobi returned no optimum for 31 of 252 biomarkers (Harvey: 15), so `runIEM_HH`'s own figure falls to 75.8%. HiGHS solves all of them; Python is right on 27 of the 30 scored ones.
- **Caps.** Most "no change" errors are capped values (16 of 18).
- **Flux ranges.** Comparing minima as well as maxima recovers none of the capped calls on either model. The full-range rule loses 10 (Harvey) and 13 (Harvetta) correct calls.

## Runs

| Run | What | Fingerprint | IEMs | Wall time |
|---|---|---|---|---|
| A | Harvetta, maxima (the v0.3 protocol) | `51d7bda7…` | 57, all complete | 12.9 h |
| B | Harvey, minima | `92385e32…` (LP bounds hash `693af5c9…`, as v0.3) | 57, all complete | 6.1 h |
| C | Harvetta, minima | `5b6c8d07…` | 57, all complete | 6.8 h |
| M | MATLAB `runIEM_HH` on Harvetta (Gurobi 12) | job `iem_ref_harvetta_20261005_step2b_runiem` | 57 | 3.05 h |

Notes:
- **One process each.** A, B and C each ran as a single process; no shards were needed. `logs/watch_v04.sh` started C when B finished and would have resumed any killed run; none was killed.
- **M took two attempts.** The first submission (`…_step2_runiem`, 13:40Z) stopped at 14:34Z, when a Docker status check on Tim's Mac timed out and the queue stopped the container. MATLAB itself reported no error. After inspection it was submitted once more under a new name, as the plan allows. That attempt ran from 15:28Z to 18:31Z. Both records are in `results/wbm_iem/matlab_reference/harvetta_step2/`.

## Q1: the Harvey findings on Harvetta

| Finding | Harvey (v0.3) | Harvetta (v0.4) | Reading |
|---|---|---|---|
| (a) Setup bounds identical to MATLAB | 81,094 / 81,094 | 83,521 / 83,521 (setup and global) | replicated |
| (b) Same call wherever MATLAB has values | 237 / 237 | 221 / 221 | replicated |
| (c) Biomarkers without a MATLAB optimum | 15 of 252 | 31 of 252 | replicated |
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

The published figure for Harvetta is 84.9% (Thiele et al. 2020, other model and code versions).

**Opposite-direction errors.** Twelve are shared with Harvey: citrulline in CPS1, NAGS and OTC, AADC L-DOPA and 3-methoxytyrosine, PC glucose, GMT creatinine, FED cholesterol, OXOP 5-oxoproline, LNS folate, NAGS orotate and HCYS ornithine. Harvetta adds three:
- ASNSD blood asparagine;
- MMA blood carnitine;
- BTD urinary 3-hydroxypropionate. In Harvetta its healthy maximum is 1.85 × 10⁻⁶, just above the 10⁻⁶ threshold below which the protocol sets a value to zero, so the call becomes "decreased". In Harvey both values were zero and the call was "no change". The difference between the two models is numerical noise at the threshold, not biology.

### Python against MATLAB on Harvetta (A2)

- **Bounds.** Identical after the setup (checked before the freeze) and after the global constraints.
- **Calls.** Same call on all 221 biomarkers where MATLAB returned values. All 30 differences are biomarkers that MATLAB recorded as NaN, so `runIEM_HH` scores them "unchanged"; Python calls 27 of these 30 correctly.
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
| Harvey | 217 | 207 (−10) | 217 (0) | 10 right → conflicting, 0 fixed | 1 of 18 |
| Harvetta | 218 | 205 (−13) | 218 (0) | 13 right → conflicting, 0 fixed | 3 of 16 |

**Readings (pre-declared).**
- **R1 is harmful:** it worsens on both models.
- **R2 is not useful:** it fixed no error on either model. It turned one error into another on Harvey and three on Harvetta.

**Why.**
- **The capped biomarkers.** Their minimum is zero in both states (17 of 18 in Harvey, 13 of 16 in Harvetta). The model can send the metabolite elsewhere, so a blocked pathway forces nothing into urine or blood.
- **The full-range rule's losses.** All are "conflicting": the disease maximum rises while the disease minimum falls to zero. The protocol's healthy state forces maximal flux through the IEM's reactions, which forces some downstream metabolites to be made and excreted. The disease state forces nothing.

The healthy reference is therefore an extreme state, and Shlomi-style flux-range logic does not carry over to it.

## Post hoc (not in the plan, descriptive)

**Reporting capped values as indeterminate.** A biomarker whose healthy and disease maxima are equal and positive cannot be called by a maximum-only comparison. Reporting such cases as "indeterminate" instead of "no change" would remove only errors:
- Harvey: 18 indeterminate; 217 of the remaining 233 calls correct (93.1%).
- Harvetta: 16 indeterminate; 218 of 235 correct (92.8%).

This is a reporting change, not an improvement of the model. Coverage falls from 251 to 233 and 235 biomarkers.

## Deviations and corrections

1. **MATLAB job M ran twice.** See Runs.
2. **The plan quotes Harvey's older GFR for Harvetta.** For Harvetta the older value is 128.64 ml/min (20% of renal plasma flow with the female parameters), not 129.75. The port check used the correct value; only the plan's text is affected.
3. **Analysis scripts.** Two were written after the freeze, as implementations of declared definitions:
   - `scripts/iem_error_anatomy.py` (A3) reproduces the Harvey v0.3 counts;
   - `scripts/v04_report.py` applies the frozen reading rules.

## Limits

- **Not blind.** The IEM labels and the published figures were known. Harvetta shares structure, labels, diet and constraints with Harvey.
- **One MATLAB run per model, with one solver.**
- **Direction only,** against the protocol's literature labels.

## Next

1. Independent verification of this note and of the paper 2 draft.
2. Paper 2: finish with these numbers.
3. **The cap problem needs another approach.** Options, each to be planned and frozen before testing:
   - report capped values as indeterminate;
   - define a typical healthy state rather than one with maximal IEM flux;
   - relax the binding cap of the biomarker itself.
