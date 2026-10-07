# IEM disease ranking: results

7 October 2026.

- **Plan.** `wbm-iem-ranking-plan.md` (committed 347b21f, amended before the run in 34ab349, a6e98ef and def01b0).
- **Matrix.** `results/wbm_iem/Harvey_1_03d_iem_cross_ranking_v1.json` (sha256 11779e84…), joined from three shards.
- **Ranking.** `results/wbm_iem/ranking/Harvey_1_03d_ranking_v1.json`.
- **Checks.** `results/wbm_iem/ranking/Harvey_1_03d_ranking_v1_checks.json`, the cross-solver recheck
  `Harvey_1_03d_ranking_v1_recheck_highs.json`, and the independent check `results/wbm_iem/ranking/independent_check_results/`.

One deviation after the run started: the planned rechecks did not re-solve anything and were replaced (see
`wbm-iem-ranking-deviations.md`). Analyses marked *post hoc* were not in the plan.

## Summary

- **Primary result.** For each of the 57 simulated inborn errors of metabolism (IEMs), we ranked all 57 knockouts by
  how well their predicted biomarker changes match that disease's known biomarkers. The model ranks the right disease
  well above chance:
  - mean reciprocal rank (MRR) 0.25 against 0.08 by chance; no random draw out of 100,000 reached it;
  - first in 8.8 of the 57 profiles in expectation (1.0 by chance);
  - in the top five in 17.4 (5.0 by chance).
- **Why the effect is modest in absolute terms (post hoc, descriptive).** The protocol predicts "increased" for most
  readouts in most knockouts, a median of 80% of the 186 readouts. An unrelated knockout predicts the direction of a
  disease's known biomarker 50% of the time; the true disease does so 86% of the time.
  - Diseases whose profiles contain a known decrease are often ranked uniquely first (6 of 19).
  - Diseases whose profiles are all increases never are (0 of 38); the true disease ties with a median of 27.5
    others.
- **Pre-specified adjustment.** Discounting each knockout's overall tendency to raise readouts raises the MRR by about
  two-thirds, to 0.42: first in 18.5 profiles and in the top five in 29, with a median rank of 5. The gain is uneven:
  diseases whose own knockout changes few readouts gain, and those whose knockout raises almost everything lose.
- **Post hoc: accuracy alone is uninformative.** Because 85% of the known biomarkers are increases, a rule that
  calls every biomarker "increased" would score 85.3% accuracy. The protocol scores 86.5% on Harvey (v0.3), and the
  published figures are about 85%. The protocol does carry disease-specific information, as the ranking shows, and it
  gets 25 of the 37 known decreases right. Agreement percentages cannot show either.

## The prediction matrix

**Size.** 57 IEMs × 187 readouts, each maximised in the healthy and disease state of every IEM: 21,204 solves. On top
of that are the protocol's 171 set-up solves, repeated in each shard.

**Computation.** Gurobi 13.0.1 ran on Tim's MacBook in 3 parallel shards of 3 threads, in readout order with
concurrent warm starts. It took 19.4 hours, from 19:37 on 6 October to 15:00 on 7 October (Irish time), with a median of 9.3 s per solve.

The checks before ranking:

| Check | Result |
|---|---|
| Completeness | All 57 IEMs complete, 187 readouts each. Every solve ended optimal, and no barrier fallback was needed. The only unavailable readout is `EX_25aics[u]`, which is not in Harvey. |
| Own biomarkers against v0.3 (HiGHS) | **252 of 252 calls the same** (stop rule: more than 5 differing). Largest absolute value difference 0.002; largest relative difference 6e-5 among values above 1e-3. |
| Rechecks | **Not performed as planned.** The 197 hash-selected rechecks returned the stored solution without re-solving (0.000 s each), because the model was not reset; this was found by the independent check. They were replaced by re-solving the same 197 LPs from scratch with HiGHS (interior point with crossover), a different solver: RECHECK_RESULT. |
| HIS against the one-IEM-at-a-time feasibility run | 161 of 161 calls the same; largest difference 1.5e-6. |
| Re-join of the three shards in the cloud | Identical to the joined file written on Tim's Mac. |

## Primary analysis: lab profiles, protocol calls, plain score

| | Observed | Chance |
|---|---|---|
| Mean reciprocal rank | **0.252** | 0.081 |
| Expected number ranked first | 8.8 of 57 | 1.0 |
| Expected number in the top five | 17.4 of 57 | 5.0 |
| Median expected rank | 15.5 | 29 |
| Monte Carlo p (MRR, top 1, top 5) | < 0.0001 each | |

No random draw reached the observed values, so p is (1 + 0)/(1 + 100,000).

By profile size:

| Biomarkers in the profile | Profiles | MRR | Expected first | Expected top five | Median rank |
|---|---|---|---|---|---|
| 1 | 10 | 0.18 | 0.8 | 2.3 | 18.0 |
| 2–3 | 16 | 0.19 | 1.5 | 3.6 | 14.8 |
| 4–6 | 21 | 0.35 | 5.3 | 8.9 | 14.5 |
| ≥ 7 | 10 | 0.22 | 1.3 | 2.6 | 14.8 |

**Score ties.**

- In 47 profiles no knockout scores higher than the true disease, but most profiles end in large ties.
- The true disease is uniquely first in 6 profiles: AGAT, STAR, GMT, GACR, LTC4S and ASNSD. All six contain a known
  decrease.
- **Ranked in the bottom half:** MMA and HCYS (rank 41), FED (33.5), HLYS2 (32.5) and DPYR (31.5).
  - For MMA, HCYS and DPYR, the true disease's own calls miss most of its biomarkers: 2 of 7, 2 of 7 and 1 of 4.
  - FED and HLYS2 miss half.

## Secondary and sensitivity analyses (pre-specified)

| Profiles | Calls | Score | Profiles | MRR | Expected first | Expected top five | Median rank |
|---|---|---|---|---|---|---|---|
| lab | protocol | plain | 57 | **0.252** | 8.8 | 17.4 | 15.5 |
| hpo | protocol | plain | 39 | 0.204 | 4.5 | 9.7 | 18.0 |
| hpo_frequent | protocol | plain | 38 | 0.231 | 5.5 | 10.5 | 17.8 |
| lab_hpo_corroborated | protocol | plain | 29 | 0.248 | 4.6 | 8.2 | 17.0 |
| lab | protocol | adjusted | 57 | 0.423 | 18.5 | 29.0 | 5.0 |
| hpo | protocol | adjusted | 39 | 0.351 | 10.1 | 16.4 | 16.5 |
| hpo_frequent | protocol | adjusted | 38 | 0.382 | 11.0 | 17.1 | 15.5 |
| lab_hpo_corroborated | protocol | adjusted | 29 | 0.421 | 9.5 | 15.0 | 5.0 |
| lab | material | plain | 57 | 0.272 | 9.9 | 18.4 | 15.0 |
| hpo | material | plain | 39 | 0.272 | 7.0 | 12.7 | 15.5 |
| lab | material | adjusted | 57 | 0.462 | 20.0 | 33.1 | 3.5 |
| hpo | material | adjusted | 39 | 0.409 | 12.6 | 18.1 | 9.0 |

- Every analysis is far above chance: the chance MRR is 0.081, and Monte Carlo p < 0.0001 throughout. The full table,
  with hpo_frequent and the corroborated set under every rule, is in the ranking file.
- **HPO profiles.** They come from an independent source, Orphanet's HPO annotations, and give the same picture as the
  lab's own biomarker list, a little weaker.
- **The adjusted score** discounts each knockout's base rate of increases and decreases.
  - It helps profiles made only of increases when the true knockout is specific: 8 of the 38 such lab profiles become
    uniquely first, against none with the plain score. Their knockouts call only 1% to 19% of readouts increased
    (HIS, 3MCC, HYCARO, FIGLU, GA1, MSUD, PKU, TYR3).
  - The gain is uneven. The adjusted rank follows the true knockout's own increase share (Spearman 0.70; 0.06 for the
    plain score). 27 profiles rank better and 22 worse, for example DESMO 21.5 → 49, GA2 6 → 26.5 and BTD 26 → 38.
- **The material-change rule** (a change above 1e-3 and above 5%) helps a little. It removes the smallest changes but
  not the generic increases.

## Confusions (descriptive)

- **Plain score.** The median knockout ties with a positive-scoring true disease in 31 profiles (range 0 to 40). No
  knockout outranks the true disease in more than 7 profiles. The confusion is generic rather than between related
  diseases.
- **Adjusted score.** The knockouts that most often outrank the true disease are:
  - AADC (22 profiles);
  - BTD, SUCLA and TETB (21 each);
  - AGAT and GMT (20);
  - LTC4S (19);
  - PC (18);
  - the urea-cycle disorders ARG, CPS1 and OTC (18 each, also tying in 2) and CIT1 (16);
  - FED (16).

  Only AADC, TETB, LTC4S and FED change fewer than half the readouts (55 to 87 of 186). The others change 154 to 168,
  about the median of 154. They outrank true diseases whose own knockouts raise even more readouts, which the
  adjustment discounts more. Some confusions are biologically expected: within the urea cycle, ammonia, glutamine and
  orotate are shared readouts.

## Why the plain ranking is limited (post hoc, descriptive)

**Generic increases.**

- Under the protocol's rule, the median knockout calls 80% of the 186 readouts increased (range 1% to 95%) and 1%
  decreased.
- A readout is called increased by a median of 32 of the 57 knockouts. 31 readouts are increased by at least 40.
- With the material-change rule, the median falls to 68%.

**Matching a disease's known biomarkers.**

- Across the 251 scored lab tuples, the true disease predicts the known direction 86.5% of the time; any other
  knockout does so 50.4% of the time (7,080 of 14,056).
- Known decreases are where the predictions are specific, because predicted decreases are rare.

**Where the increases come from.**

- **What follows from the protocol.** In the healthy state, the disease pathway is forced to run at its maximum flux;
  in the disease state it is off. Without the forced flux, the disease state's feasible set would lie inside the
  healthy state's: zero lies within the bounds of all 1,541 IEM reactions after the Toolbox set-up, as the
  independent check confirmed. So no maximum could rise, and every predicted increase comes from the forced healthy
  state, specific or not.
- **What is not tested.** Why forcing one pathway to its maximum lowers the maxima of so many unrelated readouts,
  presumably through shared resources such as diet uptake, organ coupling and energy.
- **How to test it.** Relaxing the pin completely is not a test, because the outcome is fixed. A useful test would
  force only part of the maximum flux, or compare against a typical healthy flux state.

**Accuracy against a trivial baseline (post hoc).**

| | Harvey v0.3 | Harvetta v0.4 |
|---|---|---|
| Protocol, correct of 251 scored | 217 (86.5%) | 218 (86.9%) |
| "Always increased", correct | 214 (85.3%) | 214 (85.3%) |
| Protocol on known increases | 192 / 214 | 194 / 214 |
| Protocol on known decreases | 25 / 37 | 24 / 37 |

The published 84.9% and 85.3% are likewise at the base rate. Agreement percentages should be reported with this
baseline. Specificity, as measured here, should be reported with them.

## Limits

- **Not a held-out or clinical test.** The lab profiles are the biomarker lists the protocol was built around, and
  have been used in our development since v0.2. Profiles are literature biomarker lists, not patient measurements.
- **HPO profiles.** They are independent of the lab's list, but are mapped by our own rules. 21 of the 39 profiles have
  only one or two readouts, and hyperammonaemia appears in 15 of the 39.
- **One model, one context.** Harvey with the protocol context. The minimal-context and Harvetta replications were
  optional and were not run.
  - The context matters for the forced healthy state, because each IEM's maximum flux is computed with its own sinks
    open.
  - In the independent check's runs, LTC4S's maximum flux fell from 125.5 to 0.0087 without its sinks. Its own calls,
    and those of STAR and CYP21D, did not change.
  - Whether the ranking holds without the sinks is untested.
- **The 1e-6 call threshold.** The protocol's threshold makes the smallest differences count. The material-change
  sensitivity addresses this only partly.

## Per-profile ranks (lab profiles)

The table gives the expected rank of the true disease among 57. "Own" is the number of the profile's biomarkers whose
known direction the true disease predicts.

| IEM | Biomarkers (up/down) | Own | Plain | Adjusted |
|---|---|---|---|---|
| 2OAA | 3 (3/0) | 3 | 9.0 | 2.0 |
| 3MCC | 3 (3/0) | 3 | 12.0 | 1.0 |
| 3MGA | 1 (1/0) | 1 | 13.0 | 2.0 |
| AADC | 10 (8/2) | 8 | 16.0 | 1.0 |
| ADSL | 1 (1/0) | 1 | 18.0 | 30.5 |
| AGAT | 2 (0/2) | 2 | 1.0 | 1.0 |
| AKGD | 5 (5/0) | 5 | 19.5 | 34.0 |
| AMA1 | 2 (2/0) | 2 | 21.5 | 40.0 |
| ARG | 6 (6/0) | 6 | 13.0 | 3.0 |
| ASA | 7 (7/0) | 7 | 14.0 | 16.5 |
| ASNSD | 4 (2/2) | 4 | 1.0 | 1.0 |
| BTD | 8 (8/0) | 6 | 26.0 | 38.0 |
| CD | 3 (3/0) | 3 | 15.5 | 14.5 |
| CIT1 | 6 (6/0) | 6 | 15.5 | 15.5 |
| CMO1 | 1 (0/1) | 1 | 1.5 | 1.0 |
| CPS1 | 6 (5/1) | 5 | 14.5 | 17.5 |
| CYP21D | 6 (5/1) | 4 | 16.0 | 1.0 |
| DESMO | 1 (1/0) | 1 | 21.5 | 49.0 |
| DGK | 1 (1/0) | 1 | 20.5 | 29.5 |
| DPYR | 4 (4/0) | 1 | 31.5 | 44.0 |
| EF | 1 (1/0) | 1 | 21.0 | 36.0 |
| EP | 1 (1/0) | 1 | 25.5 | 39.0 |
| FED | 2 (1/1) | 1 | 33.5 | 33.5 |
| FIGLU | 1 (1/0) | 1 | 16.0 | 1.0 |
| GA1 | 2 (2/0) | 2 | 7.0 | 1.0 |
| GA2 | 15 (15/0) | 15 | 6.0 | 26.5 |
| GACR | 8 (4/4) | 6 | 1.0 | 1.0 |
| GMT | 5 (3/2) | 4 | 1.0 | 1.0 |
| HCYS | 7 (6/1) | 2 | 41.0 | 6.0 |
| HIS | 5 (5/0) | 5 | 16.5 | 1.0 |
| HLYS1 | 3 (3/0) | 3 | 17.5 | 1.5 |
| HLYS2 | 4 (4/0) | 2 | 32.5 | 2.5 |
| HMET | 4 (4/0) | 4 | 2.0 | 16.0 |
| HMG | 2 (2/0) | 2 | 13.0 | 20.5 |
| HPC | 4 (4/0) | 3 | 20.0 | 2.0 |
| HPII | 6 (6/0) | 6 | 16.0 | 17.5 |
| HYCARO | 1 (1/0) | 1 | 15.5 | 1.0 |
| HYPRO1 | 4 (4/0) | 4 | 16.0 | 23.0 |
| HYPVLI | 3 (3/0) | 3 | 16.5 | 3.0 |
| IVA | 3 (3/0) | 2 | 20.0 | 3.0 |
| LNS | 3 (2/1) | 2 | 20.0 | 34.0 |
| LTC4S | 6 (0/6) | 6 | 1.0 | 1.0 |
| MMA | 7 (6/1) | 2 | 41.0 | 3.0 |
| MSUD | 6 (6/0) | 5 | 4.0 | 1.0 |
| NAGS | 6 (4/2) | 4 | 24.0 | 41.0 |
| OTC | 9 (8/1) | 8 | 15.5 | 15.0 |
| OXOP | 8 (7/1) | 7 | 11.5 | 25.0 |
| PC | 15 (14/1) | 14 | 6.5 | 11.0 |
| PHOX1 | 3 (3/0) | 3 | 12.0 | 27.0 |
| PKU | 3 (3/0) | 3 | 9.0 | 1.0 |
| SSADHD | 5 (5/0) | 5 | 6.0 | 10.0 |
| STAR | 5 (0/5) | 5 | 1.0 | 1.0 |
| SUCLA | 2 (2/0) | 2 | 20.5 | 17.0 |
| TETB | 1 (1/0) | 1 | 18.0 | 5.0 |
| TYR1 | 4 (4/0) | 4 | 14.0 | 2.0 |
| TYR3 | 3 (3/0) | 3 | 14.0 | 1.0 |
| XAN1 | 5 (3/2) | 5 | 1.5 | 1.0 |
