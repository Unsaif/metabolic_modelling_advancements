# Independent check: IEM disease-ranking results (Harvey 1.03d)

7 October 2026. Written by an independent checking agent after the ranking had been computed. It returned the report as
a message, and the main session saved it here verbatim apart from formatting. Its scripts are in this directory.

## Summary

I found **0 critical, 2 major and 6 minor problems.**

- **Numbers reproduce.** Every ranking metric, per profile and in summary, matches exactly when recomputed with my own
  code. The joined matrix is sound.
- **The rechecks never re-solved anything.** Plan check 3 was therefore not performed, although the write-up reports it
  as passed.
- **One explanation is contradicted by the data.** The account of the adjusted-score confusers is wrong for most of the
  knockouts it names.

The repository is unchanged.

## Problems

### 1. Major: the 197 rechecks were no-ops

Check 3 was therefore not performed, and "There were no deviations" is wrong.

- **Cause.** `scripts/run_wbm_iem_cross.py` (lines 431–434) re-solves right after the warm solve with nothing changed.
  `GurobiWBM.solve` only sets `Params.Method` and calls `optimize()`, with no `model.reset()`, so Gurobi returns the
  solution it already has.
- **Evidence in the matrix.** 196 of the 197 recheck times are 0.000 s and one is 0.001 s; a barrier solve here has a
  median of 7.9 s. All 197 `abs_diff` values are exactly 0.
- **Reproduced.** In gurobipy 13.0.1 on a test LP, setting Method=2 and calling `optimize()` after a concurrent solve
  gave IterCount 0, BarIterCount 0 and 1.5e-5 s. With `reset()` it re-solved.
- **What independent checking remains.** Only about 814 of the 21,204 values (3.8%, not a random sample) are
  reproduced independently: own biomarkers against HiGHS and HIS against the feasibility run.
- **Fix.**
  - Add `reset()` before the recheck solve.
  - Re-solve the 197 hash-selected LPs from scratch (about 30 minutes).
  - Record the deviation in `wbm-iem-ranking-deviations.md`.
  - Reword the check table.

### 2. Major: the explanation of the adjusted-score confusers is wrong

- "These knockouts change few readouts" holds only for AADC (80 of 186 readouts changed), TETB (87) and LTC4S (74).
- It fails for BTD (166), SUCLA (160), PC (168), CIT1 (166), and AGAT, GMT, ARG, CPS1 and OTC (154 each). The median
  knockout changes 154.
- These knockouts outrank true diseases whose own increase share is even higher.
- FED (16) is left out of the list, although CIT1 (16) is in it.

### 3. Minor: another explanation is partly wrong

Under "Ranked below most knockouts … own calls miss most of its biomarkers":

- BTD has 6 of 8 own matches and an expected rank of 26.0, so it is not below most knockouts.
- FED (1 of 2) and HLYS2 (2 of 4) miss only half.
- The statement fits only MMA, HCYS and DPYR.

### 4. Minor: number slips

- "Ties with a median of 27 others" should be 27.5 (38 profiles; middle values 26 and 29).
- "87%" should be 86% (217/251 = 86.45%).
- "Almost doubles" overstates the change from 0.252 to 0.423, which is ×1.68.
- "One or two mapped terms (21 of 39)": 21 counts readouts; by distinct HPO terms it is 20.
- "Across the 252 lab tuples": 251 are scored.
- "Largest relative difference 6e-5" holds only for values above 1e-3, which is not stated. Below that, for example,
  HMET `EX_hcys_L[u]` healthy is 2.9e-6 with Gurobi and 9.7e-5 with HiGHS.

### 5. Minor: the adjustment's gain is uneven, and the write-up does not say so

- The adjusted rank tracks the true knockout's increase share (Spearman 0.70; 0.06 for the plain score).
- The 8 all-increase profiles that become uniquely first have increase shares of 1–19%, against a median of 80%.
- 22 of the 57 profiles rank worse under the adjustment, for example DESMO 21.5→49, GA2 6→26.5 and BTD 26→38.

### 6. Minor: the mechanism is stated more weakly than it can be

The logic is correct, and the section is labelled untested.

- I checked that 0 lies within the bounds of all 1,541 IEM reactions after the Toolbox set-up. Bound tweaks and global
  constraints only set bounds to 0 (or upper bounds to 100).
- So without the pin, the disease state's feasible set lies inside the healthy state's, and no increase is possible
  beyond tolerance. That part is provable, not "likely".
- What remains untested is why the pin lowers unrelated readouts' maxima.
- A fully relaxed pin is not a test, because the outcome is already fixed. A useful test would use a partial pin (a
  fraction of v_max) or a typical healthy flux state.

### 7. Minor: the dependence on context is not disclosed

Each IEM's v_max is computed with its own sinks open.

- I ran the minimal context with HiGHS for STAR, LTC4S and CYP21D. All 17 own calls were unchanged.
- But LTC4S's v_max fell from 125.5 to 0.0087, so the strength of the forced healthy state depends on each IEM's sinks.
- Whether the ranking holds without them is untested, because the replication was not run. The "One model, one
  context" limit should say this.

### 8. Minor: post hoc numbers are not labelled

The Summary's "Why the effect is modest" bullet (80%, 50% vs 87%, 6 of 19, 0 of 38) is post hoc but not marked as
such.

## Verified correct (my own recomputation)

### Matrix

- The joined file equals my independent interleaving of the three shards (0 differing entries). My re-join is
  byte-identical (sha256 11779e84…).
- 57 IEMs × 187 readouts in the same order, all complete. 10,602 solve pairs ended optimal; only `EX_25aics[u]` is
  absent; no fallbacks.
- All 10,659 calls follow the rule: 5,891 Increased, 783 Decreased, 3,928 Unchanged, 57 NA.
- The 197 recheck positions match the hash selection.

### Own biomarkers

- 252 of 252 calls are the same as v0.3. The largest absolute difference is 0.00195 (XAN1 `DM_urate[bc]`); the largest
  relative difference (values above 1e-3) is 6.1e-5; v_max agrees within 6e-7.
- HIS against the feasibility run: 161 of 161 calls the same, largest difference 1.53e-6.

### Ranking

- All 16 analyses are identical to the ranking file: per profile, the number scoring higher, ties, expected rank,
  expected reciprocal rank, P(top 1) and P(top 5), score, candidate lists, own matches and confuser counts.
- Primary: MRR 0.2521, expected first 8.83, top five 17.37, median rank 15.5, 47 profiles with no candidate higher.
  Exact null MRR H₅₇/57 = 0.08121.
- Adjusted (lab): MRR 0.4232, 18.5 first, 29.0 in the top five, median rank 5.0.
- My Monte Carlo (100,000 draws, my own seed): no draw reached the observed MRR, top-1 or top-5 count in any of the 48
  tests. The largest null MRR for the primary analysis was 0.146.
- The strata table, the secondary table and all 57 rows of the per-profile table match.

### Post hoc claims

- Median 80.1% of readouts called increased (range 1.1–94.6%), 1.1% decreased.
- A readout is increased by a median of 32 knockouts; 31 readouts are increased by 40 or more. Under the material rule
  the median share is 68.3%.
- True disease 86.45% against other knockouts 50.37% (7,080 of 14,056).
- Uniquely first: 6 of 19 and 0 of 38 (plain), 8 of 38 all-increase profiles under the adjustment.
- Plain confusions: median 31 ties (range 0–40); no knockout outranks the true disease in more than 7 profiles.
- Baseline: protocol 217 and 218 of 251, "always increased" 214; increases 192/214 and 194/214, decreases 25/37 and
  24/37. The comparison is fair: both use the same 251 scored tuples, and the one excluded tuple is an increase.

### Timing and HPO limits

- 19:37:29 plus 69,782 s gives 15:00:31 (19.4 h); median 9.30 s per solve.
- 21 of 39 HPO profiles have 1–2 readouts; `DM_nh4[bc]` (hyperammonaemia) is in 15 of 39.

### Plan

- All amendments precede the run: def01b0 at 19:30:45 Irish time, run start 19:37:29.
- The strata, confusions and secondary analyses the plan requires are present.
