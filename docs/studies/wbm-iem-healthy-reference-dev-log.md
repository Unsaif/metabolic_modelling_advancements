# Healthy reference state: development log

Every development look, in order, so that the confirmatory plan can say what was tried. Plan:
`wbm-iem-healthy-reference-dev-plan.md`. Only development IEMs are scored here: 19 profiles, with the 19 development
knockouts as candidates, so the chance mean reciprocal rank (MRR) is 0.187.

## 7 October 2026

**1. Scoring variants on the current matrix** (no new LPs; development profiles only).

| Score | Lab, protocol calls | HPO, protocol calls | Lab, material rule | HPO, material rule |
|---|---|---|---|---|
| Plain | 0.442 | 0.404 | 0.456 | 0.463 |
| Adjusted (pre-specified in the ranking plan) | 0.623 | 0.547 | 0.715 | 0.592 |
| Information content of each predicted change (as in phenotype-matching tools) | 0.440 | 0.368 | 0.446 | 0.459 |
| Likelihood ratio, with own-call rates from development profiles | 0.539 | 0.540 | 0.573 | 0.573 |

Neither new score beats the adjusted score, so no further scoring variants will be tried. The adjusted score's
held-out results are already known from the ranking study, so it cannot be tested on the held-out set again.

**2. Development knockouts under the current protocol.**

- The development IEMs fall into two groups:
  - 12 IEMs have v_max of 79 or less, and call 1% to 34% of readouts increased.
  - 7 IEMs have v_max of 545 or more, and call 79% to 93% increased.
- The development set's median v_max (67) is lower than that of all 57 IEMs (793). It therefore holds fewer of the
  knockouts that raise almost everything. The split was stratified by profile, not by v_max.

**3. Scale of v_max** (development IEMs; EU average diet as set up for the protocol).

- HYPRO1's v_max (92,160) is about 1,400 times its dietary proline supply (43 to 65 mmol/day).
- HIS's v_max (66.7) is 4 to 6 times its histidine supply (11 to 17 mmol/day).
- Whole-body O2 uptake is bounded to 15,000–25,000 mmol/day.
- The flux distribution at v_max, from interior point with crossover, carries loop fluxes of ±800,000 at their
  bounds. It cannot be read without removing loops first.

## 8 October 2026 (night)

**4. Reaction sets** (set-up information only; no prediction was looked at). The name patterns, copied from
`runIEM_HH`, leave part of two deficiencies active:

- FED's list skips `_LCAT6e` to `_LCAT9e`: 7 reactions (liver and adrenal gland).
- HYPRO1's patterns (`_r1453`, `_PROD2m`, `_PRO1xm`) leave 9 cytosolic `PROD2` reactions active in the disease state.

FED is a held-out IEM; only its reaction list was read. Any fix to these sets would be a separate, declared change.

**5. Pin sweep timing.** The sweep started at 23:11 UTC on 7 October, as 3 shards of 3 threads. The solves after a
change of pin level do not warm-start faster than barrier, at about 8 s each. One readout takes about 950 s per shard,
so the 187 readouts should take about 16 to 17 hours.
