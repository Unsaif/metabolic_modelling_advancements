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

## 8 October 2026 (morning): partial look at the pin sweep

At 06:31 UTC, 73 of the 187 readouts were finished, and all 8,322 solves had ended optimal. This look is partial and
descriptive, on development IEMs only. No variant is chosen from it.

- **Check.** α = 1 reproduces the matrix: 1,387 values, largest difference 2.8e-5, no call differs.
- **Share of finished readouts called increased** (median over the 19 development knockouts): 0.26 at α = 1, 0.08 at
  0.5, 0.03 at 0.1, 0.01 at 0.01 and 0.001, and 0 at 0. For the 7 knockouts with v_max of 545 or more, it falls from
  0.95 to 0.10 at α = 0.5 and to 0 to 0.04 at α = 0.01.
- **Own lab biomarkers with a finished readout** (50 tuples):

  | α | Correct | Increases | Decreases |
  |---|---|---|---|
  | 1 | 43 | 39 of 43 | 4 of 7 |
  | 0.5 | 17 | 11 of 43 | 6 of 7 |
  | 0.1, 0.01, 0.001 | 15 | 9 of 43 | 6 of 7 |
  | 0 | 6 | 0 of 43 | 6 of 7 |

- **Two kinds of increase at α = 1:**
  - *Capped readouts.* The healthy maximum equals the disease maximum at every α below 1, and falls to about zero
    only at α = 1. Examples: HIS urinary histidine (29.0), HYPRO1 and GA2 urinary proline (889.9), and GA2 blood C4
    and C10 carnitines (50). These increases come from draining the healthy state, which any knockout with a large
    forced flux also does.
  - *Coupled readouts.* The healthy maximum falls in proportion to the pin. Examples: HIS blood histidine
    (66.7 − 66.7α) and CPS1 blood glutamine (793 − 793α). These stay increased at every α above 0, and look specific.
- **Implication, not yet tested.** A single gentler pin removes the generic increases, but also the capped own
  increases. Once the sweep is complete, two rules are worth comparing on the development set:
  - calling increases from coupling at a small pin;
  - reporting capped readouts as indeterminate.
