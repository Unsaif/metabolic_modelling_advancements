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

Note on items 2 and 6: "v_max of 545 or more" counts AGAT, whose v_max is 544.8. The 7 knockouts are those with v_max
above 500.

## 8 October 2026 (evening): the complete pin sweep, and the decision

The sweep finished at 18:57 UTC: 187 readouts × 19 development IEMs × 6 pin levels, all solves optimal. The joined
file equals a re-join of the three shards. α = 1 reproduces the matrix: 3,534 values, largest difference 9.3e-5, no
call differs. Analyses: `scripts/iem_pin_sweep_analysis.py` and `scripts/iem_pin_sweep_rules_dev.py`; outputs are in
`results/wbm_iem/pin_sweep/`. The ranking is among the 19 development knockouts (chance MRR 0.187).

| Rule | Share called increased (median) | Own lab correct | Lab MRR, plain | Lab MRR, adjusted | HPO MRR, plain | HPO MRR, adjusted |
|---|---|---|---|---|---|---|
| α = 1 (current protocol) | 0.30 | 65 of 84 | 0.442 | 0.623 | 0.404 | 0.547 |
| α = 0.5 | 0.06 | 30 of 84 | 0.619 | 0.659 | 0.604 | 0.684 |
| α = 0.1 | 0.01 | 23 of 84 | 0.509 | 0.577 | 0.545 | 0.592 |
| α = 0.001 (increases only where coupled) | 0.01 | 23 of 84 | 0.522 | 0.576 | 0.561 | 0.610 |
| α = 0 | 0.00 | 8 of 84 | 0.286 | 0.363 | 0.420 | 0.554 |
| α = 1 with capped increases indeterminate | 0.10 | 33 of 84 | 0.639 | 0.587 | 0.538 | 0.565 |
| Pin capped at 100 mmol/day | 0.08 | 35 of 84 | 0.603 | 0.622 | 0.614 | 0.589 |

**Paired comparison with α = 1** (19 lab profiles; 16 HPO profiles; bootstrap 95% intervals):

- **α = 0.5, plain score:** lab +0.18 [−0.02, +0.37], 11 profiles better and 6 worse; HPO +0.20 [0.00, +0.40].
- **α = 0.5, adjusted score:** lab +0.04 [−0.15, +0.23], 5 better and 5 worse; HPO +0.14 [−0.05, +0.33].
- **Capped increases indeterminate, plain:** lab +0.20 [+0.03, +0.37]; adjusted: lab −0.04.

**Where the change comes from.**

- At α = 0.5, the low-v_max diseases move to the top: CYP21D, DPYR, HIS, HLYS1, HYCARO, IVA and TYR3 go from
  0.1–0.5 to 1.0, because the large-v_max knockouts stop matching every profile.
- The large-v_max diseases fall: GA2 from 0.75 to 0.07, SSADHD 0.52 to 0.09, HMG 0.46 to 0.19, HYPRO1 0.37 to 0.16.
  Their own matches were drained, capped readouts.

**Decision: stop without a held-out test.** The rule agreed with Tim on 8 October was to stop unless the development
data showed a clear win. They do not:

- Every interval for the plain score includes or touches zero.
- Against the adjusted score, which is the current best method, the lab gain is +0.04, with 5 profiles better and 5
  worse.
- The gain comes from the low-v_max diseases, at the expense of the large-v_max ones. The held-out set has 24 of 38
  IEMs with v_max above 500 (20 above 10,000), against 7 of 19 in development, so the gain would probably not carry
  over.
- Own-biomarker agreement falls from 65 to 30 of 84.

**What the sweep shows.** The forced flux sets a trade-off between sensitivity and specificity. At the full pin, most
apparent hits on known increases are capped readouts drained to zero, which any large forced flux produces. Real,
coupled increases are fewer and survive any pin. Fixing this needs model work on caps, loops and disease definitions,
not another protocol setting. That work is the next direction (see `docs/briefs/` and the ideas document).
