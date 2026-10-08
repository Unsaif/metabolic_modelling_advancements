# Healthy reference state: development plan

7 October 2026. Follows the disease-ranking study (`wbm-iem-ranking-results.md`). This plan and the
development/held-out split were committed before any variant of the protocol was computed. A separate confirmatory
plan will be checked independently and committed before the held-out test.

## Question

The protocol forces the healthy state to carry the IEM reactions' maximum summed flux, v_max. Every predicted
increase comes from that forcing, and the median knockout calls 80% of readouts increased.

Can a less extreme healthy reference make the predictions more disease-specific, as measured by the disease ranking,
without losing the known biomarkers?

## What is known before this plan

- **The ranking study, for all 57 IEMs, under the current protocol:** per-profile ranks, own-biomarker calls and each
  knockout's share of increased calls (`wbm-iem-ranking-results.md`).
- **v0.4:** minima were computed in the protocol's healthy and disease states. The full-range rule lost correct calls.
  Capped biomarkers had zero minima in both states in 17 of 18 cases on Harvey.
- **Looked at while writing this plan:** from the existing matrix, in aggregate over all 57 knockouts; no profile was
  scored.
  - v_max ranges from 0.004 to 92,160 (median 793), and 23 IEMs have v_max above 10,000.
  - Across knockouts, v_max and the share of readouts called increased have a Spearman correlation of 0.75.
  - The healthy maximum is zero while the disease maximum is positive for a median of 7.5% of readouts.
- **The Toolbox already has the parameter.** `checkIEM_WBM` takes the forced fraction as `minRxnsFluxHealthy`, and
  `runIEM_HH` sets it to 1.

## Split

`data/iem/iem_ranking_split_v1.json`, made by `scripts/make_iem_ranking_split.py` (seed 20261007). It has 19
development and 38 held-out IEMs, stratified by whether the lab profile contains a decrease and by profile size.

Development IEMs: 3MGA, AGAT, CMO1, CPS1, CYP21D, DPYR, EP, GA2, GACR, HCYS, HIS, HLYS1, HMG, HYCARO, HYPRO1, IVA,
MMA, SSADHD, TYR3.

Rules until the confirmatory plan is committed:

- New computations for a variant are made for development IEMs only.
- Development outcomes (calls, own-biomarker matches, ranks) are computed for development profiles only, with the 19
  development knockouts as the candidates.
- No prediction of a variant for a held-out IEM is computed or looked at.
- The held-out IEMs' results under the current protocol are already known. They are not used to choose a variant.

## Candidates

1. **Partial pin.** This is the Toolbox's `minRxnsFluxHealthy`: the healthy state needs summed IEM flux of at least
   α·v_max, truncated to six decimals. α takes the values 0.5, 0.1, 0.01 and 0.001. The disease state is unchanged.
2. **Scoring changes on the current matrix,** with no new LPs: for example, weighting each readout by how specific its
   change is across knockouts. These are developed on development profiles.
3. **If neither helps,** a rule without a pin: increases from what the block forces (minimum total disposal), and
   decreases from what it prevents (maximum). Another option is a pin scaled to the pathway's substrate supply. Either
   would be planned separately.

## First computation

`scripts/run_wbm_iem_pin_sweep.py` computes, for each development IEM and each of the 187 readouts, the healthy maximum
at α = 1, 0.5, 0.1, 0.01, 0.001 and 0, in that order.

- **The α = 1 values are a check.** They should reproduce the matrix.
- **α = 0 is a reference.** The pin is then only sum ≥ 0, which the disease state meets, so no readout can be higher
  in the disease state. These values show what the block removes: the decrease side of a rule without a pin.
- **Only the pin's bound changes between the levels,** so the solves after the first warm-start.
- **The disease maxima come from the matrix,** because the disease state does not depend on α.
- **Everything else is as in the matrix:** the model, setup, protocol context, sinks and tolerances. The IEM set-ups
  are rebuilt from the matrix's records (v_max, pin and reactions), as the rechecks did. Gurobi 13.0.1 runs on Tim's
  Mac, as for the matrix.

## What is looked at in development

For each α, on development IEMs only:

- each knockout's share of readouts called increased and decreased;
- own-biomarker matches, for the lab and HPO profiles;
- the ranking among the 19 development knockouts: mean reciprocal rank, plain and adjusted, for lab and HPO
  profiles;
- how each own increase and each other increase changes with α.

## After development

A confirmatory plan, checked independently and committed before its run, will fix:

- the variant or variants to test;
- the full computation, for all 57 IEMs;
- the primary outcome: mean reciprocal rank on the 38 held-out lab profiles, the variant against the current
  protocol, with a paired test;
- the secondary outcomes: HPO profiles and own-biomarker balanced accuracy;
- the decision rules.

## Amendments before the first development run

Committed with the sweep script, before the development run.

- **α = 0 added** to the sweep (see above).
- **Smoke test of the launcher.** HiGHS in the cloud, HIS only, two readouts, α = 1 and 0.5, two shards.
  - α = 1 reproduced the matrix: urinary histamine 0.000285 and blood histamine 50.0 in the healthy state.
  - At α = 0.5, urinary histamine's healthy maximum rose to the disease value, 29.0, so that own increase would become
    "unchanged". Blood histamine's rose from 50.0 to 83.4, still below the disease value of 116.7.
  - HIS is a development IEM, so this is within the rules.

## Outcome

8 October 2026: the development sweep finished, and the study stopped without a held-out test. The development data
showed no clear win over the current protocol with the adjusted score. The held-out set is dominated by the large-v_max
diseases that a gentler pin hurts. The results and reasoning are in `wbm-iem-healthy-reference-dev-log.md`. The
held-out IEMs remain unused by any variant.
