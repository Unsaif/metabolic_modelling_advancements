# IEM disease ranking: deviations from the plan

Plan: `wbm-iem-ranking-plan.md`. Amendments made before the run are in the plan itself and its git history
(347b21f, 34ab349, a6e98ef, def01b0). Listed here are the changes made after the run started (6 October 2026,
19:37 Irish time).

| # | Deviation | When decided | Why | Effect |
|---|---|---|---|---|
| 1 | The 197 planned rechecks (barrier from scratch, about 1 in 100 warm solves) were not performed. The same 197 LPs were re-solved from scratch with HiGHS 1.15.1 (interior point with crossover) in the cloud, by `scripts/recheck_cross_matrix.py`. | 7 October, after the ranking was computed. Found by the independent check of the results. | The runner called the barrier solve without resetting the model. Gurobi then returned the solution it already had (0.000 s each, difference exactly 0). The fault is fixed: `fresh()` before every recheck, with tests that rechecks really solve. | The planned check was replaced by a stronger one, a different solver. Its criterion is the plan's: whether any call changes, with the value differences reported. Result: RECHECK_RESULT. |

Not deviations, but recorded:

- **Post hoc analyses** in the results are labelled as such:
  - the per-knockout share of readouts called increased;
  - the matching rates of true and other knockouts;
  - the uniquely-first counts;
  - the "always increased" baseline;
  - the Spearman correlation between adjusted rank and the true knockout's increase share.
- **Explanations corrected after the independent check:**
  - the adjusted-score confusers;
  - the profiles ranked in the bottom half;
  - number slips: 86% rather than 87%, 27.5 rather than 27 ties, "about two-thirds" rather than "almost doubles", and
    251 rather than 252 scored tuples.
