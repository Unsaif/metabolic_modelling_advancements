# IEM disease ranking: deviations from the plan

Plan: `wbm-iem-ranking-plan.md`. Amendments made before the run are in the plan itself and its git history
(347b21f, 34ab349, a6e98ef, def01b0). Listed here are the changes made after the run started (6 October 2026,
19:37 Irish time).

| # | Deviation | When decided | Why | Effect |
|---|---|---|---|---|
| 1 | The 197 planned rechecks (barrier from scratch, about 1 in 100 warm solves) were not performed during the run. After the ranking, the same 197 LPs were re-solved from scratch twice by `scripts/recheck_cross_matrix.py`: with Gurobi 13.0.1 barrier as planned (Tim's Mac), and with HiGHS 1.15.1 interior point with crossover (cloud). Two HiGHS rechecks that were infeasible with the matrix's pin were retried with HiGHS's own pin. | 7 October, after the ranking was computed. Found by the independent check of the results. | The runner called the barrier solve without resetting the model. Gurobi then returned the solution it already had (0.000 s each, difference exactly 0). The fault is fixed: `fresh()` before every recheck, with tests that rechecks really solve. | The planned check was done late, and a different solver was added. The criteria are the plan's: agreement within 1e-6 (absolute, or relative above 1), and whether any call changes. **No call changed in either recheck**, so no IEM had to be recomputed. Gurobi: all 197 optimal, 189 within the tolerance; the other 8 are healthy-state maxima below 0.001 that differ by up to 1.6e-5. HiGHS: 195 optimal at the first attempt, and the 2 PC healthy states after the retry with HiGHS's own pin; 187 within the tolerance, the other 10 healthy-state maxima below 0.001 that differ by up to 4.2e-5. |

Not deviations, but recorded:

- **Post hoc analyses** in the results are labelled as such:
  - the per-knockout share of readouts called increased;
  - the matching rates of true and other knockouts;
  - the uniquely-first counts;
  - the "always increased" baseline;
  - the Spearman correlation between adjusted rank and the true knockout's increase share.
- **Observed in the HiGHS recheck:** the healthy pin, the maximum truncated to six decimals, leaves too small a margin
  to be reused across solvers. PC's pin from Gurobi sits 4.5e-7 above HiGHS's maximum (see the results' limits).
- **Explanations corrected after the independent check:**
  - the adjusted-score confusers;
  - the profiles ranked in the bottom half;
  - number slips: 86% rather than 87%, 27.5 rather than 27 ties, "about two-thirds" rather than "almost doubles", and
    251 rather than 252 scored tuples.
