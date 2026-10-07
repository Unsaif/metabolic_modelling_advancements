# Main session's response to the results check

7 October 2026.

| # | Problem | Response |
|---|---|---|
| 1 | Rechecks were no-ops | **Accepted, with an apology to Tim, who had been told they matched.** Both backends now have `fresh()`, called before every recheck (`model.reset(0)` for Gurobi; `clearSolver()` for HiGHS), and tests check that rechecks do solver work. The same 197 LPs were re-solved from scratch with HiGHS, a different solver, by `scripts/recheck_cross_matrix.py`. The deviation is recorded in `docs/studies/wbm-iem-ranking-deviations.md`, and the check table is reworded. |
| 2 | Adjusted-score confusers explained wrongly | **Rewritten from the data.** Only AADC, TETB, LTC4S and FED change fewer than half the readouts. The others change about the median number, and outrank true diseases whose own knockouts raise even more, which the adjustment discounts more. FED is added to the list. |
| 3 | "Ranked below most knockouts" | **Rewritten** as "ranked in the bottom half": MMA, HCYS, FED, HLYS2 and DPYR. The own-miss explanation is applied only to MMA, HCYS and DPYR; FED and HLYS2 miss half. |
| 4 | Number slips | **Fixed:** 86%, 27.5 others, "about two-thirds", 251 scored tuples, the relative difference qualified as "among values above 1e-3", and "21 of the 39 profiles have only one or two readouts". |
| 5 | Uneven gain from the adjustment | **Added:** Spearman 0.70 (0.06 plain), 27 profiles better and 22 worse with examples, and the increase shares (1% to 19%) of the 8 all-increase profiles that become uniquely first. |
| 6 | Mechanism | **Restated.** What follows from the protocol is now stated as such: zero lies within the bounds of every IEM reaction, so without the pin no maximum could rise. Only the reason for the loss on unrelated readouts is untested. A partial pin or a typical healthy state is given as the test. |
| 7 | Context | **Added to the limits:** the LTC4S maximum flux without sinks, the unchanged own calls of the three IEMs tested, and that the ranking without sinks is untested. |
| 8 | Post hoc labels | **Added** to the summary bullet. |
