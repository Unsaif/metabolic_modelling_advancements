# IEM disease ranking (Harvey 1.03d): result files

Plan: `docs/studies/wbm-iem-ranking-plan.md`. Results: `docs/studies/wbm-iem-ranking-results.md`. Deviations:
`docs/studies/wbm-iem-ranking-deviations.md`.

| File | What it is |
|---|---|
| `../Harvey_1_03d_iem_cross_ranking_v1.json` | The prediction matrix: 57 IEMs × 187 readouts, healthy and disease state, Gurobi on Tim's Mac. Joined from the three `_shardKof3` files beside it. |
| `Harvey_1_03d_ranking_v1_checks.json` | The checks before ranking (`scripts/iem_ranking_checks.py`). **Its `rechecks` block is void:** those 197 rechecks did not re-solve (deviation 1). |
| `Harvey_1_03d_ranking_v1_recheck_gurobi.json` | The same 197 LPs re-solved from scratch with Gurobi barrier after the ranking, as the plan specified (`scripts/recheck_cross_matrix.py --backend gurobi`, Tim's Mac). |
| `Harvey_1_03d_ranking_v1_recheck_highs.json` | The same 197 LPs re-solved from scratch with HiGHS interior point (`--backend highs`, cloud). Two PC healthy-state LPs, infeasible with the matrix's pin, were retried with HiGHS's own pin (`--retry-infeasible-pins`, recorded under `first_attempt`). |
| `Harvey_1_03d_ranking_v1.json` | The ranking: all 16 analyses (4 profile sets × protocol or material calls × plain or adjusted score), per profile and in summary (`scripts/iem_disease_ranking.py`). |
| `run_logs/` | The three shard logs of the matrix run. |
| `independent_check_inputs/` | The two independent checks of the inputs before the run (report, response and the checker's scripts). |
| `independent_check_results/` | The independent check of the results (report, response and the checker's scripts). |
