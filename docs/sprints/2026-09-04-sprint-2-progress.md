# Sprint 2 progress — whole-body models in Python

**4–5 September 2026 · Status: the whole-body solver work is complete; the full 63-IEM run is in progress after a bug fix (see §2). Numbers below are from completed solves.**

## 1. Harvey and Harvetta solve with open-source solvers

The COBRA Toolbox's test-model repository (github.com/opencobra/COBRA.models) ships Harvey 1.03d and Harvetta 1.03d as MATLAB files including their coupling constraints (`C`, `d`, `dsense`), which COBRApy's reader silently drops. `gembench/wbm.py` loads them directly (Harvey: 81,094 reactions, 56,452 metabolites, 103,600 coupling rows; Harvetta: 83,521 / 58,851 / 107,351), builds the LP with the coupling rows stacked under the stoichiometric rows, solves with HiGHS or GLPK, and certifies every solution by recomputing all residuals in double precision — a solver-independent check that the answer satisfies the model.

| Model | Problem | Solver | Time | Max residual (S rows / coupling rows) |
|---|---|---|---|---|
| Harvey | feasibility (shipped bounds) | HiGHS 1.15.1 dual simplex | 138 s | 3e-9 / 9e-8 |
| Harvey | feasibility | HiGHS interior point + crossover | 19 s | 2e-9 / 1e-7 |
| Harvey | feasibility | GLPK 5.0 dual simplex | 1,003 s | 6e-8 / 1e-5 |
| Harvey | max urinary phenylalanine exchange | HiGHS simplex | 66 s | 1e-9 / 5e-8 |
| Harvey | same | HiGHS interior point | 15 s | 1e-9 / 1e-7 |
| Harvey | same | GLPK dual simplex | 1,167 s | 3e-8 / 1e-7 |
| Harvetta | feasibility | HiGHS simplex | 194 s | 2e-9 / 1e-7 |
| Harvetta | feasibility | HiGHS interior point | 20 s | 2e-9 / 1e-7 |

All three solver runs on the phenylalanine problem return the same optimum, 0.011802906525, agreeing to twelve digits. The models as shipped carry bounds of ±1,000,000, coupling coefficients of 20,000 and stoichiometric coefficients down to 1e-6, so the coefficient range spans twelve orders of magnitude; HiGHS handles it at a 1e-7 feasibility tolerance without any scaling on our side, on one core of a two-core machine. This answers the roadmap's open question Q3 in the affirmative for the open-solver path: the lab's flagship models are solvable in Python without CPLEX or MOSEK, with interior point plus crossover the right default for cold starts (about twenty seconds) and dual simplex a poor choice for warm starts here (bound and objective changes led to multi-minute re-solves, so the IEM protocol below uses interior point for every solve).

What this does not yet cover: the physiological constraints (`physiologicalConstraintsHMDBbased`, `standardPhysiolDefaultParameters`), the diet (`EUAverageDietNew`, `setDietConstraints`) and the solver-specific scaling in `optimizeWBModel` that the Toolbox applies before any published analysis. Those files (about 2,000 lines of MATLAB plus the HMDB concentration, organ-weight and blood-flow input tables) have been pulled from the Toolbox and are the next port.

## 2. The IEM biomarker protocol in Python

`gembench/wbm_iem.py` ports `checkIEM_WBM` and the direction rule from `runIEM_HH.m`: maximise the summed flux through the IEM reactions in the healthy model; pin it there; close the reactions in the disease model; check the whole-body objective is still feasible; for each biomarker exchange or demand reaction, maximise it in both states; call increased or decreased when the difference exceeds 1e-6. The per-IEM protocol (reaction patterns to knock out, exclusions, per-block bound tweaks, biomarker reactions with expected direction) is extracted from the script into `iem_protocol_v0.json` by `scripts/parse_iem_protocol.py`, and the global reaction constraints set at the top of the script are applied. Demand reactions for blood biomarkers are added on the fly as in the original.

On the shipped bounds (no physiological or diet constraints yet), the IEMs completed so far give: histidinemia 5/5 biomarkers in the expected direction, AGAT deficiency 2/2, argininemia 6/6, CPS1 deficiency 3/6, 3-methylglutaconic aciduria 2/3, AMA1 2/2, AADC deficiency 8/10 — 28 of 34. The first full pass (5 September) silently skipped every IEM after the fourth: HiGHS applies its time limit to the cumulative run time of a persistent model object, so once the first four IEMs had used 30 minutes of solver time every later solve returned "time limit reached" in zero seconds; the protocol treated that as the Toolbox does (IEM not simulated) and moved on. The limit is now re-based before each solve and the remaining IEMs are re-running. Each IEM takes 7 to 13 interior-point solves at roughly 17 to 30 seconds each on this machine, so the full 63-IEM pass takes several hours here and is running in the background; its per-biomarker table and overall accuracy will be appended when it completes. The published figure (85 percent over 252 biomarkers) was obtained with the physiological and diet constraints in place, so agreement or disagreement with it is not yet meaningful — the comparison that matters comes after the constraints port.

## 3. Housekeeping

The cross-solver FROG pass stopped after iML1515 (four E. coli models done, all agreeing to 1e-6) and will be resumed when the machine is free. The E. coli and yeast benchmark results from Sprint 1 are unchanged. Everything is in the updated `gembench` archive.

## 4. Next

In order: port the physiological and diet constraints and re-run the IEM protocol under the published conditions; then read `optimizeWBModel` to document the Toolbox's solver settings for the conformance suite; then, once Tim's GitHub repository exists, move the code there and set up the Fitness Browser and DepMap ingestion.
