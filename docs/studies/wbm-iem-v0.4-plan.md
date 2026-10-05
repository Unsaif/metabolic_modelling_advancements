# Whole-body IEM protocol v0.4: Harvetta and flux ranges (plan)

Declared 5 October 2026 by Claude (Opus 5.5). No IEM simulation has been run on Harvetta, and no biomarker minimum has been computed on either model.

What is already known:
- the Harvey v0.3 results (217 of 251 scored directions correct) and the MATLAB reference run for Harvey ([v0.3 results](wbm-iem-v0.3-results.md));
- the published accuracies: Harvey 85.3% and Harvetta 84.9% (Thiele et al. 2020, other model and code versions).

This is a development study for the second paper: an open Python pipeline validated against the COBRA Toolbox, and how robust the IEM biomarker predictions are. The freeze fixes three things before the new numbers exist: what is run, how it is scored, and how the results will be read.

## 1. Questions

**Q1. Do the Harvey findings hold on Harvetta 1.03d?**
- (a) The port reproduces the Toolbox's model setup exactly.
- (b) Python and MATLAB make the same call wherever MATLAB returns an optimum.
- (c) `runIEM_HH`'s own accuracy figure is lowered by solver failures that it counts as "no change".
- (d) Most "no change" errors are biomarkers whose healthy and disease maxima are equal and positive (a cap).
- (e) Some correct calls rest on small relative differences.

**Q2. Do flux ranges fix the cap problem?** The protocol compares only maxima. Comparing the minimum as well might recover correct calls where both states reach the same maximum. The question is whether it does so without losing more correct calls than it gains.

## 2. Established before the freeze (no IEM outcome involved)

- **Model.** `Harvetta_1_03d.mat` from opencobra/COBRA.models commit `75c070d` (sha256 `98810012918a5df1e46499fe2bd59cfc9730a570c10d1d4e146ab82ef671fad3`), 83,521 reactions. Harvey 1.03d from the same commit has the hash recorded for v0.3.
- **Port check** (`results/wbm_iem/Harvetta_1_03d_constraint_port_check_v0.3.json`):
  - With the two older parameter values (GFR 129.75 ml/min, CSF export from 0.35 ml/min), the port reproduces every stored Harvetta bound except three blood–brain-barrier uptakes, as in Harvey.
  - Re-applying the current Toolbox code changes 1,000 lower and 562 upper bounds (Harvey: 1,001 and 571).
- **MATLAB setup.** Runner job `iem_ref_harvetta_20261005_step1_setup` ran on 5 October at 13:32 UTC on Tim's Mac: MATLAB R2024b, Gurobi, and none of the 19 checked Toolbox files differ from commit 67c790d. Its setup bounds equal the port's for all 83,521 reactions, lower and upper (`results/wbm_iem/matlab_reference/harvetta_step1_setup_bounds_comparison.json`). Q1(a) therefore holds already.
- **Code.**
  - `run_iem` gains a `senses` option ("max", "min"); `scripts/run_wbm_iem.py` exposes it as `--senses`.
  - The maximum path is unchanged: the tests compare it with and without minima, and a max-only run has the same provenance fields as v0.3.
  - `scripts/iem_range_calls.py` implements the range rules below, with tests.

## 3. Runs

Fixed for every Python run:
- **Protocol:** `data/iem/iem_protocol_v0.2.json` (57 IEMs, 252 biomarkers).
- **Scoring:** `checkIEM_WBM` logic as in v0.3 (minRxnsFluxHealthy 1, tolerance 1e-6, per-IEM state isolation).
- **Solver:** HiGHS 1.15.1, interior point with crossover, tolerances 1e-7, 1,800 s per solve.
- **Global constraints:** as in `runIEM_HH.m`, with the 28-reaction bile-duct list.
- **Model setup:** `toolbox`, meaning the current Toolbox constraints and EU average diet re-applied, as in v0.3.

| Run | Model | Optima | Command |
|---|---|---|---|
| **A. Harvetta v0.4** (Q1) | Harvetta | maxima (the protocol) | `python3 scripts/run_wbm_iem.py Harvetta_1_03d --model-setup toolbox --bile-duct toolbox --out-suffix _v0.4` |
| **B. Harvey minima** (Q2) | Harvey | minima | `python3 scripts/run_wbm_iem.py Harvey_1_03d --model-setup toolbox --bile-duct toolbox --senses min --out-suffix _v0.4_min` |
| **C. Harvetta minima** (Q2) | Harvetta | minima | as B, for `Harvetta_1_03d` |
| **M. MATLAB reference** (Q1) | Harvetta | `runIEM_HH` unchanged except the model file | queue job `iem_ref_harvetta_20261005_step2_runiem` (`tools/matlab/runner/iem_reference_harvetta_step2_runiem.m`) |

Notes on the runs:
- **Harvey maxima.** For Q2, the Harvey maxima are those of the frozen v0.3 run (`results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json`). Run B must report the same `lp_bounds_sha256` as v0.3 (`693af5c9…`).
- **Shards.** Runs B and C may be split into disjoint IEM shards with `--iems` and suffixes `_v0.4_min_a`, `_v0.4_min_b`. They are merged only if every record carries the same run fingerprint.
- **Order.** The two cores of the workspace run A and B first, then C.

## 4. Analyses and how they will be read

### Q1, on Harvetta, with Harvey v0.3 alongside

1. **Accuracy (A1).**
   - Correct calls among scored biomarkers, coverage, and errors split into "opposite direction" and "no change".
   - IEMs fully correct.
   - The published 84.9% is shown for context only. It is not a target.
2. **Python against MATLAB (A2), with `scripts/compare_matlab_reference.py --model Harvetta_1_03d`.**
   - Setup and global bounds; per-biomarker calls; non-finite MATLAB optima.
   - `runIEM_HH`'s formula, and accuracy among biomarkers with finite MATLAB values.
3. **Error anatomy (A3).** Among the "no change" errors, the number with equal positive maxima.
4. **Effect size (A4), with `scripts/iem_effect_size_sensitivity.py`.** Correct calls at relative thresholds τ = 0, 0.001, 0.01, 0.05 and 0.10. This was post hoc for Harvey; it is declared in advance for Harvetta.

Reading for Q1: each finding is reported as "replicated" or "not replicated" on Harvetta, with its numbers.

| Finding | Replicated if |
|---|---|
| (a) | the setup bounds are identical (already met) |
| (b) | the calls are identical for every biomarker where both MATLAB optima are finite; otherwise the differing biomarkers are listed |
| (c) | MATLAB returns at least one non-finite optimum for a directional biomarker, and `runIEM_HH`'s formula gives a lower figure than Python's accuracy among scored |
| (d) | more than half of the "no change" errors have equal positive maxima |
| (e) | descriptive only: the number of correct calls below τ = 5% (Harvey: 10 of 217) |

### Q2, both models

5. **Range calls (A5), with `scripts/iem_range_calls.py`.**
   - **Scored set:** directional biomarkers whose four optima (healthy and disease, maximum and minimum) are all Optimal and finite.
   - **R1, full range** (in the spirit of Shlomi, Cabili and Ruppin 2009). Compare the minima and the maxima of the two states, each with the protocol tolerance of 1e-6:
     - Increased if at least one rises and neither falls;
     - Decreased if at least one falls and neither rises;
     - Unchanged if both are equal;
     - Conflicting if one rises and the other falls, counted as wrong.
   - **R2, tie-break.** Keep the protocol call unless it is Unchanged; in that case call from the minima alone.

Reading for Q2, per model:
- A rule **improves** if it has more correct calls than the maximum-only protocol on the same scored set, **worsens** if it has fewer, and shows **no net change** otherwise.
- **R1** is **supported** if it improves on both models, **harmful** if it worsens on either, and **inconclusive** otherwise.
- **R2** can change only Unchanged calls, so it cannot lose a correct call. It is **useful** if, on both models, it turns more errors into correct calls than into opposite-direction calls.
- Also reported:
  - the transitions;
  - Conflicting calls;
  - results by biofluid;
  - the subset with equal positive maxima.

## 5. Rules

- **No changes after the freeze** to the protocol, the port, the solver settings, the scoring or the reading rules.
- **Technical failures.** A killed process resumes from its results file under the same fingerprint, and shards are merged only under one fingerprint.
- **Bugs.** A bug found after the runs start is reported. A corrected run gets a new suffix and is shown alongside the original, never instead of it.
- **MATLAB.** Job M is submitted once. If it fails or times out, its records are inspected; it is not resubmitted automatically. Any second attempt uses a new filename and is reported.
- **Verification.** A separate agent recomputes every reported number from the result files.
- **Time order.** The plan commit, the freeze manifest and an external record in the Claude project come before runs A, B, C and M start.

## 6. Limits

- **Not blind.** The IEM labels, the Harvey results and the published Harvetta figure are known. Harvetta shares structure, labels, diet and constraints with Harvey, so Q1 is a replication on a sibling model, not an independent validation.
- **The range rules were motivated by Harvey's capped errors.** They are fixed in code before any minimum exists on either model, so on both models they are tested on values nobody has seen.
- **Directions only.** The ground truth is the protocol's literature labels, with one diet and one model version per sex.
- **Small numerical differences.** Repeated interior-point solves can differ by about 1e-5 relative, as the v0.3 certification showed. Calls that rest on differences of that size are flagged by A4.
