# Whole-body IEM protocol v0.3: plan

Declared 4 October 2026 by Claude (Opus 5.5), before any v0.3 or v0.2b simulation was run. The v0.2 results (220 of 251 scored biomarker directions correct) are known, so this is not a blind evaluation. The point of the freeze is narrower: the protocol, the constraint port and the solver settings are fixed before the new numbers exist, so they cannot be adjusted towards the published 85 percent or towards v0.2.

## 1. Why a v0.3

Porting the model preparation of `runIEM_HH.m` (Harvey branch: `standardPhysiolDefaultParameters`, `physiologicalConstraintsHMDBbased`, `EUAverageDietNew`, `setDietConstraints`) gave three findings.

1. **Harvey 1.03d already carries the diet and physiological constraints.** Its `SetupInfo` records the EU average diet. When the port is run with two parameter values from the older Toolbox function (`physiologicalConstraintsHMDBbased_old.m`), it reproduces the stored bound of every reaction except three. The two older values are a glomerular filtration rate computed as 20 percent of renal plasma flow (129.75 ml/min, where the current code uses the parameter 90 ml/min) and CSF export upper bounds from the CSF flow rate (0.35 ml/min, where the current code uses the CSF-to-blood flow of 0.52 ml/min). The three exceptions are the blood-brain-barrier uptakes of tryptophan, kynurenine and kynurenate, which are unconstrained in the stored model. **The v0.2 caveat ("shipped bounds without the physiological constraints and the EU average diet") was therefore wrong.** v0.2 ran with the constraints as they were set when the model was released.
2. **Re-applying the current Toolbox code changes 1,572 bounds.** `runIEM_HH.m` re-applies these steps after loading the model, so this is what the Toolbox does today:
   - 1,011 kidney bounds scale by 0.694 (GFR 90 instead of 129.75 ml/min);
   - 517 CSF export upper bounds scale by 1.486 (0.52 instead of 0.35 ml/min);
   - the three BBB uptakes become constrained;
   - 31 diet uptakes open at 0.1 mmol/day (AGORA-essential metabolites added to the Toolbox list in 2022 and 2024);
   - 5 bile acids enter the diet at 10 mmol/day ±20 percent.

   Source: `results/wbm_iem/Harvey_1_03d_constraint_port_check_v0.3.json`.
3. **v0.2 deviated from `runIEM_HH.m` in one global step.** The Toolbox sets the upper bound to 100 on a list of 28 bile-duct exits (`Rnew`). The v0.2 code set it on every `BileDuct_EX_*[bd]_[luSI]` reaction (261 in Harvey). That opened 226 exits that are closed in the model and raised 7 more.

The port was checked by an independent agent that wrote its own literal re-implementation from the Toolbox files. Its bounds were bit-identical to the port's on Harvey and Harvetta, from stored, relaxed and randomised starting bounds. It found only latent differences (duplicate identifiers in input tables); these were fixed to the MATLAB first-match behaviour.

## 2. Runs

Fixed for both runs:
- **Model:** `Harvey_1_03d.mat` (sha256 `10e6cb6d02736fb88ae266e0c901e90766176bb80db2753e694e3a9330368c9c`).
- **IEM protocol:** `data/iem/iem_protocol_v0.2.json`, unchanged (57 IEMs, 252 biomarkers).
- **Scoring logic:** the `checkIEM_WBM` logic of v0.2: minRxnsFluxHealthy 1, tolerance 1e-6, per-IEM state isolation.
- **Solver:** HiGHS interior point with crossover, tolerances 1e-7, threads 0, 1,800 s per solve.
- **Global constraints:** as in `runIEM_HH.m`, with the 28-reaction bile-duct list.

| Run | Model bounds | Command |
|---|---|---|
| **v0.3 (primary)** | Shipped model with the constraints and the EU average diet re-applied by the port (`gembench/wbm_constraints.py`; inputs `data/iem/wbm_constraint_inputs_v0.3.json` from Toolbox commit 67c790d). This is what `runIEM_HH.m` does at that commit. | `python3 scripts/run_wbm_iem.py Harvey_1_03d --model-setup toolbox --bile-duct toolbox --out-suffix _v0.3` |
| **v0.2b (secondary)** | Stored bounds (the constraints as released with the model), with the corrected bile-duct step | `python3 scripts/run_wbm_iem.py Harvey_1_03d --model-setup shipped --bile-duct toolbox --out-suffix _v0.2b` |

v0.2 → v0.2b isolates the bile-duct correction. v0.2b → v0.3 isolates the re-applied constraints.

## 3. What will be reported

- For each run:
  - IEMs simulated, by status;
  - biomarkers scored (optimal, finite healthy and disease values);
  - correct calls, accuracy among scored, and coverage;
  - errors split into "no change predicted" and "opposite direction".
- Per biomarker and per IEM, compared with v0.2: the same call, or a changed call (and whether it became right or wrong).
- The published figure (85 percent, Thiele et al. 2020, with the model version and code of that time) is shown for context only. It is not a target.

## 4. Rules

- After this freeze, nothing changes in the protocol, the port, the solver settings or the scoring.
- **Technical failures.** A killed process resumes from its results file under the same run fingerprint. IEMs without a terminal status are attempted again, as in v0.2.
- **Bugs.** A bug found after the runs start is reported. A corrected run gets a new suffix and is reported alongside the original, never instead of it.
- **Verification.** A separate agent recomputes every summary number from the results files and checks a sample of records against the protocol.

## 5. Limits

- **No MATLAB here.** The port is checked against the stored model bounds and an independent re-implementation, not against a MATLAB run. One run of `runIEM_HH.m` in MATLAB by the lab would be the decisive check.
- **Model version.** `loadPSCMfile('Harvey')` loads the newest Harvey on the MATLAB path, which may be newer than 1.03d.
- **Toolbox quirks, reproduced as they are.** Two urine-table entries in other units are used as µmol/mmol creatinine: pyridoxine (nmol/mmol creatinine) and creatinine (µM). The urinary creatinine parameters (0.5–1.2 mg/dL) are plasma-like values.
- **The 85 percent figure** comes from a different model version and code. A like-for-like comparison needs the lab's MATLAB run.
