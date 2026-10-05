# MATLAB reference run of the IEM protocol (Harvey 1.03d)

Purpose: run the COBRA Toolbox's own `runIEM_HH.m` on Harvey 1.03d, then compare its bounds and biomarker predictions with the Python port. The Python port is `gembench/wbm_constraints.py` and `scripts/run_wbm_iem.py` (run v0.3, plan `docs/studies/wbm-iem-v0.3-plan.md`).

## What the script changes

`run_iem_reference_harvey_1_03d.m` runs a copy of `runIEM_HH.m` with exactly two changes:

1. **Model loading.** It loads `./Harvey_1_03d.mat`. Unmodified, `runIEM_HH.m` loads the newest Harvey it finds on the MATLAB path.
2. **Stray command.** It drops the stray `edit` command at the start of line 9.

It also checks:
- the model file's sha256;
- the Toolbox files the protocol depends on, against commit `67c790d`;
- the LP solver in use.

It writes all three results to `out/run_info.txt`.

## Run

Put this folder (the two `.m` files and `Harvey_1_03d.mat`) anywhere, then run:

```
matlab -batch "initCobraToolbox(false); changeCobraSolver('gurobi','LP'); cd('<this folder>'); run_iem_reference_harvey_1_03d"
```

You can also run the three statements at the MATLAB prompt. The script must run in the base workspace, because `runIEM_HH` reads its results back with `evalin('base', ...)`. Do not wrap it in a function.

## Outputs (`out/`)

| File | Contents | When |
|---|---|---|
| `run_info.txt` | Toolbox commit, file checks, solver, MATLAB version, timings | throughout |
| `matlab_setup_bounds_Harvey_1_03d.mat` | `rxns`, `lb`, `ub` after the diet and physiological setup | minutes |
| `matlab_global_bounds_Harvey_1_03d.mat` | bounds after the unified reaction constraints (`modelO`) | end |
| `matlab_iem_results_Harvey_1_03d.mat` | `IEMSol_*` of all IEMs, `Table_IEM`, `Accuracy`, ... | end (hours) |
| `Results_IEM_Harvey_1_03.mat` | the workspace `runIEM_HH` saves itself | end; large, not needed |

## Compare

Copy `out/` back into the repository. Then run:

```
python scripts/compare_matlab_reference.py --matlab-dir <path to out> --out results/wbm_iem/matlab_reference_comparison.json
```

The comparison checks two things:
- **Bounds.** Every bound against the Python setup; with the same Toolbox version they should be identical.
- **Biomarker calls.** Every biomarker call against the Python v0.3 results. MATLAB calls are made from the printed values, as in `runIEM_HH`.
