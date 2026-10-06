<!-- Copied on 6 October 2026 from Tim's MATLAB job runner (/Users/timhulshof/Documents/matlab-agent-runner/VERIFICATION.md), written by the agent that set the runner up. It is the record of the MATLAB, COBRA Toolbox and Gurobi versions used for the MATLAB reference runs. -->

# Local verification — 5 October 2026

All checks below passed using the restored image and existing licences. No model API calls were made. These are engineering checks, not a new reconstruction experiment or a scientific acceptance decision.

| Check | Elapsed time | Result |
|---|---:|---|
| First MATLAB/COBRA/Gurobi start | 11.264 s | Passed |
| COBRA LP and direct Gurobi LP | 0.331 s | Passed |
| Load the frozen reconstruction | 0.872 s | Passed |
| Baseline optimizeCbModel call | 0.590 s | Passed |
| Intentional error and pre-error output capture | 0.333 s | Expected failure correctly reported |
| State persists after error; file saved | 0.225 s | Passed |
| Intentional timeout and cleanup | 1.375 s | Expected failure correctly reported |
| Confirm stopped session | 0.083 s | Passed |
| New ready session for handoff | 9.550 s | Passed |

Eight offline lifecycle tests also passed. The warm requests shared one MATLAB process. Both small optimisation problems returned the known optimum of 1.

The loaded model contains 13,637 metabolites, 22,417 reactions and 3,402 genes. The baseline call returned optimal solver status with objective 0; this checks execution only and does not establish growth or physiological performance. The original source model hash was unchanged.

MATLAB reported R2024b Update 9 (24.2.0.3212159). Gurobi 12.0.0 and COBRA revision 67c790dbac809d9d891fdbafc33e18c21fc9bddc were restored from the existing pinned recipe.

The final session was left ready, with an empty user workspace. It exits after 30 minutes idle. Use start if it has expired. The timeout test used the earlier session; no failing job was retried.

Image ID: `sha256:0ff43ebcfadc780d7642a1ed7304991c50f371b88a2767c574c2b1d2ad663c12`

Source model SHA-256: `bdd361f148a57153f2f4440bb112d8fd00ac3c85abb2eea131d252d260029772`

Machine-readable measurements: [verification.json](verification.json). Detailed execution records are retained under this runner’s runtime folder and the project’s outputs/matlab-runner-2026-10-05 folder.

## Background queue verification

The installed login service started MATLAB automatically and processed three file submissions in sequence: a COBRA optimisation with a 14,400-second timeout setting, an intentional MATLAB exception, and a following successful script. Result JSON appeared next to each script, and generated output was readable. The short jobs took about 0.20–0.26 seconds each; this did not test four hours of continuous execution.

Restarting the watcher preserved the results and did not rerun the submissions. All 23 offline runner/queue tests passed, including timeout validation, busy-session handling, crash recovery and duplicate prevention. No model API calls were made.

[Queue verification measurements](queue-verification.json).
