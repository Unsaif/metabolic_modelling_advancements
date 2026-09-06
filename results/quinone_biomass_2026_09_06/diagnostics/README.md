# Post-run quinone diagnostic

These checks diagnose the completed development study. They do not select a biological patch or establish independent validation. The frozen study script, inputs, and matrices were unchanged.

Run from the repository root, choosing a new output directory:

```sh
.venv/bin/python scripts/diagnose_quinone_producibility.py --out /tmp/quinone-diagnostic-new
```

The script verifies the main input hashes before computing. It refuses an existing output directory; a repeated invocation was confirmed to fail without changing existing output bytes. `producibility.json` records the script hash, main manifest hash, solver versions, primal residuals, and the post-run hash of the source-card metadata file omitted from the main freeze.

The combined `q8_c` and `q8h2_c` mass-balance rows sum to exactly zero in the prepared baseline. Ten reactions interconvert this pool; none supplies it. Adding either positive biomass requirement therefore forces zero growth by mass balance. All 42 diagnostic solves agree with the expected outcomes across GLPK at feasibility tolerances 1e-7 and 1e-9 and HiGHS at 1e-9. Baseline glucose growth is 0.9193705348; both ubiquinol demands block growth, and an ideal diagnostic source restores baseline growth. Demanding growth of at least 0.1 is infeasible. The maximum mass-balance residual among the 36 optimal solutions is 1.13e-11. The ideal source is a diagnostic intervention, not a proposed biochemical reaction.

`independent_summary.json` independently recomputes the stored confusion counts and MCCs. The baseline has 34/43 growing conditions; both ubiquinol arms have none; the restored menaquinol control grows only on benzoate. Its three changed gene predictions remain zero under all diagnostic solvers and tolerances.

Prior matrices are aligned by exact gene and condition identifiers before comparison. The true development sprint cycle7 baseline reproduces all binary predictions, all wild-type rates exactly, and knockout rates within 4.85e-11. The incomplete first-attempt baseline also has identical binary predictions. The model-v0.3 archive is correctly identified as cycle6 and differs at 135 binary knockout predictions.

Two completed diagnostic attempts are preserved in sibling `diagnostics_attempt01` and `diagnostics_attempt02` directories with their source snapshots. Attempt 1 mislabeled the cycle6 archive as cycle7; attempt 2 corrected the label but compared axis order without aligning labels. Neither issue affected the physical model checks. This directory contains the final aligned comparison.
