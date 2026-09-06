# Numerical sensitivity of growth-ranking scores

A fresh rerun of the exact saved quinone-demand parent reproduced every binary
growth prediction. Wild-type growth, fitness values and all matrix labels were
identical. Mutant growth values differed by at most **6.002 × 10⁻¹¹**. Nevertheless,
the standard average-precision score changed from **0.406846 to 0.401588**.
Its ranking uses negated absolute mutant growth against measured fitness labels.
Tiny differences can reorder tied or nearly tied predictions without changing
any growth/no-growth decision. The primary MCC was unchanged at 0.574384.

This observation prompted a separate **post-outcome numerical diagnostic**, not a
change to the frozen experiment or its reported scores. The diagnostic compares
two policies over the fixed grid 10⁻¹², 10⁻¹⁰, 10⁻⁹, 10⁻⁸, 10⁻⁷ and 10⁻⁶: setting
values within that distance of zero to zero, and rounding all finite growth
values to that quantum. Missing and nonfinite values are preserved.

| Treatment | Previous parent AP | Fresh parent AP |
|---|---:|---:|
| Raw saved predictions | 0.406846334 | 0.401587529 |
| Zero band of 10⁻⁹ | 0.416703138 | 0.416871963 |
| Rounding quantum of 10⁻⁹ | 0.414810130 | 0.414810130 |

Rounding at 10⁻⁹ also makes AUROC equal (0.761461358) for the two saved runs.
Zero-only treatment leaves some sensitivity among nonzero predictions. None of
the 24 transformed run/policy/precision combinations changes a binary prediction.
These examples demonstrate numerical rank fragility; they do not identify a
biologically correct precision or justify selecting whichever score is higher.

The raw arrays, primary cards and intervals remain unchanged. The sensitivity
report contains point estimates only; it does not recompute bootstrap intervals.
The source code, input hashes and all grid results are preserved under
[`analysis/rank_stability_parent/`](../../results/ppnp_repair_2026_09_06/analysis/rank_stability_parent/).
The direct array comparison is
[`parent_reproduction.json`](../../results/ppnp_repair_2026_09_06/analysis/parent_reproduction.json).

Before relying on small ranking-score differences in future studies, define and
version a precision policy from solver accuracy and scientific resolution, test
it on development examples, and freeze it before independent evaluation. A policy
must also address values near the classification threshold; automatic rounding
is not universally guaranteed to preserve classifications. Report binary and
ranking metrics separately. Numerical reproducibility cannot establish biological
validity.

Reproduce this diagnostic in a fresh directory with
`scripts/audit_growth_rank_stability.py --runs <previous-parent> <fresh-parent> --out <new-directory>`.
