# Fixed numerical comparison after the curated-model failure

6 September 2026. Development software conformance, explicitly post-failure.

The v2 medium-preparation run stopped at the first curated iJN1463 condition
because the returned HiGHS simplex vector exceeded the original absolute
mass-balance tolerance. The primary run did not preserve that failed vector.
A separately frozen one-condition diagnosis reproduced the failure. A bounded
follow-up found that primal simplex and tighter default-simplex tolerances could
meet the unchanged acceptance gate on that condition. These observations motivate
this comparison; it is not a prospectively independent numerical benchmark.

The complete recipe is [medium_curated_numerics_v1.json](../../data/studies/medium_curated_numerics_v1.json).
It fixes all 43 mapped Putida carbon conditions for the curated iJN1463 model,
in the frozen source inventory's order. Every condition receives all three methods,
giving exactly 129 solves:

| Method | Settings |
|---|---|
| Native GLPK | Original feasibility tolerance 1e−9; optlang's existing auto-scaling and other native settings |
| HiGHS primal simplex | Strategy 4; primal/dual feasibility tolerances 1e−9 |
| HiGHS default simplex | Original strategy; primal/dual feasibility tolerances tightened to 1e−10 |

Both HiGHS methods use presolve `choose`, parallel `off`, automatic thread count
and a 60-second limit, as explicitly recorded. Supported option values are checked
before optimization. Native GLPK effective tolerances are saved with each result.
The [official HiGHS options](https://ergo-code.github.io/HiGHS/dev/options/definitions/)
document strategy 4 as primal simplex. Changing strategy does not guarantee that
large internal circulations disappear; they persisted in the one-condition tests.

Preparation uses the frozen v2 helper and the original iJN1463 source loader.
The existing `EX_AEP_e -> 2ameph_e` identity exception is preserved. The complete
medium union, exclusions, stoichiometry, reaction bounds, objective, GPRs and
chemical metadata must exactly match the historical fresh-condition formulation
before each case is solved. There is no reaction repair, new nutrient, changed
bound or adaptive solver setting. No numeric fitness values are parsed.

Every raw returned method record is written before post-solve algebra checks or
acceptance decisions. Optimal vectors that fail numerical acceptance remain
available, as do time limits, infeasibility and exceptions. All cases continue
after numerical rejection; unsupported inputs, altered LPs or changed frozen files
stop the run. Each completed case includes the full prepared signature, medium
report, all three solver records and a separate assessment.

The original sparse primal gate remains separately identifiable: finite maximum
absolute mass-balance residual and bound violation must both be at most 1e−8.
The new `accepted_for_comparison` field additionally requires an optimal status,
a complete finite reaction vector, finite compensated-summation residual at most
1e−8, and consistency of the reported objective with that vector. Objective
comparisons retain the original absolute and relative tolerance of 1e−8.

Cross-solver evidence requires accepted native GLPK and at least one accepted
HiGHS record. **All reported optimal objectives must agree, including those from
rejected primal vectors.** An optimal classification conflicting with infeasible,
unbounded or an ambiguous infeasible-or-unbounded status prevents that evidence
label. A third method's timeout or capture exception can coexist with two accepted
methods, but is always reported separately and never counted as accepted.
Mutual infeasibility remains an explicit status, never a fabricated zero growth
rate. Primal certificates and objective agreement are numerical evidence; they
are not a mathematical optimality proof or validation of the biochemical model.

The proposal, exact input chain, all completed numerical-diagnostic records and
method settings are committed and pushed before execution. The runner requires
a clean checkout and writes a new manifest before its first solve. Explicit file
lists avoid absorbing active checker outputs. It refuses an existing output
directory, rechecks frozen inputs after completion, and retains failed attempts.

The original v1/v2 primary attempts remain incomplete and failed. This additive
comparison does not rename them as successful, revise phenotype scores or supply
independent biological validation. Its outcome will be reported for the full
declared panel, without further method selection inside the run.

After the proposal has been committed and pushed, reproduce in a fresh directory:

```sh
.venv/bin/python scripts/compare_curated_medium_solvers.py --out <fresh-directory>
```
