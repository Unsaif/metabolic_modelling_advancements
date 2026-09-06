# Medium preparation that does not depend on prior use

6 September 2026. Software conformance using exposed development models.

The previous precursor audit exposed a concrete software risk: some nutrient
exchanges disappear from COBRA's inferred exchange list when a model mixes
compartment labels. The old helper can then leave nutrients available after
switching to a medium that excludes them. A tiny synthetic model demonstrates
the consequence: growth persists after the sole nutrient should have been
withdrawn. The new explicit preparation path prevents this dependency on history.

This is a software correction, not a new biological pathway or a fitness-score
improvement. The inspected historical runners used fresh preparation and per-condition state
restoration, which avoid the demonstrated sequential-use leak. Their frozen
inputs remain unchanged. The new API is opt-in; the historical benchmark entry
points still use their original implementation and require an explicit future
migration before receiving this behavior.

## Why the old behavior fails

The draft models primarily use `C_e`/`C_c` for extracellular/cytoplasmic
compartments. Earlier patches and medium completion also introduce `e`/`c`.
COBRA's exchange classification uses compartment inference and annotations.
Consequently, iterating only over `model.exchanges` can miss actual environmental
boundaries. Adding enough metabolites in the second compartment can also change
which compartment the heuristic chooses.
[Official boundary API](https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/medium/boundary_types/index.html).

The [source inventory](../../results/medium_preparation_2026_09_06/evidence/medium_structure_inventory.md)
covers the exact saved Putida PpnP parent, reconstructed historical cycle-7
Btheta/MR1/Smeli parents, and the existing Putida iJN1463 reference. The raw
drafts already contain omitted exchanges; their initial uptake bounds are zero.
The problem demonstrated here concerns subsequent state handling. iJN1463 uses
consistent `c`/`p`/`e` compartment labels. Keio/iML1515 is absent locally and
was not silently substituted or downloaded.

Other synthetic cases expose related errors: a positive-oriented exchange can
be turned into an unintended source; malformed exchange annotations can cause
internal reactions to be reset; a deleted exchange can be reopened while its
gene remains inactive; and a colliding carrier name can be silently ignored.
Invalid late inputs can also leave the caller's model partially modified.
These are demonstrated API failure modes, not claims that each occurred in
the historical biological evaluations.

## What the new implementation guarantees

`gembench.medium_preparation_v3.prepare_medium` validates a declared canonical
BiGG exchange model and returns an isolated copy plus a report. It:

- Requires explicit compartment roles and records approved alias normalization.
- Validates single-metabolite, coefficient −1 exchanges and their annotations.
- Resets every validated environmental exchange independently of COBRA's
  compartment heuristic and records every boundary bound.
- Rejects malformed nutrient bounds, inactive gene flags, conflicting completion
  carrier definitions and known formula/charge disagreements between metabolite
  pairs inspected by the completion recipe. Other existing transport reactions
  are not chemically validated by this helper.
- Reports missing components only when the caller explicitly allows this;
  otherwise it raises without changing the input model.
- Retains and discloses internal demands/sinks, including any able to supply
  material, and reports exchanges lacking a network connection.

“Boundary” here means COBRA's single-metabolite boundary definition. This report
does not claim to identify every artificial multi-metabolite source or sink.

The export cap, uptake-only carriers and excluded completion components remain
explicit assumptions. A matching formula or placeholder charge does not validate
an annotation. Exchange presence does not prove physiological importability.
Prepare media before applying gene deletions. The helper does not optimize,
infer transport biology or convert concentrations into uptake capacities.
[Official medium documentation](https://cobrapy.readthedocs.io/en/latest/media.html).

## Declared test and interpretation

The [definition](../studies/medium-preparation-v1.md) and recipe specify all
156 mapped conditions across five model preparations. Proposal `64d03eb` and
the added withdrawal-artifact records in `dc23b1f` were pushed before model
optimization. The clean manifest at `dc23b1f` freezes 148 inputs with fingerprint
`d6424a04a0a2fbaf6fd0a6eaf2197ac698c33915995aae24ca9f444c22ca1981`.

That first run completed 113 draft-model cases (226 solver vectors), then
correctly refused a naming exception in the curated model before solving it.
iJN1463's `EX_AEP_e` contains `2ameph_e`. An exhaustive source check finds only
this mismatch among its 348 EX reactions. The single metabolite, coefficient −1
and external compartment are unambiguous; neither identifier is requested by
the mapped carbon conditions or base media.

The [v2 definition](../studies/medium-preparation-v2.md) declares the exact
existing-reaction binding `EX_AEP_e -> 2ameph_e` for iJN1463 only. Nothing is
renamed, supplied or added by this binding. Its tested API extension and the
preserved first attempt were pushed as `69355bb` before the complete rerun.
The v2 manifest freezes 159 inputs from a clean checkout, fingerprint
`f9a4464851589820ffd524d9ccaf39be51ff0bd159c1e13724f7281dff546e97`.
The incomplete v1 run remains separate and has no completed panel summary.

For each fresh condition, compare all stoichiometry, reaction bounds, objective
coefficients, gene rules and formula/charge metadata exactly before solving.
Compartment labels and descriptive names are the intended representational
changes. Equality of the LP and GPRs entails identical reaction/gene-deletion
problems under matching subsequent restrictions; it does not rely on measuring
only one coincident growth rate. Native GLPK and separately translated HiGHS
then check the new LP and save complete primal vectors.

The separate withdrawal test opens newly completed components in a synthetic
stress medium and removes them again. It checks boundary limits, not growth or
experimental fitness. The stress medium must not be presented as an actual
experimental recipe or evidence of physiological nutrient carryover.

## Structural results and the retained numerical failure

The v2 attempt again saves all 113 draft conditions and 226 accepted solver
vectors, then stops at iJN1463's first condition, L-arginine in MOPS. Its
`failure.json` records a numerical certificate failure. Neither attempted
primary panel has a completed summary; neither is counted as 156 solved cases.

A separate [structural completion](../../results/medium_preparation_2026_09_06/structural_completion/report.json)
disables optimization and finishes the declared comparisons:

| Preparation | Fresh conditions with identical LP/GPR/chemistry | Old helper retains withdrawn supplements in synthetic test | New helper |
|---|---:|---:|---|
| Putida PpnP candidate | 43 | 5 | All withdrawn |
| Btheta cycle 7 | 25 | 2 | All withdrawn |
| MR1 cycle 7 | 12 | 5 | All withdrawn |
| Smeli cycle 7 | 33 | 4 | All withdrawn |
| Putida iJN1463 | 43 | 0 | All withdrawn |

All five new reuse preparations equal fresh preparations. The four draft
prefixes have 34/16/11/27 conditions above the existing growth threshold,
respectively. These are software conformance counts, not new measured phenotypes.
Exact equality compares numbers, so an integer bound `0` and floating-point
bound `0.0` are equal even though their serialized JSON hashes can differ.

The [one-condition numerical diagnosis](../../results/medium_preparation_2026_09_06/numerical_diagnostic/summary.json)
reproduces the failure and preserves its newly returned rejected vector; the
primary run's failed vector was not saved. It predeclares four numerical methods. Native
GLPK has maximum mass-balance residual 7.13e-14. HiGHS simplex reports an optimum
within 3.00e-9 of GLPK but has residual 7.18e-8, above the unchanged 1e-8 gate.
An independent compensated sum also rejects two metabolite rows. Its solution
uses flows at the inherited magnitude 999999, with about 1.60e7 total absolute
terms in the proton balance. This identifies large cancelling flows as a
numerical concern; it does not establish a biochemical error from that vector.
GLPK with an advanced basis also passes; HiGHS interior point times out at the
declared 60 seconds. The failed primary run remains preserved.

The bounded [four-option follow-up](../../results/medium_preparation_2026_09_06/numerical_diagnostic_attempt02b/summary.json)
tests primal simplex, dual simplex without presolve, tighter default-simplex
tolerances, and simplex without scaling. Primal simplex at the original solver
tolerance passes with residual 1.75e-9 and objective agreement 1.39e-12. Tightening
the default method's solver tolerance to 1e-10 also passes; disabling presolve
or scaling does not. All methods retain flows reaching 999999: the successful
methods improve numerical feasibility, not the physiological interpretation of
those flows. The supported settings follow the
[official HiGHS options](https://ergo-code.github.io/HiGHS/dev/options/definitions/).
This is a one-condition result. An earlier setup attempt stopped before solving
because an SBML round trip flattened 64 GPR strings. Its failure remains saved;
the successful diagnostic rebuilds from the original sources to preserve the
full declared signature.

Independent verification passes both saved draft prefixes separately: 113 cases,
226 vectors and 421,954 flux values per attempt. Maximum independently summed
mass-balance error is 8.62e-9, maximum bound violation 1.83e-9, and maximum
between-solver objective difference 5.72e-11. The checker also verifies all
156 structural comparisons, the individual serialized hashes, 166 complete
boundary reports and all five withdrawal controls. Numerical failures remain
failures. Checker development attempts stopped on relative-path handling and
integer/float serialization assumptions; their source and failures are retained
beside the passing verification, rather than overwritten.

## Full curated numerical comparison

The [post-failure comparison definition](../studies/medium-curated-numerics-v1.md)
was pushed as `54b145d` before execution. It applies native GLPK, HiGHS primal
simplex and tighter default simplex to every one of the 43 curated conditions,
with no adaptive choice inside the run. The clean manifest freezes 204 inputs,
fingerprint `80d563f29315ee764d39abf24c76582e0500912e323c6165a2fdaa157bbb354b`.
Its acceptance check retains the original 1e-8 bounds/balance gate and adds a
compensated-sum check at the same tolerance. All reported optimal objectives
must agree, including those whose vectors fail feasibility checks. Rejected
vectors are written before any acceptance or post-solve integrity check.

The [complete results](../../results/medium_preparation_2026_09_06/curated_numerical_comparison/summary.json)
contain all 43 conditions and 129 method records. Every method reports an optimum
in the same 36 conditions and infeasibility in the same seven. All 36 optimal
objective triples agree under the declared tolerance; the largest absolute
difference is 7.92e-9. Numerical feasibility is less consistent:

| Method | Reported optimal | Accepted optimal vectors | Rejected optimal vectors | Reported infeasible |
|---|---:|---:|---:|---:|
| Native GLPK | 36 | 36 | 0 | 7 |
| HiGHS primal simplex | 36 | 26 | 10 | 7 |
| HiGHS default simplex, tighter tolerances | 36 | 9 | 27 | 7 |

There is accepted GLPK-plus-HiGHS evidence in 28 conditions; all three methods
pass in seven. Eight reported-optimal conditions remain without that
cross-solver certificate: L-histidine, vanillin, hydroxy-L-proline,
4-hydroxybenzoic acid, hexanoic acid and oleic acid in MOPS; glucose and
potassium acetate in RCH2. GLPK passes all eight, while both HiGHS methods fail
at least one residual gate. In L-histidine, tighter HiGHS passes the sparse
gate (8.46e-9) but fails the compensated sum (1.09e-8).

Tighter solver tolerances alone therefore do not provide a reliable remedy
across this panel. The largest absolute optimal flow is 231.889 for GLPK,
while both HiGHS methods still reach the inherited 999999 bounds. Maximum
compensated residuals are 9.64e-11, 2.11e-7 and 2.73e-7, respectively. These
observations support continued numerical investigation without changing the
biochemical formulation to force acceptance.

The seven unanimous infeasible reports concern sodium butyrate, 1,2-propanediol,
3-methyl-2-oxobutanoate, 4-methyl-2-oxovalerate, butanol, 3-methyl-1-butanol and
3-methyl-2-oxopentanoate in MOPS. They are retained as solver-reported
infeasibility, with no independent ray proof or fabricated zero-growth value.

The [independent numerical checker](../../results/medium_preparation_2026_09_06/curated_numerical_comparison/independent_verification/result.json)
passes all 43 cases and 129 records, including rejected results. It reconstructs
337,525 flux values and 248,055 balance rows in 115 complete finite vectors,
verifies 172 raw files, 16,555 boundary records and all 204 frozen inputs, and
recovers every acceptance decision above. None of the saved vectors associated
with infeasible statuses is independently feasible. This verification confirms
the classifications and retained limitations; it does not turn rejected vectors
into accepted solutions or prove biological validity.

## Metadata isolation found in review

Independent review found that COBRA's `Model.copy()` does not detach all nested
annotations and notes. V1/v2 preserve the caller's mathematical model during
preparation, but a later metadata edit in their returned copy can affect the
original. This is a real limitation of their claimed isolation contract; their
frozen source and results remain preserved. The additive v3 API deep-copies the
caller model before delegating the unchanged preparation to v2. This also
isolates the caller's active context while retaining correct group ownership
in the prepared copy. New code should use v3. This revision does not
change compartment policy, nutrient rules or model equations.

The [independent v3 comparison](../../results/medium_preparation_2026_09_06/metadata_isolation_comparison_attempt02/report.json)
confirms exact v2/v3 algebra, gene rules, chemistry and report equality across
all 156 conditions and all five withdrawal controls. It checks 332 versioned
reports against the independent boundary oracle without optimization. The
35 new isolation tests cover metadata edits in both directions, groups,
compartment descriptions, active contexts, failures and solver independence.
The maintained suite passes 447 tests plus six subtests at this checkpoint.
Test discovery now explicitly selects `tests/`, because archived source
snapshots are evidence and must not be collected as duplicate test modules.

After adding the numerical comparison, its adversarial acceptance tests and the
independent checker tests, the final maintained suite passes **500 tests plus
six subtests**. The 18 existing scikit-learn warnings remain unchanged.

## Reproduction and next use

Use a clean committed checkout, fresh directories and `requirements-audit.txt`:

```sh
.venv/bin/python -m pytest -q
.venv/bin/python scripts/run_medium_preparation_v2.py --out <fresh-run>
.venv/bin/python scripts/audit_medium_structures.py --run results/medium_preparation_2026_09_06/runs/main_v2 --out <fresh-structural-audit>
.venv/bin/python scripts/diagnose_medium_numerics.py --out <fresh-numerical-diagnosis>
.venv/bin/python scripts/compare_curated_medium_solvers.py --out <fresh-curated-comparison>
.venv/bin/python scripts/verify_curated_medium_numerics.py --run results/medium_preparation_2026_09_06/curated_numerical_comparison --out <fresh-curated-verification>
```

The primary v2 command reproduces the declared attempted panel, including its
numerical gate; it is not advertised as a completed run.

The original helpers remain available to reproduce frozen studies. Subsequent
protocols should explicitly select the new preparation API and declare their
missing-component, compartment and completion policies before evaluation.
Broader namespace/orientation support requires an explicit adapter and tests;
this implementation refuses unsupported structures rather than guessing.

This work improves the reliability of the research infrastructure. It does not
resolve the remaining biological precursor-sharing hypothesis, revise the
31 reported Putida disagreements, or create independent validation data.
