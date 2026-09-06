# Explicit medium preparation v1

6 September 2026. Software conformance on exposed development models.

The precursor diagnostic exposed mixed compartment identifiers: inherited
models contain `C_e`/`C_c` alongside `e`/`c`. The old completion helper creates
new extracellular metabolites in `e`. Some completed exchanges then disappear
from COBRA's inferred exchange list. Sequential medium switching can leave
their uptake open. Synthetic examples reproduce this leak; the preceding
studies' fresh preparation and per-condition restoration protect them against
that demonstrated reuse pattern. No historical growth error is inferred from
the structural finding alone.

COBRA documents boundary classification as heuristic. Its annotations can
override classification based on stoichiometry and compartment; consequently,
replacing the helper with assignment to `model.medium` alone would not resolve
the mixed-label issue. [Official boundary API](https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/medium/boundary_types/index.html).

## Opt-in API contract

New `gembench.medium_preparation.prepare_medium` returns an isolated model copy
and report. The historical helper and every preceding frozen input stay intact.
The caller declares external/cytoplasmic compartment identifiers and any aliases.
Alias normalization is restricted to matching `_e`/`_c` identifiers and records
every change. New completion metabolites retain the declared external compartment.

This version supports canonical BiGG `EX_<component>_e` exchanges with exactly
one metabolite and coefficient −1. Unsupported orientation, shape, identity,
compartment or demand/sink annotation raises. Nutrient bounds must be finite,
non-boolean and nonpositive. Missing components raise unless the caller explicitly
requests reporting them. The export cap is explicitly 1000 by default, matching
the historical protocol. Preparation precedes gene deletions; inactive gene flags
raise to avoid reopening a deleted exchange.

Completion retains the inherited exclusions (`pnto__R`, `fol`, `hco3`) and
gene-less inward-only carrier assumption. Existing exchanges are not evidence
of importability. Missing network connections are reported; no new carrier is
added just because an existing exchange cannot import. Carrier name collisions
and known formula/charge mismatches raise. Unknown or placeholder chemistry does
not become validated chemistry by matching another compartment's metadata.

The helper resets every structurally validated environmental exchange, without
using COBRA's exchange heuristic. All boundary bounds are serialized. Internal
demands/sinks remain unchanged and any that can supply material are identified;
this routine does not establish the physiological validity of those sources.
Caller bounds, stoichiometry, metadata, objective and gene flags remain unchanged
on success and failure. No optimization occurs inside the helper.

## Declared comparison

The source inventory records five local model contexts: the exact saved PpnP
Putida parent; reconstructed historical cycle-7 Btheta, MR1 and Smeli parents;
and the existing Putida iJN1463 reference. Keio/iML1515 is absent locally and is
recorded as unavailable. The four unknown-exposure organisms remain quarantined.
There is no numeric-fitness access or scoring in this study.

Use all 156 mapped carbon conditions (43/25/12/33/43). Metadata and fitness
column headers define eligibility; measured fitness values are not read. For
each condition, independently prepare a legacy copy and a new copy from the
same parent and inherited medium-completion set. Require exact equality of
stoichiometry, reaction bounds, objective, gene rules and chemical metadata
immediately before that condition's solve. Only declared compartment identifiers
and descriptive names may differ. Exact algebra and gene-rule equality entails
the same reaction/gene-deletion LPs under matching later restrictions; this is
stronger than observing a few coincident growth values.

Solve the new prepared LP using native GLPK and a separately translated HiGHS
matrix. Preserve full vectors, primal residuals and objectives. Require agreement
within 1e−8 and feasible residuals below 1e−8, with solver feasibility tolerance
1e−9. Mutually infeasible cases remain explicitly infeasible with null objectives,
not zero growth. Other statuses, unsupported constraints, numerical disagreement
or algebra differences preserve a failed attempt for review. The structural gate
is per condition, not a whole-panel gate before the first solve.

For a separate bounds-only stress test, choose the first glucose condition per
model, or the first mapped condition if glucose is absent. Open every newly
completed component at −0.001 in a synthetic medium, then remove those components.
Compare legacy reuse with a fresh legacy target. Require new reuse to equal new
fresh preparation exactly. This is a software test, not a supplied experimental
medium, growth measurement or biological supplement hypothesis.

## Interpretation and provenance

Commit and push this definition, recipe, source inventory and tested helper
before real-model optimization. Freeze all explicitly named runtime/source
inputs; keep failed attempts and source versions. Independent verification must
check saved full vectors against the source model and both preparation outputs,
without assuming that an equal growth objective implies equal feasible spaces.

No model improvement score or biological transport validation is claimed. A
correct medium reset removes a software dependency on history; it cannot turn
concentrations into fluxes or resolve the biological uncertainty of medium
completion. [Official medium documentation](https://cobrapy.readthedocs.io/en/latest/media.html).
Any future migration of the benchmark to this API requires a separately declared
evaluation; this study does not silently alter earlier protocols.
