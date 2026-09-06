# Hypothetical quinone precursor supply diagnostic

The latest provisional PpnP model retains a PP_5317/chorismate-lyase dependency
that disagrees with 31 exposed pooled-fitness conditions. The gene identity and
UbiC assignment are supported. Current evidence does not justify a replacement
gene rule or an unconstrained alternative reaction. Pooled relative fitness and
isolated complete-deletion growth are different experimental targets.

This diagnostic quantifies the consequences of an explicitly hypothetical external
4-hydroxybenzoate (4HBZ) supply. It does not establish that sharing, carryover or
an endogenous alternative occurred, and will not rescore or patch the model.

## Exact model premise

The parent is the exact saved `both_forward` PpnP arm. All old gene rules,
stoichiometry, biomass coefficients and media definitions remain unchanged.
Source extraction in `scripts/map_quinone_precursor.py` verifies the identity

`CHRPL − EX_4hbz_e − EX_T4hcinnm_e = ε Growth + 4HBHYOX + sink_2ohph_c`,

where ε is the stored quinone demand, 9.68418998275772e-5 mmol/gDW. The
non-growth terms on the right have nonnegative lower bounds. Therefore a
CHRPL-deficient model without coumarate uptake needs net 4HBZ uptake of at least
ε times growth. That necessary bound alone does not prove sufficiency. The
pre-existing 4HBZtex and UHBZ1t_pp transport reactions permit both directions;
permitted secretion is distinct from secretion demonstrated in living cells.

## Declared calculations

The machine-readable definition is `data/studies/precursor_supply_v1.json`.
The runner reads experiment metadata and gene identities, not numeric fitness
values. All 43 mapped carbon-source conditions receive fresh wild-type and
PP_5317-deletion growth calculations under the inherited protocol.

For a growing wild type, fix mutant growth successively at 10%, 50% and 95% of
the parental optimum. Permit 4HBZ uptake up to 0.001 mmol/gDW/h, disallow net
4HBZ secretion, and maximize its signed exchange flux toward zero. Its negative
is the minimum required uptake. Skip this calculation if the parent does not
grow or if 4HBZ itself supplies the main carbon: a trace-only cap in the latter
case would confound precursor availability with carbon starvation. Coumarate
conditions remain included as natural precursor-supply controls.

In the declared glucose/MOPS context, measure mutant maximum growth at uptake
caps of 0, 0.1, 0.5, 1 and 2 times ε times parental growth. Separate controls
use the fixed 0.001 cap: WT; PP_5317 deletion; combined PP_5317/PP_1376 deletion;
inner or outer 4HBZ transport closure; combined PP_5317/PP_5318 deletion; and
4HBZ catabolism closure. These distinguish precursor use, represented uptake and
downstream quinone synthesis from use of the supplement as a carbon source.

For a model donor, retain at least 95% of parental glucose growth and maximize
4HBZ secretion with net import closed. Repeat with inner 4HBZ transport closed.
Neither calculation includes a measured donor population, regulation, release
kinetics or recipient competition. A positive optimum is a possible model flux,
not a prediction that donors actually release that amount.

Solve every declared LP with native COBRA/GLPK and independently translated
HiGHS, guarded against unsupported solver constraints. Require optimal, finite
solutions, agreeing objectives and full primal residuals below 1e-8. Save exact
case specifications, medium contexts, complete primal vectors from both solvers,
selected flux witnesses, solver versions,
input hashes and completed-case records. Preserve any failure in a fresh attempt
directory; do not turn failed solves into zero growth or change an arm to make
it pass. Check the exact precursor equation in each witness. The standard five
energy probes are also run; their boundary closure removes the hypothetical
supply, so repeated energy contexts are not independent sharing experiments.

## Interpretation and provenance

The 0.001 cap is a diagnostic capacity, not a measured supply. Fluxes in
mmol/gDW/h must not be reported as micromolar concentrations. A concentration
estimate would additionally require population biomass, time, relative donor
and mutant abundance, dilution and losses. These quantities are unavailable.
No clinical, organism-wide essentiality or independent predictive validation is
claimed. The question arose after exposed development errors; committing and
freezing it before these solves makes choices auditable, not independent.

Source evidence, including unresolved assay details and enzyme specificity, is
under `results/quinone_precursor_2026_09_06/evidence/`. Existing model-source
notices continue to apply. No new reaction or artificial transporter is introduced.
