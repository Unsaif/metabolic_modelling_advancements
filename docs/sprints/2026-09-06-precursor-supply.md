# Quinone precursor supply: a diagnostic, not a model correction

6 September 2026. Exposed Putida development data only.

The remaining 31 PP_5317 disagreements have a precise possible explanation:
the provisional model's mutant grows when it can import a small amount of
4-hydroxybenzoate (4HBZ), the missing quinone precursor. The amount required
matches the model's material balance in all 31 disputed conditions. The relevant
experiments measured mutants together in a pooled library, whereas the benchmark
simulates isolated complete deletions. This difference deserves investigation.
Neither precursor sharing nor another biological explanation is established.
The model, historical medium definitions and benchmark scores remain unchanged.

## Evidence before intervention

The [identity audit](../../results/quinone_precursor_2026_09_06/evidence/pp5317_identity.md)
confirms that PP_5317, NP_747418.1, WP_010955810.1 and Q88C66 refer to the same
185-residue protein. UbiC activity has strong homology support. A published
P. putida enzyme construct has activity and matching primer ends, but its donor
strain and complete sequence are not specified. It is not an exact KT2440
protein assay. The unresolved PP_3784 family annotation does not establish an
alternative source. There is no basis here to change the gene–reaction rule.
[Kitade et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5835730/).

The [fitness-context audit](../../results/quinone_precursor_2026_09_06/evidence/fitness_context.md)
traces the 31 disagreements to 67 exported experiment-level gene scores. Every
score is above the fixed −2 defect threshold; the range is −1.053 to +0.553.
Averaging does not explain the classification. These are normalized pooled
relative-fitness values, not isolated-mutant growth rates. The assays and model
also lack a calibrated conversion between their respective thresholds.
[Wetmore et al., 2015](https://journals.asm.org/doi/10.1128/mbio.00306-15),
[Price et al., 2018](https://doi.org/10.1038/s41586-018-0124-0).

Neither relevant minimal-medium recipe supplies 4HBZ or nitrate. A vitamin named
4-aminobenzoate is a different compound. Rich defined MOPS includes 4HBZ, but no
mapped condition uses that medium. Recovery/wash protocols suggest possible
carryover questions, without establishing sufficient precursor or stored quinone.
Exact generations, mutant insertion coverage, gene-specific t-statistics and
extracellular precursor measurements remain unavailable. A prior GEM comparison
discusses carryover and cross-feeding as possible explanations in another
organism; it does not demonstrate either here.
[Bernstein et al., 2023](https://link.springer.com/article/10.15252/msb.202311566).

## Fixed question and exact balance

The [study definition](../studies/precursor-supply-v1.md), recipe, source audits,
runner and tests were committed and pushed as
[`836e5fe`](https://github.com/Unsaif/metabolic_modelling_advancements/commit/836e5fe2408c362c851fe98cb379ab6d81028f55)
before these calculations. The parent is the exact saved PpnP `both_forward`
model, SHA256 `b053283adf1257675b771439e99bb0dcaf2ee044c89e8df34096888c26d7f812`.
No new precursor pathway or precursor transporter is introduced. Inherited
medium preparation adds seven vitamin/ion exchange–carrier pairs, all disjoint
from the audited pools. This clarifies the frozen definition's broader
“no new reaction” wording: the precursor intervention adds none, while ordinary
medium completion still applies. The diagnostic changes declared exchange
bounds, deletion controls, growth constraints and objectives temporarily.

Adding the existing precursor, coumarate and quinone-pathway metabolite balances
gives

`CHRPL − EX_4hbz_e − EX_T4hcinnm_e = ε Growth + 4HBHYOX + sink_2ohph_c`,

where `ε = 9.68418998275772e−5 mmol/gDW` is the inherited quinone biomass demand.
Both disposal terms have nonnegative flux. With CHRPL disabled and coumarate
import closed, net 4HBZ uptake must therefore be at least `ε × growth`. This
necessary bound comes from the source model, rather than the observed fitness.
The calculation then asks whether the rest of the network can attain it.

## Results

The [main run](../../results/quinone_precursor_2026_09_06/runs/main/summary.json)
completed 199 cases, each with GLPK and HiGHS: 398 optimization runs in 77.45 s.
Fresh WT calculations grow in 34 of 43 conditions. Unsupplemented PP_5317
deletion still blocks 31; the three natural model rescues use 4HBZ or coumarate
as the already supplied carbon substrate.

All 31 disputed conditions attain 10%, 50% and 95% of their parental model
growth with hypothetical 4HBZ uptake. Across these 93 targets, the minimum
equals `ε × target growth` within 1.37e−15 mmol/gDW/h. Required fluxes range
from 2.40e−6 to 3.01e−4 mmol/gDW/h. The two coumarate conditions require no
resolvable additional 4HBZ uptake at all three targets. Nine non-growing WT
conditions are excluded from this calculation. The 4HBZ carbon-source condition
is also excluded because restricting its uptake would remove its main carbon
supply as well as the precursor.

On glucose/MOPS, the WT optimum is 0.919264794 h⁻¹. The supply needed for
95% of this growth is 8.45721817e−5 mmol/gDW/h. A dose cap expressed as a fraction
of `ε × WT growth` produces the following result:

| Cap multiplier | Mutant growth, h⁻¹ |
|---|---:|
| 0 | Below numerical resolution |
| 0.1 | 0.0919264794 |
| 0.5 | 0.4596323971 |
| 1 | 0.9192647942 |
| 2 | 0.9192829958 |

The slight increase beyond the original WT at the largest cap is permitted
because the supplement also supplies material and can be catabolized. It is
not evidence of a new energy source without input.

At the separately declared uptake cap of 0.001 mmol/gDW/h:

| Glucose/MOPS control | Growth, h⁻¹ |
|---|---:|
| Supplemented WT | 0.9193463964 |
| Supplemented PP_5317 deletion | 0.9193463964 |
| Also delete PP_1376/PcaK | Below numerical resolution |
| Also close inner 4HBZ transport | Below numerical resolution |
| Also close outer 4HBZ transport | Below numerical resolution |
| Also delete PP_5318/UbiA | Below numerical resolution |
| Close 4HBZ catabolism in the supplemented mutant | 0.9192761291 |

Rescue therefore depends on represented uptake and downstream quinone synthesis,
and survives closure of 4HBZ catabolism. Imported precursor can spare biosynthetic
cost even when its catabolism is closed, so exact equality to unsupplemented WT
growth is not required.

The donor optimization permits at most 0.360995851 mmol/gDW/h of 4HBZ secretion
while retaining 95% of parental glucose growth. Closing inner transport removes
resolvable secretion. This is an optimized capacity under existing reversible
transport assumptions. It is not a measured secretion rate or a prediction of
what wild-type cells choose to release. It supplies no donor/recipient abundance,
uptake kinetics, regulation or extracellular concentration.

## Verification and limits

The manifest froze 118 input files from a clean checkout at `836e5fe`; fingerprint
`769a9a4cb904ebbaac48caf4a15b43d0bef78f170c519a3b21c22a546f8813e4`.
All inputs remained unchanged. Complete primal flux vectors from both solvers,
43 medium contexts and every case specification are saved. The largest difference
between solver objectives is 5.51e−12; the largest reported mass-balance residual
is 2.93e−10 and bound violation 1.03e−11, below the declared 1e−8 tolerance.
Values at or below that tolerance are numerically unresolved, rather than
biological measurements of zero. The five closed-boundary energy probes return
zero in each context; their boundary closure also removes the hypothetical
supply, so these are not independent tests of metabolite sharing.

An independent XML/Fraction audit reproduces six exact metabolite-row sums
without importing the model library or runner. The largest saved residual for
those rows is 1.01e−11. It also checks that the seven inherited medium-completion
exchange/carrier pairs do not touch these pools. Its initial conservative guard
rejected a nonempty completion list; that checker source and failed check are
retained. The revised checker reconstructs the additions explicitly. This is
checker development, not a failed or altered primary experiment. Since ε is
small, a precursor-row tolerance of 1e−8 alone corresponds to approximately
1.03e−4 h⁻¹ in the derived growth bound; full-vector and growth checks also apply.

The [independent full-vector verifier](../../results/quinone_precursor_2026_09_06/runs/main/independent_verification_v2/report.json)
reconstructs all 199 cases and 398 vectors from the source model, metadata,
medium tables and deletion rules without solving again. It confirms all input
hashes, reaction bounds, mass balances, objectives and reported case totals.
Its largest recomputed mass-balance residual is 2.927e−10. The
[exact-balance audit](../../results/quinone_precursor_2026_09_06/runs/main/analysis/precursor_balance.md)
provides the separate XML-only derivation. These are additive checks of saved
outputs, not independent biological validation or saved dual optimality proofs.

The full verifier also caught a reporting limitation in inherited preparation:
five new extracellular metabolites use compartment `e`, while the parent's
external compartment is `C_e`. Consequently, COBRA's recognized exchange list
and the saved `exchange_bounds` context omit `EX_4abz_e`, `EX_btn_e`, `EX_nac_e`,
`EX_ribflv_e` and `EX_thm_e`. The reactions remain present, with their actual
medium bounds applied. Verification explicitly reconstructs and checks every
added boundary and carrier. Energy probes use `reaction.boundary`, so they close
these reactions too. The first failed checker attempt and revised source are
preserved under `independent_verification/` and `independent_verification_v2/`.
The experiment is unchanged. A future preparation revision should preserve
compartment identifiers and record all boundary bounds explicitly, followed by
a newly frozen evaluation.

The full test suite passes 287 tests plus six subtests, including 33 independent
checker tests covering malformed or corrupted artifacts and distinct valid
solver solutions. Eighteen existing single-label metric warnings remain.
The [software validation record](../../results/quinone_precursor_2026_09_06/runs/main/analysis/software_validation.json)
preserves the test output and source hashes.

Flux is **not concentration**: mmol/gDW/h cannot be converted to micromolar
without biomass, elapsed time, population composition, dilution and losses.
The inherited zero ATP-maintenance lower bound and provisional PpnP/quinone
assignments remain limitations. No benchmark rescore was performed, and the
31 disagreements remain reported.

## Next decision

This result supports investigating whether the experimental target is being
represented faithfully before adding another biosynthetic pathway. Useful
discriminating evidence would include gene-specific insertion coverage and
uncertainty, isolated-mutant growth and complementation, and measured precursor
availability. These distinguish incomplete disruption, transient stores,
external supply and an endogenous alternative. A source-supported alternative
enzyme remains a separate question. The present computation cannot choose among
these biological explanations.

For future evaluations, record pooled versus isolated assay design, inoculum
history and supplement evidence explicitly. Any numerical precision or medium
policy must be declared before independent evaluation; exposed errors must not
be silently removed by assigning hypothetical nutrients.

## Reproduction

With `requirements-audit.txt` installed, use a fresh output directory:

```sh
.venv/bin/python scripts/run_precursor_supply.py --out <fresh-directory>
.venv/bin/python scripts/verify_precursor_supply.py --study <completed-run> --out <fresh-verification-directory>
.venv/bin/python scripts/audit_precursor_balance.py --run <completed-run> --out <fresh-audit-basename>
```

The runner refuses an existing directory and preserves failed attempts. Source
audit reproduction commands and input hashes are in the evidence documents.
Existing source-model distribution notices continue to apply.
