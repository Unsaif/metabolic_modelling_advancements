# Physical mechanism interpretation of the provisional PpnP intervention

Date: 2026-09-06. These are interpretations of the seven saved physical arms from the frozen development proposal. This audit read their model equations, gene mappings, physical results and selected flux witnesses. It performed **no optimization, phenotype-value inspection or model amendment**. The [companion JSON](mechanism_interpretation.json) preserves raw values, reaction equations, exact row checks and hashes of 33 inputs. Biological evidence and its limits are in [evidence note 11](../../../docs/evidence/11-putida-ppnp-candidate.md).

The forward PpnP intervention removes the model's ribokinase requirement in the declared glucose condition. Removing the added chemistry restores that requirement. The result has a clear algebraic explanation, while the biological identity and capacity of the proposed routes remain provisional.

## Saved conditional growth predictions

All arms use glucose with MOPS minimal medium, glucose uptake capped at 10 mmol/gDW/h, GLPK feasibility tolerance 1e-9 and unchanged quinone biomass demand. ATP maintenance has lower bound **zero**. Growth below is rounded; “0” denotes a raw magnitude below 1e-8, not solver failure. Every classified zero is also below the prespecified growth threshold of 0.001. The JSON retains small signed numerical values.

| Arm | WT | RBK deletion | RBK + PpnP deletion | RBK + AhcY deletion | RBK + ADA deletion | RBK + PPM deletion |
|---|---:|---:|---:|---:|---:|---:|
| Parent | 0.919263665 | 0 | Unrepresented PpnP | 0 | 0 | 0 |
| Adenosine forward | 0.919264794 | 0.919264794 | 0 | 0 | 0.919264794 | 0.919263665 |
| Inosine forward | 0.919263665 | 0.919263665 | 0 | 0 | 0 | 0.919262537 |
| Both forward | 0.919264794 | 0.919264794 | 0 | 0 | 0.919264794 | 0.919263665 |
| Parent, RHCYS closed | 0.919262537 | 0 | Unrepresented PpnP | 0 | 0 | 0 |
| Both forward, RHCYS closed | 0.919264794 | 0.919264794 | 0 | 0 | 0.919264794 | 0.919263665 |
| Both reversible | 0.919264794 | 0.919264794 | 0 | 0 | 0.919264794 | 0.919263665 |

Closing the added reactions while deleting RBK returns numerical zero in every PpnP arm. Deleting inosine hydrolase alongside RBK preserves growth in every PpnP arm. PpnP and RBK single deletions each preserve growth in the principal candidate arm, while the double deletion does not. That is a conditional model prediction of an interaction; it is not experimentally established synthetic lethality.

The parent has no model representation of PP_4248. Its absent PpnP deletion result cannot be interpreted as biological dispensability. PP_5317/CHRPL deletion still blocks growth in all seven glucose contexts; this intervention does not resolve that separate upstream dependence.

## Why the main rescue and its rollback follow from the balances

Let `g` be Growth flux and `ε = 9.68418998275772e-5` the unchanged quinone coefficient. The four internal tail metabolites each occur in exactly two opposing tail reactions. The quinone row sum has only DMQMT production and Growth consumption. Therefore, without an additional diagnostic demand:

`v_OHPHM = v_OMPHHX = v_OMBZLM = v_OMMBLHX = v_DMQMT = εg`.

Each of the three methylation steps produces one SAH. Independently summing the saved SBML rows for SAH, S-ribosylhomocysteine, adenosine, inosine and cytosolic ribose yields:

`v_RBK + v_PUNP1 + v_PUNP5 ≥ 3εg`.

These are signed fluxes; absent PUNP reactions have zero flux. The inequality follows because all 14 positive terms in this row sum have nonnegative flux bounds, and three are the obligatory methylations. RBK and the added PUNP reactions are the only negative terms. Blocking all three therefore forbids positive exact steady-state growth. Closing RHCYS does not remove this constraint: its weighted coefficient is zero.

Two smaller row sums explain the other controls:

- SAH + S-ribosylhomocysteine + ribose gives `v_RBK + v_AHCi ≥ 3εg`.
- Adding adenosine to that sum gives `v_RBK + v_ADA + v_PUNP1 ≥ 3εg`.

Thus AhcY is required for the observed rescue when RBK is absent. ADA is also required when only the inosine phosphorolysis reaction is supplied. Adenosine phosphorolysis bypasses ADA. The JSON records these exact rational row sums and verifies their signs and bounds for all seven saved models. These are necessary model constraints, not proofs of enzyme function in cells.

In the principal arm's saved RBK-deletion witness, growth is 0.919264794239 and each terminal methylation carries about 8.90233e-5. AhcY, adenosine phosphorolysis and PPM each carry about 2.67070e-4, exactly the three-methylation disposal amount within rounding. Their combined equation is:

`SAH + H2O + Pi → homocysteine + adenine + ribose-5-P`.

The inosine-only RBK-deletion witness instead uses AhcY → ADA → inosine phosphorolysis → PPM. Its WT witness still uses the old ribokinase branch, illustrating why an unused reaction in one optimum does not show that the reaction lacks a feasible role. Across growing saved witnesses, the checked tail, demand-coupling and disposal equalities differ by at most 4.75e-17. These arithmetic checks concern selected flux relationships; separate numerical verification addresses the full optimization results.

## PPM deletion exposes an inherited alternative

PPM deletion does **not** abolish the rescue. In the principal RBK/PPM deletion witness, PPM is zero and three older reactions carry equal flux magnitude, about 2.67070e-4:

| Reaction and direction | Saved model equation | Gene representation |
|---|---|---|
| PNP, reverse | nicotinamide + ribose-1-P + H+ → nicotinamide riboside + Pi | No GPR |
| RNMK, forward | nicotinamide riboside + ATP → NMN + ADP + H+ | PP_4218 |
| NMNN, forward | NMN + H2O → nicotinamide + ribose-5-P + H+ | No GPR |

Nicotinamide, its riboside and NMN cancel. The exact net reaction is:

`ribose-1-P + ATP + H2O → ribose-5-P + ADP + Pi + H+`.

All three directions are permitted by existing bounds. This is an ATP-consuming alternative to the mutase in the model, with nicotinamide recycled internally. It explains the selected witness without requiring ribose-1-P export or a new artificial source. It also introduces another biological uncertainty: two reactions have no gene assignment, and this audit has not established their KT2440 activity. Persistence of rescue after PP_1777 deletion therefore validates neither PP_1777's disputed substrate specificity nor the alternative route. The tiny optimized growth differences are not measured enzyme capacities or evidence of a physiological preference.

## Closure and reversal sensitivities

Closing RHCYS alone leaves WT growth through the alternative AhcY → ADA → inosine hydrolase → RBK chain. RBK dependence remains. With PpnP added, rescue persists when RHCYS is closed. The candidate rescue consequently does not require the questionable historical RHCYS chemistry, although that chemistry remains part of the open parent arm. AhcY single deletion becomes growth-blocking in both RHCYS-closed contexts.

Allowing reverse PpnP chemistry produces the same displayed glucose growth outcomes and a forward PUNP1 witness. Reverse activity is unnecessary for this observed rescue. These results do not resolve the conflicting primary reverse-assay evidence, establish equivalence in other media or exclude other optimal flux distributions. All seven saved tests report zero generation for the five tested energy currencies with nutrient uptake closed; this limited check is not a general thermodynamic validation.

## Most discriminating biological follow-up

1. **Resolve PpnP independently of everything downstream.** Assay the sequence-verified KT2440 PP_4248 protein with adenosine and inosine separately, identify ribose-1-P and the released base, and include phosphate omission, no-enzyme and free-ribose controls. Test reverse synthesis separately. This addresses the candidate enzyme itself.
2. **Resolve the PPM assignment separately.** Assay PP_1777 with ribose-1-P and its annotated mannose-phosphate substrates, distinguish product isomers and measure rates. A successful PP_4248 assay does not establish this downstream activity or sufficient capacity.
3. **Distinguish the two downstream routes.** In cell extracts, trace ribose-1-P conversion while separately restricting PPM and RNMK activity and measuring nicotinamide intermediates. Establish the identities of the gene-less PNP/NMNN activities before assigning genetic controls. An isolated PP_1777 deletion would be inconclusive in the presence of an alternative route; supplied ribose-1-P should not be assumed to enter intact cells.
4. **Then test physiological contribution.** Compare matched WT, rbk, ppnP and combined perturbations with complementation and targeted SAH/nucleoside/ribose-phosphate measurements in the declared condition. This separates plausible chemistry from sufficient, expressed activity in cells.

These are proposals only. No experiments or outreach were performed. Quinone demand, upstream quinone chemistry, enzyme capacity, maintenance assumptions and species-specific functional assignments remain unresolved. The physical rescue is a well-explained model result; its biological validity and predictive value require separate evidence.
