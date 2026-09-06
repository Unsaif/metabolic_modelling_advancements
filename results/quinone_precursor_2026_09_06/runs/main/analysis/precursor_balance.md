# Independent precursor balance and interpretation audit

The saved results satisfy independently derived precursor balances. This audit uses the frozen SBML XML and both saved solver flux vectors; it imports no model library or runner helper and performs no optimization or phenotype analysis.

The exact stored coefficient is ε = 7145670966924535/73786976294838206464 = 9.6841899827577e-05 mmol/gDW.

`CHRPL − EX_4hbz_e − EX_T4hcinnm_e = ε Growth + 4HBHYOX + sink_2ohph_c`.

Additional rows establish `DMQMT = ε Growth`, `HBZOPT = ε Growth + sink`, and `EX_4hbz_e = −UHBZ1t_pp = −4HBZtex`. With CHRPL disabled and no coumarate import, the nonnegative disposal terms require uptake ≥ ε times growth. Transport closure prevents net supply; HBZOPT closure prevents growth independently of the external supply.

Verified 398 saved witnesses across 199 cases and 43 conditions. The maximum absolute residual among these six independently summed rows was 1.002e-11; the maximum difference between saved solver objectives was 5.501e-12. The inherited protocol adds seven exchange/carrier pairs for vitamins and ions. Their independently reconstructed stoichiometry has zero contribution to all six pools; bounds and deletions also preserve the rows.

## Minimum supply and dose caps

Across 93 minimum-supply cases without coumarate import, uptake differed from ε times target growth by at most 1.368e-15 mmol/gDW/h. Required uptake ranged from 2.39926e-06 to 0.000301286. The 6 coumarate cases required at most 0 absolute 4HBZ uptake. The full condition/fraction records are in the JSON.

A cap of d × ε × parental growth imposes mutant growth ≤ d × parental growth. It does not guarantee equality. At larger caps, the supplement can also alter carbon/precursor availability; an exact plateau at unsupplemented growth is not required.

| Dose multiplier | Uptake cap (mmol/gDW/h) | Growth (h⁻¹) | Fraction of unsupplemented WT |
|---:|---:|---:|---:|
| 0 | -0 | 0 | 0 |
| 0.1 | 8.9023349e-06 | 0.09192647943 | 0.1 |
| 0.5 | 4.4511675e-05 | 0.4596323971 | 0.5 |
| 1 | 8.9023349e-05 | 0.9192647942 | 1 |
| 2 | 0.0001780467 | 0.9192829958 | 1.0000198 |

## Controls

| Control | Growth (h⁻¹) | 4HBZ uptake (mmol/gDW/h) |
|---|---:|---:|
| supply_wild_type | 0.9193463964 | 0.001 |
| supply_pp5317 | 0.9193463964 | 0.001 |
| supply_pp5317_pcaK | -7.505767469e-16 | -0 |
| supply_pp5317_inner_transport_closed | 0 | -0 |
| supply_pp5317_outer_transport_closed | 0 | -0 |
| supply_pp5317_ubiA | 0 | -0 |
| supply_pp5317_catabolism_closed | 0.9192761291 | 8.9024447e-05 |

Positive rescue under a declared supply is a feasibility result. Loss of rescue under the transport/downstream controls identifies the represented route; preserved rescue after catabolism closure separates precursor incorporation from the need to degrade the supplement.

## Deliberately maximized donor secretion

| Donor control | Maximum secretion (mmol/gDW/h) | Growth (h⁻¹) | Fraction of WT |
|---|---:|---:|---:|
| donor_wild_type | 0.3609958506 | 0.8733015545 | 0.95 |
| donor_transport_closed | -1.52915505e-16 | 0.9192647942 | 1 |

The balance requires secretion ≤ CHRPL flux − ε times growth. These calculations maximize secretion subject to a growth floor. They do not predict spontaneous secretion, demonstrate release by cells, or show that a donor population supplies enough material to recipients. The model includes no sharing dynamics, biomass fractions, accumulation time or extracellular losses.

## Numerical and biological limits

Absolute values ≤1e-8 are treated as operational numerical zeros, not biological zeros. Because ε is small, a precursor-balance residual allowance of 1e-8 translates to 0.0001033 h⁻¹ in the derived growth inequality; the actual saved row residuals are reported above. The audit checks saved equations and agreement, not a dual optimality certificate.

Fluxes in mmol/gDW/h cannot be converted into a medium concentration without a biomass/time balance. Neither successful rescue nor a positive donor optimum establishes carryover, cross-feeding or an alternative native enzyme. All inherited quinone, PpnP, transport and zero ATP-maintenance assumptions remain. No biological correction or new benchmark score is accepted here.

The preceding experiment froze 118 files with fingerprint `769a9a4cb904ebbaac48caf4a15b43d0bef78f170c519a3b21c22a546f8813e4`. This additive auditor is separate from those experiment inputs. Its own source hash and all directly read artifact hashes are recorded in `precursor_balance.json`; the frozen experiment input hashes were also rechecked.

The first checker attempt conservatively stopped on the nonempty medium-completion inventory before checking fluxes. Its source and failure record were retained. The revised checker explicitly reconstructs that inherited completion and verifies pool disjointness. This refinement does not represent a failed primary run or a change to the experiment.
