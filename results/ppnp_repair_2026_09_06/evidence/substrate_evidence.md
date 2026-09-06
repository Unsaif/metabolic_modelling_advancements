# PpnP substrate audit: forward chemistry supported, reverse evidence conflicts

**2026-09-06.** Both proposed forward reactions have direct purified-protein evidence in **E. coli YaiE/PpnP**, whereas the **KT2440 PP_4248/Q88F51** assignment remains based on homology. The original publication contradicts itself about reverse catalysis. A primary forward-only model hypothesis and a separately declared reversible sensitivity can preserve that uncertainty; neither constitutes experimental validation of KT2440 activity. No model was changed or optimized, and no phenotype values were inspected.

## What the original substrate experiments establish

Sévin's Supplementary Figure 8b/f shows adenosine/inosine conversion in the presence of phosphate, with a product ion assigned as ribose 1-phosphate. Assays used purified enzyme at **50 µg/mL**, **10 mM each substrate**, **three replicates**, and heat-inactivated controls; plotted values are ion abundance, not calibrated rates. Figure 8k documents affinity and anion-exchange purification. Supplementary Table 3 identifies the protein as **b0391/YaiE**. [Original figure, PDF page 9](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.4103/MediaObjects/41592_2017_BFnmeth4103_MOESM62_ESM.pdf), [protein table, row 270](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.4103/MediaObjects/41592_2017_BFnmeth4103_MOESM65_ESM.xlsx).

| Candidate chemistry | Evidence boundary |
|---|---|
| Adenosine + phosphate → adenine + ribose 1-phosphate | Direct E. coli assay; provisional KT2440 gene assignment. |
| Inosine + phosphate → hypoxanthine + ribose 1-phosphate | Direct E. coli assay; provisional KT2440 gene assignment. |
| Reverse nucleoside synthesis | Figure 8i/j and its reversible scheme support it; Table 9 row 2 explicitly reports the reverse assays as negative. Unresolved. |
| Free ribose + phosphate → ribose 1-phosphate | Reported negative in Table 9 and crossed out in Figure 8l; no separate time course is displayed. |

Table 9 also reports no nucleoside cleavage when phosphate was omitted. It supports phosphorolysis over hydrolysis under those assay conditions. The conflicting reverse statement must remain visible; it cannot be silently treated as a typographical error. [Original reaction-validation table, row 2](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.4103/MediaObjects/41592_2017_BFnmeth4103_MOESM71_ESM.xlsx).

The experimental amounts and product-ion traces establish a catalytic candidate, not native flux capacity. Mass alone cannot distinguish isomeric sugar phosphates. The authentic ribose-1-phosphate reverse assays provide additional support, subject to the directionality conflict. Do not infer kinetic constants, substrate ranking, cellular equilibrium, or full knockout rescue from these data. Other reported nucleosides do not justify expanding this bounded intervention; in particular, the generic ribose-phosphate labeling for thymidine warrants separate deoxy-sugar verification.

## What Wen's structures add

The accessible author manuscript describes crystallography, size-exclusion chromatography and dimer-interface stability experiments; it attributes substrate activity to Sévin rather than presenting new substrate-conversion assays. It reports unsuccessful attempts to capture bound substrates. The deposited **7EYJ** structure maps to E. coli **P0C037** and contains sulfate; **7EYP** maps to **P. aeruginosa Q9I3E3**, not KT2440. These are direct structural observations, not Pseudomonas substrate-specificity measurements. [Author manuscript, v1](https://doi.org/10.22541/au.163768486.60510150/v1), [7EYJ](https://www.rcsb.org/structure/7EYJ), [7EYP](https://www.rcsb.org/structure/7EYP).

The manuscript supports a homodimer and describes the absence of a modeled catalytic metal and corresponding metal-binding residues. Those observations do not establish activity after metal depletion. Purified YaiE activity requires no accessory enzyme named in the inspected assay descriptions, but its complete buffer/metal-omission protocol was unavailable. A homodimer does not imply two distinct GPR genes. Do not add ATP consumption, a redox donor, an invented partner gene, or free-ribose phosphorylation on the strength of the family name. The published Wiley full text was inaccessible; the audit distinguishes the full author preprint from the [published abstract](https://pubmed.ncbi.nlm.nih.gov/35094440/).

## Consequences for a controlled model experiment

Retain both forward substrate hypotheses independently of benchmark agreement. Treat PP_4248 as the same provisional gene for both candidate activities only after the separate sequence identity check. Record that opening reverse flux is a sensitivity assumption motivated by inconsistent primary evidence. Do not replace a phosphate-dependent nucleoside-cleavage reaction with a ribose kinase or a free-ribose outlet.

The proposed route also depends on upstream adenosine/inosine supply and downstream ribose-phosphate handling. Evidence for PpnP alone does not validate AhcY, ADA, PPM, their gene assignments, or their physiological capacity. Those dependencies should remain explicit in the structural interpretation and targeted model controls.

The [machine-readable evidence ledger](substrate_evidence.json) records exact figure/table locations, short conflicting extracts, evidence grades, source URLs and SHA-256 fingerprints. Full copyrighted read copies remain outside the repository. This bounded audit does not establish that no further evidence exists.
