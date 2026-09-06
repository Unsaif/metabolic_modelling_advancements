# Post-outcome trace of the two new gene dependencies

This analysis was prompted by the completed development benchmark: the parent analysis reported 65 adverse binary changes in each pathway-plus-quinone-demand arm, involving PP_2458 in 34 conditions and PP_5317 in 31. It explains constraints in the saved model; it does not validate the repair or select changes to recover agreement. No optimization, model modification, new fitness calculation, or biological literature search was performed. All observations concern the already exposed Putida development organism.

The companion [JSON](postrun_dependency_trace.json) records input hashes, exact row sums, complete relevant metabolite adjacency, reaction directions, gene mappings and condition metadata. The model inspected is [curated_path_template/model.xml.gz](../runs/main/curated_path_template/model.xml.gz). Its exchange bounds precede benchmark medium application; original open exchanges must not be interpreted as assay supplements.

## Ribokinase is the only exit from a byproduct ledger

The gene annotation calls **PP_2458** ribokinase. The saved mapping connects it to **NP_744606_1**, whose sole associated reaction is **RBK**:

`ATP + D-ribose → ADP + H+ + D-ribose 5-phosphate`

Each new quinone molecule requires three methylations, **OHPHM**, **OMBZLM** and **DMQMT**. Each consumes one S-adenosylmethionine (`amet_c`, SAM) and produces one S-adenosylhomocysteine (`ahcys_c`, SAH). The four internal metabolites connecting the five transferred reactions each occur in exactly two reactions with coefficients +1 and −1. Their balances therefore make all five tail fluxes equal to a common value `t`.

Two modeled routes dispose of SAH:

| Route | Stored reactions | Products relevant to ribose disposal |
|---|---|---|
| Direct cleavage branch | AHCYSNS (PP_3254), then gene-less RHCYS | SAH → S-ribosylhomocysteine → free ribose |
| Adenosine branch | Reversible AHCi (PP_4976), ADA (PP_0591), INSH (PP_2460) | SAH → adenosine → inosine → free ribose |

These arrows summarize the branch; full water, proton and other coproduct stoichiometry is retained in the JSON. Choosing the other represented branch cannot eliminate the free-ribose disposal requirement.

For an exact certificate, sum the five zero-balance rows for:

`ahcys_c + rhcys_c + adn_c + ins_c + rib__D_c`

This is a declared algebraic ledger, not a newly established biological moiety. The internal transformations above cancel. **RBK has coefficient −1 and is the only direction-compatible exit.** Fourteen reactions have positive coefficients, and all their lower bounds are zero: ADNt2, AMMQLT8, AMMQT8_2, DMQMT, GNNUC, HCYSMT, INSt2, NTD11, NTD7, OHPHM, OMBZLM, RIBabcpp, UPP3MT and XTSNH. UPP3MT contributes +2; the others contribute +1. Thus every remaining term is nonnegative.

Independently, the sum of `q8_c` and `q8h2_c` balances contains only `+v_DMQMT − ε v_Growth`, where the template demand coefficient is `ε = 9.68418998275772e−5`. Consequently:

```text
t = v_DMQMT = ε v_Growth
v_RBK ≥ v_OHPHM + v_OMBZLM + v_DMQMT
v_RBK ≥ 3 ε v_Growth
```

Here `3 ε = 0.0002905256994827316`. With positive ε, deleting RBK forbids positive exact steady-state growth in this stored model. This implication follows from the coefficients and reaction directions; it does not depend on choosing an optimized flux distribution. The JSON uses exact rational arithmetic on the stored floating-point coefficients. Numerical solvers still have finite tolerances.

The apparent boundary alternatives do not invalidate the proof. `RIBtex` connects extracellular and periplasmic ribose, while the cytoplasmic transporter `RIBabcpp` is irreversible inward. An extracellular ribose exchange therefore provides no cytoplasmic export route. The retained `sink_2ohph_c` drains a quinone precursor, not SAH, ribose, or quinone, and cannot avoid the methylation burden. Biomass directly consumes a small amount of SAM without producing SAH; that baseline demand alone does not impose this disposal burden.

## Chorismate lyase has two substrate-dependent bypass routes

The gene annotation calls **PP_5317** probable chorismate pyruvate-lyase. It maps to **NP_747418_1**, whose sole associated reaction is **CHRPL**:

`chorismate → 4-hydroxybenzoate + pyruvate`

The quinone pathway uses 4-hydroxybenzoate through **HBZOPT → OPHBDC → OPHHX → OHPHM** and the remaining transferred tail. The complete `4hbz_c` adjacency identifies two ways to obtain this precursor without CHRPL:

1. Import supplied 4-hydroxybenzoate through **EX_4hbz_e → 4HBZtex → UHBZ1t_pp**. The represented transport genes are PP_1121 or PP_4198, and PP_1376, respectively.
2. Import p-coumarate through **EX_T4hcinnm_e → T4HCINNMtex → T4HCINNMtpp**, then use **4CMCOAS (PP_3356) → COCOAHA (PP_3358) → VNDH_2 (PP_3357)** to produce 4-hydroxybenzoate via coumaroyl-CoA and 4-hydroxybenzaldehyde.

The nominal reverse direction of **SUCBZT2** is not another sustained source: `4hbzcoa_c` appears only in this reaction, so its own balance forces net SUCBZT2 flux to zero. **4HBHYOX** consumes 4-hydroxybenzoate and cannot run in reverse under the stored bounds.

The condition metadata contains exactly three supplied-substrate cases corresponding to these routes: 4-hydroxybenzoic acid in MOPS minimal medium, and p-coumaric acid in MOPS minimal and RCH2 defined media. Their presence is consistent with 31 CHRPL changes among 34 growth conditions. This trace verifies the available directions and supplied substrate labels, not the individual benchmark calls or feasibility of these routes in isolation. The independent benchmark verification addresses those calls. Under the saved physical glucose medium, 4-hydroxybenzoate, coumarate and ribose imports are closed; glucose uptake is bounded at 10. ATP maintenance has lower bound zero, so these calculations impose no positive ATP maintenance requirement.

## Mechanistic questions left open

- **Is the SAH salvage chemistry correct for this organism?** In particular, the gene-less RHCYS step and the alternative adenosine branch need species-specific biochemical and gene evidence. A reaction missing from this model is not evidence that the organism lacks its activity.
- **Can an evidenced pathway dispose of the ribosyl products without RBK?** The stored adenosine and inosine neighborhoods have no phosphorolytic outlet, and the free-ribose neighborhood has no cytoplasmic efflux. Nucleoside salvage, alternative ribose metabolism or transport are questions for independent evidence, not instructions to add a sink.
- **Why does pooled fitness disagree with the predicted individual knockout?** Mutant identity, residual gene function, assay timing and the distinction between pooled fitness and individual deletion growth need examination. Carryover or crossfeeding remain untested possibilities; neither follows from reaction adjacency, and either requires a plausible metabolite supply and uptake mechanism.
- **Does actual CHRPL deficiency depend on precursor supply as represented?** Independent evidence from defined precursor availability, quinone measurements and gene complementation could distinguish an endogenous alternative source, uptake-mediated rescue and an incorrect reaction assignment. The model's two supplied-substrate routes give a falsifiable prediction without establishing that they explain the observed biology.

The repair remains provisional. The transferred half-oxygen hydroxylation chemistry, species-specific pathway assignments and physiological quinone dilution amount remain unresolved. Any positive demand creates the algebraic RBK dependency above; selecting a coefficient to improve the benchmark would not resolve the missing biological evidence. No reaction, coefficient, GPR, medium or threshold was changed by this analysis.
