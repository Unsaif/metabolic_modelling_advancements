# PP_5317 fitness measurements and experimental context

This is a post-outcome audit of exposed Putida development data. It adds no model reactions, runs no optimization and inspects no other organism’s numeric phenotypes. The companion JSON preserves all 43 condition groups, 92 experiment-level scores, original metadata and input hashes.

**The 31 remaining disagreements are not caused by averaging away a strongly defective replicate.** Their 67 exported PP_5317 scores range from −1.053 to +0.553; every score is above the fixed −2 defect threshold. Condition means range from −0.842667 to +0.3625. These data establish an absence of a severe pooled competitive deficit at this threshold. They do not establish monoculture viability of a complete deletion.

## What the assay measures

RB-TnSeq compares mutant barcode abundance before and after competitive growth in a pooled library. Gene scores combine insertion strains and are normalized around the typical gene. A score near zero is therefore relative performance within the pool, rather than an absolute growth rate or direct viability measurement. Wetmore’s method paper describes typically 4–6 generations. [Wetmore et al., 2015](https://journals.asm.org/doi/10.1128/mbio.00306-15).

Price’s later methods describe typically 4–8 population doublings to saturation, central 10–90% insertions and abundance filtering; replicates need not have identical concentrations. Exact generations for these Putida samples are absent from the local metadata. The benchmark’s growth threshold 0.001 and fitness threshold −2 compare different quantities and do not calibrate one into the other. [Price et al., 2018](https://doi.org/10.1038/s41586-018-0124-0).

For the Putida fatty-acid/alcohol work, JBEI-1 was recovered in LB+kanamycin to OD600 0.5, sampled for time zero, washed **once** in carbon-free MOPS and diluted 1:50 into 10 mM carbon source; cultures used 10 mL tubes at 30 °C and 200 rpm. The separate individual-growth protocol used **two** washes and 1:100. Those assays must not be conflated. This is consistent with set15/set16 metadata, but is not a unique expName linkage. [Thompson et al., 2020](https://journals.asm.org/doi/10.1128/aem.01665-20).

The aromatic-substrate paper instead describes JBEI-1 LB/kanamycin recovery, a MOPS wash and 1:50 transfer into 10 mM substrates, 600 µL deep wells at 30 °C/700 rpm, and combining two 600 µL samples before BarSeq. This fits set12’s substrate/format pattern; no explicit expName link, exact wash count, endpoint duration or generation count is supplied. Engineered production strains in that paper belong to separate assays. [Incha et al., 2020, §4.5](https://doi.org/10.1016/j.mec.2019.e00119).

A primary GEM/fitness comparison identifies cofactor carryover and cross-feeding as possible explanations for discrepancies of this kind, and warns that unsupported biosynthetic additions can create errors elsewhere. That work concerns a different organism; it supplies a methodological caution, **not evidence that PP_5317 is cross-fed in KT2440**. [Bernstein et al., 2023](https://link.springer.com/article/10.15252/msb.202311566).

## Mapping, replication and uncertainty

The raw locus and sysName are both PP_5317, described as “Probable chorismate pyruvate-lyase”. The frozen model uniquely maps it to NP_747418_1, assigned only to CHRPL. The saved benchmark axis and every condition mean match an independent reconstruction from the raw export to within 1e-12. These exported gene-level scores already combine insertion strains; they are not raw barcode counts.

The 43 mapped groups contain 34 model-WT-growing conditions. PP_5317 deletion blocks 31 and grows in the three conditions below. Nine groups with model WT growth below 0.001 are outside this gene-level comparison; their 19 scores are retained in JSON for coverage accounting. Fourteen additional eligible metadata groups are unmapped.

The 67 disagreement samples include 52 Putida_ML5_JBEI and 15 Putida_ML5 labels; 20 include DMSO. Every sample is annotated aerobic and 30 °C. Some averages combine library labels, DMSO presence, dates, vessels or concentrations. For example, RCH2 glucose combines 20/40 mM and tube/plate experiments; RCH2 acetate combines 5/20 mM. MOPS benzoate mixes 5/10 mM and both library labels. Hence “replicate average” is broader than a matched repeat under every experimental detail. This heterogeneity does not explain away the −2 classification because every contributing score remains above it.

All 67 disagreement samples pass a descriptive comparison to the historical experiment-wide quality criteria. Their median-gene counts range 68–368, half-gene correlation 0.1083–0.2979 and half-gene median absolute difference 0.1490–0.3351. These metrics do not give PP_5317-specific read depth or statistical uncertainty. Local t-scores, insertion positions/counts and raw barcode counts are unavailable; the public t-table endpoint returned HTTP 403 on 6 September 2026. The specific-phenotype export has no PP_5317 row, which is not evidence that its fitness is statistically equivalent to zero.

## Three natural model rescues

| Condition | Exported PP_5317 fitness values | Mean | Saved model knockout growth |
|---|---|---:|---:|
| 4-Hydroxybenzoic Acid / MOPS minimal media_noCarbon | +0.156, -0.223 | -0.0335 | 0.7157983 |
| p-Coumaric acid / MOPS minimal media_noCarbon | +0.365, +0.843 | +0.6040 | 1.0074059 |
| p-Coumaric acid / RCH2_defined_noCarbon | +0.285, +0.057 | +0.1710 | 1.0079711 |

These six measurements use DMSO and are all above −2. The existing model can take up 4-hydroxybenzoate or form it during coumarate catabolism, bypassing CHRPL. “Rescue” here describes the saved model prediction under an already supplied substrate, not a dedicated experimental PP_5317 complementation or dose-response assay.

## Medium, uptake and nitrate confounds

The retained FEBA source recipes for MOPS minimal and RCH2_defined_noCarbon contain **no 4-hydroxybenzoate and no nitrate**. RCH2 contains Wolfe vitamins including **4-aminobenzoate**, a chemically distinct compound. MOPS Rich Defined lists 4-hydroxybenzoate at 0.01 mM, but none of the 43 mapped groups uses that medium, and the primary recovery protocols above use LB. This audit found no evidence that Rich Defined was an inoculum or that residual LB supplied enough precursor.

The only nitrate-labelled Putida rows are set27IT016/set27IT017, classified as nitrogen-source experiments in MOPS minimal media_Glucose_noNitrogen with 10 mM sodium nitrate. They are outside this carbon benchmark. A reference to the soil isolate RCH2 or to its nitrate-reduction studies must not be imported as a condition of the KT2440 assay. All selected metadata is aerobic; dissolved oxygen histories are not available.

Existing 4HBZtex and UHBZ1t_pp are reversible in the saved model. Its glucose physical context closes 4HBZ and nitrate uptake. The carbon-condition mapping separately opens 4HBZ or coumarate uptake at −10 mmol/gDW/h when supplied. That model bound is an assumption and is not the same physical quantity as an experimental concentration in mM. A hypothetical minimal-uptake or donor-secretion calculation could test model sufficiency without adding a transporter; it could not demonstrate measured secretion, cross-feeding or carryover.

## Evidence that would distinguish the explanations

Obtain PP_5317 barcode-level insertion coverage and t-statistics, then compare a sequence-verified clean deletion with its complemented strain in separately cultured, washed, defined medium across serial transfers. Test 4HBZ addition and, separately, cell-free donor-conditioned medium while measuring 4HBZ and quinone pools. Quantify inoculum and endpoint population doublings. These controls would distinguish incomplete disruption, transient stores, extracellular precursor supply and an unrepresented biosynthetic route; none is established by the current scores.

The 31 disagreements remain valid benchmark errors under its fixed scoring rule. They should not be removed from the report, and they do not by themselves justify changing CHRPL’s GPR, adding a precursor source, or claiming a new biological pathway.

Reproduce with `python scripts/extract_quinone_precursor_fitness_context.py --out <fresh-directory>`. The JSON records 15 input hashes, including the extraction script. No old inputs are overwritten.
