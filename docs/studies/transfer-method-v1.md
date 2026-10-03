# Transfer study v1: do fitness-blind corrections of draft models transfer to new organisms?

Declared 3 October 2026 by Claude (Opus 5.5) in Tim Hulshof's project repository. This document, the code and data files named in `data/studies/transfer_v1/method_freeze_plan.json` and the commit that adds them form the **method freeze**. Nothing in it may change after panel A's fitness tables are downloaded; a change after that point starts a new development iteration that needs a new independent panel.

The study follows the [evaluation protocol](evaluation-protocol-v1.md) and the [independent panel declaration](../../data/studies/independent_panel_eligibility_v1.json). It answers the main criticism of the [scientific audit](../reviews/2026-09-06-scientific-audit.md): every earlier improvement was measured on the same four organisms whose fitness data guided it.

## 1. Question and hypotheses

Automatically reconstructed draft models (the 2017 EMBL GEMs collection, built by CarveMe) predict which genes a bacterium needs to grow on a carbon source. Fitness Browser RB-TnSeq experiments measure it. During development (Sprints 3–4, cycles 1–7) we found systematic errors in the drafts and wrote correction rules that use only the genome annotation and biochemistry, never the fitness data of the organism being corrected. The question is whether those rules improve predictions **on organisms that played no part in their development**.

- **H1 (universal rules).** Arm U′ (the draft plus five automatic, annotation-driven transforms) agrees better with fitness data than arm B0 (the draft with the minimal gap-fill only).
- **H2 (blind adjudication).** Arm M (U′ plus gene-rule decisions made by an AI curator from annotation alone, under a fixed written procedure) agrees better than U′.

Neither hypothesis is about whether the models are good in absolute terms; both are paired comparisons on the same genes and conditions.

## 2. Organisms

Panel A of the independent panel ([selection record](../../data/studies/independent_panel_v1/panel_selection.json)): **Cola** (*Echinicola vietnamensis* DSM 17526), **Dino** (*Dinoroseobacter shibae* DFL 12), **Dyella79** (*Dyella japonica* UNC79MFTsu3.2), **MycoTube** (*Mycobacterium tuberculosis* H37Rv), **PS** (*Dechlorosoma suillum* PS) and **RPal_CGA009** (*Rhodopseudomonas palustris* CGA009). Panel B (Caulo, Cup4G11, Marino, Miya, PV4, SB2B) stays untouched for a later method version. All six panel A organisms are evaluated; none is replaced or dropped because of its results.

Development organisms (Btheta, Putida, MR1, Smeli) are reported alongside for context only; their numbers are retrospective.

## 3. Inputs fixed before outcomes

Only non-outcome inputs are used before panel A's fitness tables are downloaded: the pinned EMBL draft (`models/embl_pinned/`), the CarveMe universe (`external/carveme/universe_bacteria.xml.gz`), Fitness Browser gene and experiment metadata (`data/fitness_browser_panel/<org>/genes.tsv`, `experiments.tsv`) and protein sequences, NCBI protein sequences of the draft's genes, and the FEBA medium recipes (`data/studies/transfer_v1/feba/`, bitbucket commit 803273865ea9).

### 3.1 Gene mapping
Model genes are matched to Fitness Browser loci by protein sequence (`scripts/build_panel_gene_maps.py`): exact identity; one sequence contained in the other covering at least 90% of the longer; or equal length with at least 97% identical positions. Ambiguous matches are left unmapped. This step was run before this freeze (commit 9a702ab); 95.6–100% of model proteins map.

### 3.2 Media
Each base medium is mapped by `gembench/feba_media.py` from its FEBA recipe (mixes expanded) through the component dictionary `data/studies/transfer_v1/feba_component_bigg.tsv`: inorganic salts to their ions; vitamins, cofactors, hemin, nucleobases and reductants to their metabolites limited to trace uptake; other defined organics unlimited; buffers and chelators dropped (counter-ions kept); complex ingredients (yeast extract, peptone, tryptone, casamino acids and the like) make the medium **undefined**, and its experiments are excluded. Water, protons, CO2 and the full trace-metal set (capped at 0.1 mmol/gDW/h, decision D17) are always present; O2 is present when the medium's carbon-source experiments are recorded as aerobic in the experiment metadata. Applied to the six development media, the rule reproduces the hand-curated development table component for component; it differs only in two trace limits (cysteine counted as a reductant in MOPS rich medium; methionine unlimited in Varel–Bryant medium).

Panel components missing from the dictionary are added after this freeze, in a separately marked section, by chemical identity only (the compound's own metabolite and its counter-ions, using the class definitions above). No outcome is consulted. No light or photon input is modelled: the CarveMe universe has no anoxygenic photosystem, so phototrophic experiments are simulated as dark anaerobic growth, which the draft models are not expected to support.

### 3.3 Carbon sources
Condition names are mapped to BiGG metabolites by the development table (`data/reference/fitness_browser_carbon_sources_bigg.tsv`) where it already lists the name. New names are mapped after this freeze, in a separately marked section of `data/studies/transfer_v1/carbon_sources_bigg.tsv`: the BiGG metabolite of exactly the named compound (salts and hydrates to the organic compound; a racemate to both enantiomers; an unspecified amino acid to the L-form and an unspecified sugar to its natural D- or L-form, marked medium confidence); polymers, mixtures and compounds without a BiGG identifier map to none. A mapped compound whose exchange reaction is absent from a model is reported per condition (as in development). Experiments with a second condition other than DMSO are excluded (as in development).

### 3.4 Base model B0
The pinned draft receives the minimal gap-fill (`gembench.gapfill`, CarveMe universe, energy-generating-cycle gate, minimum growth 0.05 /h, 900 s limit) for growth on one reference condition chosen from metadata by `gembench.transfer.choose_reference_condition`: the main medium is the base medium with the most carbon-source experiments; within it, the first tested substrate in the order D-glucose, D-fructose, glycerol, L-lactate, pyruvate, succinate, acetate, fumarate, L-malate, 2-oxoglutarate, citrate, L-glutamate; otherwise the condition with most experiments (ties by name). The rule reproduces the hand-chosen references of all four development organisms. Using a tested substrate assumes the organism grows on it, which any RB-TnSeq carbon-source experiment implies; no fitness value is read. A base model is written once and never overwritten.

### 3.5 Arms
Transforms are implemented in `gembench/transfer.py` and defined as ordered lists in `data/studies/transfer_v1/arms.json`.

| Arm | Transforms | Role |
|---|---|---|
| B0 | none | baseline |
| **U′ (= UNQ)** | universal reaction patches (CBMKr catabolic direction; CBPS subunit conjunctions only; MECDPDH4E irreversible) → conjunction-rule normalisation → universal additions (glycolaldehyde diffusion and exchange) → F0F1 ATP synthase if the draft has none and the annotation names at least six of the eight subunit types → menaquinone removed from biomass when the annotation shows fewer than two menaquinone pathway steps | primary H1 |
| **M** | U′ followed by the organism's blind adjudication decisions (`decisions:` transform) | primary H2 |
| U, UN, UQ, Uonly_rxn, Uonly_atp | subsets of U′ | secondary attribution |
| M_R6, M_R1R2 | U′ plus only the R6 decisions, or only the R1/R2 decisions | secondary attribution |

U′ was chosen after all development arms had been run: it has the highest MCC of the development arms in every development organism (tied with UN in MR1 and Btheta and with UQ in Putida), and each of its separately tested parts (the reaction patches, the ATP synthase rule, normalisation, the menaquinone rule) has a non-negative paired difference in every development organism except the reaction patches in Btheta (−0.006, from the CBPS patch: Btheta has a second carbamoyl-phosphate synthase that the rule cannot see). Section 8 lists the numbers.

### 3.6 Blind adjudication (arm M)
For each organism, `scripts/transfer_candidates.py --arm UNQ` builds a packet from the U′ model and metadata only: the gene-less reactions and the OR-rule reactions that are essential in at least one mapped condition where the U′ model grows, each with its equation, annotations, current rule and, for every gene, the Fitness Browser description, RefSeq definition and genomic neighbourhood; plus the full gene table and the procedure `docs/studies/transfer-adjudication-procedure-v1.md`. The packet is written outside the repository. One fresh subagent per organism (the same Claude model as the study session, no web access, instructed to open only its packet; prompt `docs/studies/transfer-adjudication-prompt-v1.txt`) writes `decisions.json`. Each adjudicator gets one attempt. A file that is not valid JSON may be repaired syntactically without changing any decision; entries naming genes absent from the genome or reactions absent from the model are not applied and are reported. Decisions are committed and included in the inputs freeze before any outcome is downloaded.

On the development organisms the same procedure, run blind in the same way (packets from the U′ models, one fresh subagent each, 3 October 2026), produced 296 decisions for 296 candidates (18 R6, 62 R1, 2 R2+R1, 214 abstentions; all applicable). The procedure was not revised afterwards (section 8).

The adjudicator's own training may contain published phenotypes (notably for *M. tuberculosis* H37Rv, whose in vitro gene essentiality is widely published). It is told not to use them and to abstain if it notices it is relying on them, but this cannot be verified. Section 5.4 therefore includes a sensitivity analysis without MycoTube.

## 4. Scoring protocol (unchanged from development)

`gembench.protocols.carbon_fitness_generic` via `scripts/run_transfer_study.py run`, with carbon uptake −10 mmol/gDW/h, growth threshold 0.001 /h, fitness threshold −2 (inclusive), no rich-medium filter, medium completion by uptake-only transporters (not for pantothenate, folate or bicarbonate), the energy-generating-cycle gate, GLPK. Conditions are carbon-source experiments grouped by compound × base medium with fitness averaged over replicate experiments. A gene is predicted important when its deletion drops growth below the threshold in a condition where the wild-type model grows; it is observed important when its fitness is at most −2.

## 5. Outcomes and analysis

### 5.1 Per organism
For each organism and each comparison (U′ vs B0; M vs U′; and the secondary pairs), `scripts/transfer_compare.py` reports wild-type growth coverage of both arms (conditions where the model grows / mapped conditions), the number of paired gene–condition observations, the MCC of each arm on those observations, the **paired MCC difference** with a gene-bootstrap 95% interval (1,000 resamples, seed 0), and the predictions that changed with their agreement with the observed fitness. Conditions are those where both arms' wild types grow.

The **primary gene set is the union** of the two arms' scored genes: a gene that is absent from one arm's model is predicted to have no effect in that arm, which is what that model implies. The common-gene comparison (genes present in both models, as in `gembench.comparison`) is secondary. The union was made primary after the development comparison of M with U′ showed why it matters: in MR1 the blind R6 assignments added genes that the common-gene comparison cannot see (common genes +0.009, union −0.009). For arms that add no scored gene the two are identical.

### 5.2 Aggregate (primary)
`scripts/transfer_aggregate.py`: the unweighted mean of the per-organism paired differences over evaluable organisms, a percentile organism-level bootstrap 95% interval (10,000 resamples, seed 0), the numbers of organisms improved, worsened and unchanged (|difference| < 0.001) and the exact two-sided sign test over organisms that changed. Every organism counts once.

### 5.3 Pre-declared reading
For each primary comparison: **supported** if the mean difference is positive and the organism-bootstrap interval excludes zero; **not supported** if the mean is zero or negative; **inconclusive** otherwise. With six organisms the bootstrap interval is coarse and the sign test can only reach p < 0.05 when all six move in the same direction; both are reported, and the per-organism table is the main display. The result is reported whatever it is.

### 5.4 Secondary and sensitivity analyses
Attribution of U′ to its parts (B0→Uonly_rxn, B0→Uonly_atp, U→UN, U→UQ, UN→UNQ) and of M to its decision kinds (U′→M_R6, U′→M_R1R2); the common-gene version of every comparison; changes in wild-type coverage; the M vs U′ comparison without MycoTube; the decisions applied per organism and their changed predictions.

## 6. Failures and missing data
- An organism whose base model cannot be gap-filled, or whose model grows in no mapped condition in either arm of a comparison, is **not evaluable** for that comparison and is reported as such; it is not replaced and its difference is not imputed.
- Conditions on undefined media, with unmapped carbon sources or with exchanges absent from the model are excluded and counted.
- Solver failures are unknown results, never essentiality calls. An undefined MCC stays undefined.
- A software defect found after unblinding is fixed only in a new, separately reported run; the original result is kept and reported first.

## 7. Time order and custody
This evaluation is **self-custodied**: the same model family designs the method, prepares the inputs and runs the evaluation. Separation rests on recorded order: (1) this method freeze, committed and recorded in a Claude project document; (2) panel A inputs prepared (media, carbon sources, base models, packets, decisions) and frozen in an inputs freeze, also recorded externally; (3) only then are panel A's fitness tables downloaded, with their hashes and download time recorded; (4) all arms are run once. Panel B's fitness tables are not downloaded.

## 8. Development evidence (retrospective, four organisms)
From `results/transfer_v1/development/` (`paired.json`, `aggregate_union.json`, `aggregate_common.json`). These organisms shaped the rules, so the numbers show what the method does where it was built, not how well it transfers.

Wild-type growth coverage and union-gene MCC (MCC of B0 and U′ from their comparison, of M from the U′–M comparison; each on the conditions where both arms of that comparison grow):

| Organism | Conditions grown B0 → U′ → M (of mapped) | MCC B0 | MCC U′ | MCC M |
|---|---|---|---|---|
| MR1 | 8 → 9 → 9 (of 12) | 0.470 | 0.492 | 0.479 |
| Btheta | 14 → 14 → 14 (of 25) | 0.497 | 0.510 | 0.551 |
| Smeli | 21 → 22 → 22 (of 33) | 0.565 | 0.617 | 0.639 |
| Putida | 28 → 28 → 28 (of 43) | 0.496 | 0.553 | 0.590 |

Paired MCC differences (union genes; gene-bootstrap 95% interval):

| Comparison | MR1 | Btheta | Smeli | Putida | Mean (organism bootstrap) | Up / down / same |
|---|---|---|---|---|---|---|
| **B0 → U′** | +0.022 [0.000, 0.051] | +0.013 [−0.014, 0.050] | +0.051 [0.017, 0.096] | +0.057 [0.017, 0.104] | **+0.036** [0.017, 0.054] | 4 / 0 / 0 |
| **U′ → M** | −0.009 [−0.035, 0.018] | +0.056 [0.001, 0.111] | +0.028 [−0.001, 0.061] | +0.037 [0.003, 0.080] | **+0.028** [0.003, 0.049] | 3 / 1 / 0 |
| B0 → reaction patches only | +0.016 | −0.006 | +0.020 | +0.025 | +0.014 | 3 / 1 / 0 |
| B0 → ATP synthase rule only | 0.000 | 0.000 | +0.006 | 0.000 | +0.002 | 1 / 0 / 3 |
| U → UN (normalisation) | +0.006 | +0.019 | +0.010 | 0.000 | +0.009 | 3 / 0 / 1 |
| U → UQ (menaquinone rule) | 0.000 | 0.000 | +0.015 | +0.026 | +0.010 | 2 / 0 / 2 |
| U′ → M_R6 (assignments only) | −0.017 | +0.050 | +0.009 | +0.001 | +0.011 | 2 / 1 / 1 |
| U′ → M_R1R2 (removals and joins only) | +0.008 | +0.007 | +0.019 | +0.037 | +0.018 | 4 / 0 / 0 |

On common genes U′ → M is +0.009, +0.027, +0.019 and +0.037 (mean +0.023, 4 / 0 / 0); the difference from the union comes from genes that M adds to the models. In MR1 the copper ABC transporter assigned to the gene-less copper uptake reaction (SO0486–SO0488) makes three genes predicted essential in all nine conditions, and none has a fitness defect in any. The added genes elsewhere were mostly right: ketol-acid reductoisomerase (BT2074, 14 of 14 conditions), ornithine carbamoyltransferase (BT3717, 14 of 14) and L-aspartate oxidase (BT3184, 11 of 14) in Btheta, and anthranilate synthase (SMc02725, 22 of 22) in Smeli. The Btheta result is not fully blind: the procedure's worked example is the Btheta dihydroorotate dehydrogenase assignment (BT0892 + BT0891) from development cycle 6, which was made with the fitness data in view, and the Btheta adjudicator made that assignment. The examples name only Btheta genes, so they cannot leak answers for panel organisms.

Two observations contradict the development note quoted in the procedure ("removing apparently false alternatives … was usually wrong"): under the procedure's restrictions, removals and joins improved all four organisms, and assignments were the mixed part. The procedure was nevertheless **not revised**, so that the panel tests exactly the procedure that the development adjudicators used; the note stays in the procedure as written.

## 9. Known limits
- Six organisms make a small panel; a true effect of the size seen in development (~0.01–0.03 MCC) could be missed, and a chance pattern could look consistent.
- The panel comes from the same experimental platform (Fitness Browser) and the same draft collection as development; transfer to other assays or to current CarveMe output is not tested.
- Correction rules were developed on four organisms and are partly specific to the defects seen there.
- Self-custody (section 7) and the adjudicator's unknown training exposure (section 3.6).
