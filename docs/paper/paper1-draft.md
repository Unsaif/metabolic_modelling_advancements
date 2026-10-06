# Paper 1 draft: do corrections to draft metabolic models transfer?

> Repository copy of the Claude Docs draft (https://claude.ai/code/artifact/58c661c0-bff9-419c-a2e4-bbf1a6e4b3c3), exported 6 October 2026 after the hand-curated comparison was added. The doc is the working version; this copy records the text at export.

Oct 5, 2026 · @Tim Hulshof

## Note for Tim

This is a first full draft of the benchmark paper, written only from the locked transfer study and its held-back second round. Apart from the hand-curated comparison added later, no number in it is new: each one comes from a result file in the repository (branch `claude/opus-continuation`), listed under Data and code availability. Numbers derived for the paper (gene counts, t-intervals) are in results/transfer\_v1/paper1\_derived\_numbers.json.

The abstract was rewritten with Tim on 6 October, and "frozen" became "locked" throughout. The abstract's 94% is 367 of the 391 changed predictions ((273 + 94) / (274 + 117)), from counts in the derived-numbers file.

Later on 6 October we added a comparison with hand-curated models and coverage figures (descriptive, after the study). These numbers are new. An independent agent recomputed them from the raw files, and its corrections, mostly to wording, are applied (`results/transfer_v1/reference_models/independent_check/`). The exact times moved to Supplementary Table S1.

Decisions that are yours:

- [ ] Target journal. The structure is journal-neutral (Introduction, Results, Discussion, Methods); candidates include Molecular Systems Biology, PLOS Computational Biology and Genome Biology.
- [ ] Authors, and the statement of how AI was used (the curation itself and most of the analysis were done by Claude).
- [ ] Whether to add a robustness check with current CarveMe drafts of the same twelve organisms (robustness only: these organisms are now exposed).
- [ ] Whether to contact Bernstein et al. to reconcile their reference numbers before citing them (roadmap Q7).
- [ ] Licence check for redistributing Fitness Browser data, EMBL GEMs drafts and the CarveMe universe.
- [ ] Confirm that the curation packets' gene descriptions are the original annotations, not the Fitness Browser's fitness-based re-annotations (this bears on the blinding of the AI curator).
- [ ] For future studies, a public registration (for example on OSF) instead of a private project record.

## Title and abstract

Working titles:

1. Do corrections to draft metabolic models transfer to new organisms? A pre-specified test on twelve bacteria
2. Blind curation of automatically reconstructed metabolic models, tested on held-out bacteria against genome-wide mutant screens

**Abstract.** Most bacteria have only automatically built metabolic models, and corrections to these models are usually judged on the same data that inspired them. We asked whether corrections developed on four bacteria also improve the models of other bacteria.

As the yardstick we used published genome-wide mutant screens (RB-TnSeq), which show which genes a bacterium needs to grow on a given carbon source. We compared them with the genes that automatically built models (2017 CarveMe drafts from the EMBL GEMs collection) predict to be needed. On the four development bacteria we wrote five automatic correction rules that use only the genome annotation, and a procedure in which an AI model corrects gene–reaction links without seeing any screening data. We locked both methods, with a time-stamped record, before downloading the screens of twelve new bacteria: all remaining bacteria in the screening database that met criteria set in advance. We tested them in two rounds of six, holding the second round back until the first had been analysed.

Ten of the twelve could be tested. Before correction, the models agreed with the screens only moderately (Matthews correlation coefficient 0.33 to 0.52 in the eight bacteria whose experiments they could represent; 0 is chance and 1 is perfect). The automatic rules raised it by 0.037 on average (95% interval 0.025 to 0.048). The gain was the same in both rounds and in development, nine bacteria improved and none got worse, and 94% of the predictions the rules changed now agree with the screens. AI curation added 0.028 (eight improved, two worse; 0.029 on the held-back round alone), almost entirely by removing genes the drafts wrongly listed as backup enzymes. Four bacteria exposed limits of the locked pipeline. For two, the locked rules could not set up a model to test. Two others had experiments the models cannot represent: one grows by photosynthesis, and one "carbon-free" medium contains organic carbon. All four are reported.

In a comparison added after the study, we scored hand-curated models with the same protocol. None approached the agreement between replicate screens (MCC 0.85 to 0.94): the extensively curated *E. coli* model iML1515 scored 0.59. In *P. putida*, on the genes and conditions both models cover, the draft scored about the same as a hand-curated model before and after correction (0.53 and 0.57 against 0.56; differences within about ±0.1), and 60% of the corrected draft's wrong calls were also made by the curated model. The curated model grew in more conditions (36 of the 43 screened conditions the pipeline can map, against 28). Across the twelve organisms whose drafts represent the experiments, the drafts grew in 150 of 267 such conditions, and the corrections barely changed this.

Simple annotation-based corrections carry over to new bacteria. Their gain is small, on a scale where neither the drafts nor the hand-curated models we scored come close to the limit set by the screens' reproducibility. A second gap, larger for the drafts, is coverage: they cannot grow on many of the carbon sources the bacteria use.

## Introduction

Corrections to automatically built metabolic models are rarely tested on organisms they were not developed on. This study runs that test.

Genome-scale metabolic models link genes to reactions and predict which genes an organism needs to grow in a given environment. Automated pipelines such as CarveMe \[1\], gapseq \[2\] and ModelSEED \[3\] now build most available models. The EMBL GEMs collection alone holds 5,587 CarveMe drafts \[1\]. Drafts from different tools disagree with each other \[4,5\] and with measured phenotypes \[2\], and curation is what turns a draft into a useful model \[6,7\].

How do we know that a correction helps? Usually by scoring it against phenotype data, often the same data that suggested it. A correction chosen because it fixes a disagreement will look good on the data that revealed that disagreement, whether or not it generalizes. Structural quality scores do not settle the question either, because they measure a model's structure rather than its predictions \[8\].

Randomly barcoded transposon sequencing (RB-TnSeq) measures the fitness of mutants in nearly every gene, for dozens of bacteria across many carbon sources \[9,10\]. Bernstein et al. showed how to use such data to evaluate *Escherichia coli* models \[11\]. Because the data cover many organisms, they allow a stronger test: develop corrections on some organisms, then score them on organisms that played no part.

AI models are entering curation as well. For Human2 (Human-GEM 2.0), GPT-4 assessed all 26,246 gene–reaction pairs; the pairs it flagged were reviewed by hand, and the updated model was evaluated against CRISPR gene essentiality \[12\]. Machine-learning gap-filling has also been shown to improve phenotype predictions of draft models \[15\]. To our knowledge, no fully automated correction procedure has been locked before any outcome was accessed, tested on several held-out organisms with an isolated, paired effect estimate, and then repeated on a held-back set.

Here we developed five annotation-based correction rules and a written AI curation procedure on four bacteria, and locked both: every file they use was fingerprinted and recorded with a timestamp, so any later change would show. We then scored them once on six new bacteria whose fitness data were downloaded only after the lock, and repeated the test unchanged on six more that had been held back. We report what transferred, what did not, and where the locked pipeline itself limited the test, and we compare the drafts with a hand-curated model scored the same way.

## Results

### Study design: develop on four bacteria, test on twelve

We developed the corrections on four bacteria, locked them, and tested them on twelve new bacteria in two rounds (Figure 1).

**Panel.** Eligibility was declared before anything was downloaded:

- Fitness Browser bacteria not used before in this project;
- an exact-strain draft in the EMBL GEMs collection;
- at least eight carbon-source experiments;
- at most three base media.

Twelve organisms qualified: every remaining bacterium in the Fitness Browser that met these criteria. A seeded random split assigned six to panel A, tested first, and six to panel B, held back until panel A had been analysed. Until each panel's inputs were locked, only non-outcome data were downloaded: gene and experiment metadata, protein sequences and the draft models.

**Arms.** Three models were built per organism:

- **B0:** the draft with a minimal gap-fill for growth on one reference condition, chosen by a written rule;
- **U′:** B0 plus five automatic, annotation-based rules;
- **M:** U′ plus gene-rule decisions made blind by an AI curator.

**Primary metric.** The paired difference in MCC between two arms, on the union of their scored genes. Per organism it carries a gene-bootstrap 95% interval. Across organisms it is the mean of the per-organism differences, with an organism-bootstrap 95% interval. The pre-declared reading: **supported** if the mean is positive and its interval excludes zero, **not supported** if the mean is zero or negative, and **inconclusive** otherwise.

**Locking.** The method, and then each panel's inputs including every AI curation decision, were locked before that panel's fitness data were downloaded. Locking means that every file was fingerprinted (hashed), committed to git and recorded with an external timestamp, so any later change would show. The full time order and its independent check are in Methods.

![Figure 1. Design and order of events of the transfer study](fig1_timeline.png)

*Figure 1. Design and order of events of the transfer study, 3 October 2026. The panel was drawn from metadata only. The organism-specific corrections behind the development arms were made in September. Their annotation-triggered form, and the arms shown here, were written and run on 3 October, after the panel was drawn and before any panel outcome was downloaded. Times of the locks and downloads are in Supplementary Table S1.*

### The drafts agree only moderately with mutant fitness

Before any correction, the drafts predict gene importance with MCC 0.33–0.52 in eight of the ten evaluable new organisms (Table 1).

- **The outliers.** *R. palustris* (0.145) and *M. tuberculosis* (−0.016) are well below that range; their models misrepresent the experiments.
- **Gap-fill.** Most drafts needed a small gap-fill to grow on their reference condition: 0 to 7 reactions.
- **Coverage.** Most models grew in only part of the mapped conditions, in as few as 7 of 21 (*D. shibae*) or 3 of 7 (*R. palustris*). *M. tuberculosis* grew in all 19.
- **Corrections and coverage.** The corrections changed the number of conditions with growth by at most one, and only in three new organisms.

*Table 1. Organisms of the transfer study. Conditions mapped: compound × medium conditions whose carbon source maps to a BiGG metabolite; some lack an exchange reaction in the model and cannot support growth. Grows: conditions in which the model grows, for B0 and U′. MCC: union-gene MCC of B0 and U′ on the H1 comparison and of M on the H2 comparison, each on its own paired observations.*

| Phase | Organism | Gap-fill reactions | Conditions mapped | Grows B0 / U′ | Union genes | MCC B0 / U′ / M | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Development | *Bacteroides thetaiotaomicron* VPI-5482 | 5 | 25 | 14 / 14 | 507 | 0.497 / 0.510 / 0.551 | evaluable |
| Development | *Pseudomonas putida* KT2440 | 2 | 43 | 28 / 28 | 1,038 | 0.496 / 0.553 / 0.590 | evaluable |
| Development | *Shewanella oneidensis* MR-1 | 3 | 12 | 8 / 9 | 713 | 0.470 / 0.492 / 0.479 | evaluable |
| Development | *Sinorhizobium meliloti* 1021 | 1 | 33 | 21 / 22 | 979 | 0.565 / 0.617 / 0.639 | evaluable |
| Panel A | *Echinicola vietnamensis* DSM 17526 | 0 | 14 | 7 / 7 | 497 | 0.440 / 0.476 / 0.533 | evaluable |
| Panel A | *Dinoroseobacter shibae* DFL 12 | 7 | 21 | 7 / 7 | 584 | 0.422 / 0.472 / 0.530 | evaluable |
| Panel A | *Dyella japonica* UNC79MFTsu3.2 | 3 | 9 | 6 / 6 | 542 | 0.502 / 0.568 / 0.597 | evaluable |
| Panel A | *Mycobacterium tuberculosis* H37Rv | 0 | 19 | 19 / 19 | 463 | −0.016 / −0.016 / −0.018 | evaluable |
| Panel A | *Rhodopseudomonas palustris* CGA009 | 5 | 7 | 3 / 3 | 679 | 0.145 / 0.176 / 0.145 | evaluable |
| Panel A | *Dechlorosoma suillum* PS | — | — | — | — | — | not evaluable: gap-fill failed |
| Panel B | *Caulobacter crescentus* NA1000 | 4 | 13 | 5 / 5 | 547 | 0.431 / 0.448 / 0.475 | evaluable |
| Panel B | *Cupriavidus basilensis* 4G11 | 1 | 36 | 15 / 16 | 958 | 0.373 / 0.435 / 0.486 | evaluable |
| Panel B | *Marinobacter adhaerens* HP15 | 4 | 29 | 16 / 16 | 615 | 0.329 / 0.373 / 0.424 | evaluable |
| Panel B | *Shewanella loihica* PV-4 | 3 | 12 | 9 / 10 | 645 | 0.483 / 0.511 / 0.515 | evaluable |
| Panel B | *Shewanella amazonensis* SB2B | 1 | 20 | 14 / 15 | 648 | 0.518 / 0.549 / 0.560 | evaluable |
| Panel B | *Desulfovibrio vulgaris* Miyazaki F | — | — | — | — | — | not evaluable: no reference condition |

### A hand-curated *P. putida* model scores about the same on shared genes but grows in more conditions

This comparison was added after the study, to put the size of the gains in context. It is descriptive, not a held-out test. Both organisms involved were used in development, and iJN1463 itself was consulted then: it was the first curated model scored on this benchmark, and development audits of the menaquinone pathway drew on its biomass composition and pathway. One of U′'s reaction patches (carbamate kinase) also fixes an error that iJN1463 shares.

We scored two hand-curated models with the study's protocol, without gap-filling them: iML1515 for *E. coli* \[17\] and iJN1463 for *P. putida* \[18\]. Of the study organisms, only *P. putida* and *M. tuberculosis* have hand-curated models in BiGG \[13\], and the pipeline cannot represent the *M. tuberculosis* screens. Published curated models of *B. thetaiotaomicron*, *S. oneidensis* and *S. meliloti* exist outside BiGG and were not scored here.

- **No model comes close to the screens' own reproducibility.** Calls from replicate screens of the same condition agree with each other at MCC 0.90 to 0.94 (0.85 to 0.88 between replicates from different experiment sets). iML1515 grows in all 32 *E. coli* conditions but scores 0.59. The limit lies in the models and the binary growth call, not in noise in the screens.
- **On shared genes, the draft scores about the same as the curated model, before and after correction.** On the 764 genes (765 for M) and 28 conditions both cover:

| Model of *P. putida* KT2440 | MCC, shared genes and conditions | Conditions with growth, of 43 |
| --- | --- | --- |
| Draft (B0) | 0.53 | 28 |
| Draft + automatic rules (U′) | 0.57 | 28 |
| Draft + automatic rules + AI curation (M) | 0.61 | 28 |
| Hand-curated iJN1463 | 0.56 | 36 |

- **The differences are within noise.** The curated model minus the draft is +0.035 \[−0.063, 0.140\] before correction, −0.006 \[−0.108, 0.103\] after the automatic rules and −0.042 \[−0.139, 0.054\] after both steps (gene bootstrap). The errors overlap: 60% of the corrected draft's wrong calls on these genes and conditions (710 of 1,190) are also wrong in the curated model, about eight times what independent errors would give.
- **On the union of genes, the drafts score higher.** On the 383 genes only iJN1463 contains, 1,235 of its 1,287 "important" calls are not confirmed by the screens. Counting every gene of either model, as the study's primary analysis does, the curated model scores 0.45 against 0.49 to 0.58 for the draft arms. The shared-gene view isolates the genes both models represent; the union view favours the drafts.
- **The curated model grows in more conditions.** Of the 57 carbon-source conditions in which *P. putida* was screened, 43 have a carbon source the pipeline can map. iJN1463 grows in 36 of these and the draft in 28, before and after correction. In 10 of the draft's 15 gaps and all 7 of the curated model's, the model has no exchange reaction for the carbon source.

**Coverage across the study.** In the twelve organisms whose drafts represent the experiments, the drafts grow in 33% to 75% of the screened conditions the pipeline can map (median 60%; 150 of 267 in total, or 47% of all 321 screened carbon-source conditions), and the corrections raise this only to 155 of 267. The gene-level scores in this paper cover only the conditions in which a model grows (Table 1).

### The automatic rules transfer, and the transfer replicates

The five automatic rules improved agreement on new organisms by the same amount as in development: +0.036 in development, +0.037 on panel A and +0.036 on panel B (Figure 2, left).

| Set | Mean MCC difference, U′ − B0 \[organism bootstrap 95%\] | Up / down / unchanged | Sign test p | Reading |
| --- | --- | --- | --- | --- |
| Development, 4 organisms (retrospective) | +0.036 \[0.017, 0.054\] | 4 / 0 / 0 | — | — |
| Panel A, 5 evaluable | +0.037 \[0.017, 0.054\] | 4 / 0 / 1 | 0.125 | supported |
| Panel B, 5 evaluable | +0.036 \[0.024, 0.051\] | 5 / 0 / 0 | 0.063 | supported: replicated |
| Panels A and B, 10 evaluable | +0.037 \[0.025, 0.048\] | 9 / 0 / 1 | 0.004 | supported |

- **No organism got worse.** The one unchanged organism is *M. tuberculosis*, whose model carries no signal in any arm (Results, last section).
- **Same on common genes.** Restricted to genes present in both models, the pooled estimate is identical (+0.037 \[0.025, 0.048\]).
- **Same with t-intervals.** With five organisms per panel the percentile bootstrap is coarse. t-based 95% intervals give the same readings: panel A \[0.007, 0.067\], panel B \[0.015, 0.057\], pooled \[0.022, 0.051\].
- **What changed.** Across the ten new organisms the gain rests on 53 genes whose predictions changed (2 to 10 per organism). Of 274 gene × condition calls that switched from important to unimportant, 273 agree with the data; of 117 that switched the other way, 94 do.

**Which rules carried it (Figure 3).**

- **Menaquinone rule.** This rule removes menaquinone from the biomass when the annotation shows fewer than two of the eight steps of its pathway. It contributed the most on both new panels: +0.020 each (panel A \[0.003, 0.041\]; panel B \[0.003, 0.040\]), against +0.010 in development. It fired in six of the ten new organisms and never lowered MCC. All 200 calls it changed went from important to unimportant, and 199 agree with the data. They involve 24 genes, mostly not menaquinone enzymes: thioesterases, β-oxidation and 2-methylcitrate-cycle enzymes and 4-hydroxybenzoyl-CoA reductase subunits, which the drafts had drawn into a route to menaquinone.
- **Universal reaction patches.** +0.005 on panel A and +0.012 on panel B, against +0.014 in development.
- **Smaller contributions.** Rule normalisation added +0.004 and +0.003, and the ATP synthase rule +0.003 and +0.001.

**Stable total, shifting parts.** The total gain was the same in all three phases, but its make-up shifted: the reaction patches and rule normalisation contributed less on the new panels than in development, and the menaquinone rule more. The gain is consistent but modest. After correction, MCC stays between 0.37 and 0.57 in the eight organisms with a working model.

![Figure 2. Per-organism paired MCC differences](fig2_paired_differences.png)

*Figure 2. Paired MCC differences per organism, union of genes. Left: automatic rules (U′ vs B0, H1). Right: blind AI curation (M vs U′, H2); the two panels have different x-axis ranges. Circles: organisms, with 95% gene-bootstrap intervals. Diamonds: means over organisms, with 95% organism-bootstrap intervals. Grey: development organisms (retrospective). Blue: panel A. Red: panel B. Black: pooled new organisms. D. suillum and D. vulgaris Miyazaki F were not evaluable under the locked rules.*

![Figure 3. Components of the corrections](fig3_attribution.png)

*Figure 3. Contribution of each component, by phase. Each row compares an arm with the arm it extends, on the union of genes. Dots: mean over organisms. Bars: 95% organism-bootstrap interval. Ticks: individual organisms.*

### Blind AI curation helps, through one kind of decision

Blind AI curation of gene rules added +0.028 MCC pooled over the ten new organisms (Figure 2, right). Almost all of the gain came from removing gene alternatives that are not catalysts.

**How the curator worked.** It never saw fitness data. For each organism, a fresh AI subagent received a packet built from the U′ model and gene annotations only. The packet listed the gene-less reactions and the alternative-gene (OR) rules of reactions essential in at least one condition. Every gene came with its Fitness Browser description, RefSeq definition and genomic neighbourhood. The subagent was the same Claude model as the study session and made one attempt, under a fixed written procedure. It was instructed to open only its packet and not to use the web; its tools were not technically restricted. Four outcomes were allowed:

- **removal** (R1): remove an alternative that is not a catalyst;
- **join** (R2): join subunits into a complex;
- **assignment** (R6): assign genes to a gene-less reaction;
- abstain.

| Set | Candidates | R1 | R2 | R6 | Abstained |
| --- | --- | --- | --- | --- | --- |
| Development (4 organisms) | 296 | 62 | 2 (with R1) | 18 | 214 |
| Panel A (5 organisms) | 265 | 76 | 2 | 25 | 162 |
| Panel B (5 organisms) | 297 | 74 | 0 | 12 | 211 |

**Results by panel.**

- **Panel A: inconclusive.** +0.027 \[−0.007, 0.060\], three organisms up and two down. The two losses were *M. tuberculosis* (−0.002) and *R. palustris* (−0.031), the two organisms whose models misrepresent the experiments. Without *M. tuberculosis* (a pre-declared sensitivity analysis, because its phenotypes are widely published), the estimate was +0.035 \[−0.008, 0.066\].
- **Panel B: supported.** +0.029 \[0.011, 0.047\], five of five up. This is the clean replication: the pooled analysis below was declared after the panel A results were known.
- **Pooled: supported.** +0.028 \[0.009, 0.046\], eight up and two down (sign test p = 0.11). On common genes: +0.026 \[0.013, 0.040\].
- **t-intervals.** Panel A \[−0.026, 0.081\], panel B \[0.001, 0.057\], pooled \[0.005, 0.051\]: the same readings, with panel B and the pooled estimate close to zero at their lower end.

**Which decisions helped.**

- **Removals** carried the gain: removals and joins together gave pooled +0.026 \[0.014, 0.038\], eight up, one down, one unchanged. On the new panels the curator joined subunits only twice, both in *M. tuberculosis*, so the effect of joins is untested.
- **Assignments** gave nothing on balance: pooled +0.003 \[−0.008, 0.013\], three up, three down, four unchanged.

**Precision of the new calls.** Curation only adds "important gene" calls: it never turned an important call into an unimportant one. So the question is how many of the new calls the fitness data confirm (Figure 4). The counts are gene × condition calls, so one gene can count many times: the 328 calls from removals and joins involve 29 genes.

- **R1/R2.** On the new panels, these decisions created 328 new calls (gene × condition) and the data confirm 174. Of the 154 misses, 95 are in *M. tuberculosis*, whose model has no signal. Without it, 174 of 233 are confirmed (75%).
- **R6.** Assignments created 85 new calls, of which 34 are confirmed; on panel B only 3 of 32.
- **For comparison.** Between 0.6% and 16% of gene × condition cells are important in the data. U′'s own important calls are 44–80% precise in the organisms whose models carry signal.

**Ion transporters.** Two thirds of the unconfirmed assignment calls (34 of 51) involve ion transporters: a zinc ABC transporter in *R. palustris*, chloride channels in *C. crescentus* and *S. loihica*, and a manganese transporter in *S. loihica*. These genes became important in every condition with growth, and none is. The other 17 involve dihydropyrimidinase, the glycine-cleavage T and H proteins, dihydroorotate dehydrogenase and a lactate permease. A secondary hypothesis added before panel B (H3), after a look at three organisms already scored, dropped assignments to inorganic-ion transport reactions. It pointed the right way (+0.003 \[0.000, 0.008\]) but was inconclusive, because the filter changed only two organisms.

![Figure 4. Precision of curation decisions](fig4_curation_precision.png)

*Figure 4. New "important gene" calls made by blind curation, by decision kind and organism. Filled: confirmed by the fitness data (fitness ≤ −2). Hatched: not confirmed. Counts are gene × condition cells on the union of genes.*

### The locked pipeline set its own limits

Applied as locked, the pipeline could not set up a model for two organisms and misrepresented the experiments of two more. None was replaced or dropped after the outcomes were seen.

| Organism | What happened | Rule responsible |
| --- | --- | --- |
| *D. suillum* PS (panel A) | All six gap-fill solutions found within the six-round limit created an energy-generating cycle, so no base model could be built. Its fitness table was never downloaded. | gap-fill with an energy-cycle gate |
| *D. vulgaris* Miyazaki F (panel B) | The carbon-source experiments of its main medium name no compound, so no reference condition could be chosen. Its fitness table was never downloaded. | reference-condition rule |
| *M. tuberculosis* H37Rv (panel A) | The "no carbon" Sauton's medium still contains asparagine, citrate (as ferric ammonium citrate) and 0.2% ethanol. As defined organics these are unlimited, so the model grows at implausible rates (5–9 per hour) in every condition and its predictions carry no signal (MCC about 0 in every arm). | media rule |
| *R. palustris* CGA009 (panel A) | The experiments are phototrophic and anaerobic, but the CarveMe universe has no photosystem. The gap-fill made the model respire the trace nitrate of the medium's cobalt nitrate, so it grows in 3 of 7 conditions with the wrong energy metabolism (MCC 0.14–0.18). | media rule; no light input |

Each problem suggests a fix:

- limit counter-ions that come only from trace-element salts;
- flag organic carbon in "no carbon" media;
- handle phototrophs;
- add more gap-fill cut rounds;
- use a reference rule that skips media without named carbon sources.

Testing such fixes needs new organisms, because every eligible organism has now been used.

## Discussion

Corrections derived from four bacteria improved predictions for nine of ten held-out bacteria and left the tenth unchanged, by the same amount in development and in both rounds (+0.036, +0.037 and +0.036). This holds within one assay platform and one 2017 draft collection.

**How large is the gain?** In absolute terms it is small: 0.037 on a scale where the drafts start at 0.33 to 0.52. For comparison, no model we scored comes close to the agreement between replicate screens (0.85 to 0.94). The extensively curated *E. coli* model scores 0.59. In *P. putida*, a hand-curated model scores about the same as the draft on shared genes, and 60% of the corrected draft's wrong calls are also the curated model's. Much of the remaining disagreement is therefore not specific to automatic drafts, and since replicate screens agree closely, it lies in what the models represent rather than in noise in the data. Against that background the gain is modest but real: no new organism got worse, the total was the same in all three phases, and 94% of the changed predictions moved toward the data. A set of rules tuned to the quirks of its development organisms would shrink on new ones; this set did not, although its parts shifted between phases (Figure 3).

**A second gap is coverage.** The drafts grow in only about half to two thirds of the screened conditions the pipeline can map, and the corrections barely change that; the hand-curated *P. putida* model grows in more of them. Filling these gaps would use condition-level growth, as the base gap-fill already does for one condition, but no gene-level outcome data. In *P. putida*, 10 of the draft's 15 gaps lack an exchange reaction for the carbon source altogether, and 14 of the 57 screened carbon sources cannot be mapped at all. A locked test of such a fix on new organisms is the natural next study.

**What transferred is unglamorous.** On the new organisms the largest single contribution came from a biomass rule: dropping menaquinone where the annotation shows fewer than two of the eight steps of its pathway. The calls it corrected were mostly not about menaquinone enzymes. They concerned enzymes that the drafts had drawn into a route to menaquinone to satisfy the biomass demand. A universal biomass that demands what an organism cannot make distorts which genes look essential. The rule's annotation trigger has not yet been checked against measured quinone types, which differ systematically between bacterial groups \[16\].

**AI curation works within limits.** Without seeing outcomes, an AI curator improved predictions on eight of ten new organisms; the clean replication on panel B gave +0.029. The gain came from one kind of decision: removing gene alternatives that are not catalysts. Assigning genes to gene-less reactions did not help on balance, and two thirds of its wrong calls involved ion transporters. That suggests a practical policy: let AI curators remove alternatives, and hold their gene assignments, especially to ion transporters, to a stricter standard. Joins of subunits were too rare here to judge.

**Pre-specification is cheap.** It took hours, not months, to do four things: lock a method, record the lock, download outcomes only afterwards and have an independent agent recompute every number. The same design suits other claims that a model change helps. Here it has three weaknesses:

- **Self-custody.** The same model family designed, prepared and scored the test, and checked it.
- **A private record.** The external time record is a private project document. A public registration would be stronger.
- **One manifest changed.** The replication plan's manifest no longer verifies at the current head, because one input file gained panel B entries before the panel B inputs were locked.

An external custodian of the outcome data, or a community challenge on newly measured organisms, would remove the first weakness.

**Limits.**

- **One assay platform and one draft collection.** RB-TnSeq fitness and the 2017 EMBL GEMs drafts; current CarveMe output and other reconstruction tools are untested.
- **Fitness is not growth.** Pooled competitive fitness is not single-mutant growth. Only genes with fitness values are scored, which excludes genes essential in the library's growth medium, so the benchmark measures carbon-source-specific importance. The threshold of −2 was used without a t-score filter.
- **The metric.** MCC on binary calls was pre-specified. Bernstein et al. recommend the area under the precision–recall curve \[11\], which we have not analysed on matched observations.
- **Small panels.** Ten evaluable organisms, two of them *Shewanella* species related to the development organism *S. oneidensis*.
- **Development exposure.** The curated *E. coli* model helped flag one reaction patch, so *E. coli* counts as a fifth exposed organism.
- **Blinding of the curator.** It rests on instructions. The packets' gene descriptions came from the Fitness Browser gene table. Whether any of them reflect the Browser's fitness-based re-annotations remains to be checked.
- **Unknown training exposure.** The AI curator's training may include published phenotypes; the pre-declared sensitivity analysis without *M. tuberculosis* did not change the reading.
- **Locked pipeline choices.** These excluded or misrepresented four of twelve organisms.
- **One curated comparison.** Only *P. putida* had a hand-curated model in BiGG whose screens the pipeline can represent; published curated models of three other study organisms were not scored. *P. putida* is a development organism, iJN1463 was consulted during development, and the comparison was made after the study.

**Next tests.**

- Growth on every screened carbon source that maps to BiGG: fill the coverage gaps from condition-level growth alone, lock the fix, and test gene predictions on new organisms.
- Hand-curated models of more study organisms (*B. thetaiotaomicron*, *S. oneidensis*, *S. meliloti*), scored the same way.
- Current CarveMe drafts of the same organisms (a robustness check, since these organisms are now exposed).
- A random-removal control for the curator's removals.
- Other phenotype types, such as growth profiles.
- Organisms whose fitness data are published after a method is locked.

## Methods

**Data.**

- **Fitness.** Fitness Browser RB-TnSeq gene fitness and experiment metadata \[9,10\]. Carbon-source experiments with no second condition other than DMSO were grouped by compound and base medium, and fitness was averaged over replicate experiments.
- **Models.** EMBL GEMs drafts \[1\] pinned at collection commit 260d0f1, and the CarveMe bacterial universe.
- **Media recipes.** FEBA medium recipes (bitbucket commit 803273865ea9).

**Panel selection.** Eligibility criteria were declared before any candidate's data were downloaded:

- in the Fitness Browser on the selection date;
- not among the nine organisms with archived downloads in the project;
- a bacterium;
- an exact-strain EMBL GEMs model, after normalising punctuation and known taxonomic renamings;
- at least eight carbon-source experiments;
- at most three base media.

Eligible organisms were ordered by identifier and shuffled with a fixed seed (20261003): the first half became panel A and the rest panel B. No organism was replaced after assignment.

**Gene mapping.** Model genes were matched to Fitness Browser loci by protein sequence:

- exact identity;
- or containment covering at least 90% of the longer sequence;
- or equal length with at least 97% identical positions.

Ambiguous matches were left unmapped. Between 91.4% and 100% of model proteins mapped (panel A at least 95.6%).

**Media and carbon sources.** Each base medium was mapped from its FEBA recipe to BiGG metabolites \[13\]:

- inorganic salts to their ions;
- vitamins, cofactors, nucleobases and reductants to trace uptake;
- other defined organics unlimited;
- buffers and chelators dropped;
- complex ingredients made a medium undefined, and its experiments were excluded.

Water, protons, CO2 and trace metals (capped at 0.1 mmol/gDW/h) were always available. O2 was available when the experiments were aerobic. Carbon sources were mapped to the BiGG metabolite of exactly the named compound. Applied to the development media, the rule reproduced the hand-curated media component by component, except for two trace limits (cysteine in MOPS rich medium, methionine in Varel–Bryant medium).

**Base model (B0).** Each draft received a minimal gap-fill from the CarveMe universe (minimum growth 0.05 per hour), with a gate that rejects solutions creating energy-generating cycles. The gap-fill targeted one reference condition chosen by a written rule: in the base medium with most carbon-source experiments, the first tested substrate in a fixed order (D-glucose, D-fructose, glycerol, L-lactate, pyruvate, succinate, …). The rule reproduced the hand-chosen references of all four development organisms.

**Automatic rules (U′).** Applied in order:

1. three universal reaction patches: carbamate kinase in its catabolic direction, subunit conjunctions for carbamoyl-phosphate synthase, and an irreversible MECDPDH4E;
2. normalisation of conjunction rules;
3. diffusion and exchange of glycolaldehyde;
4. an F0F1 ATP synthase when the draft has none and the annotation names at least six of its eight subunit types;
5. removal of menaquinone from the biomass when the annotation shows fewer than two steps of its pathway.

The rules generalise organism-specific corrections made in September. One reaction patch was found partly with the curated *E. coli* model iML1515. Their annotation-triggered form was written on 3 October, after the panel had been drawn from metadata and before any panel outcome was downloaded. U′ was chosen after all development arms had been run, as the arm with the highest development MCC.

**Blind adjudication (M).** For each organism, a packet was built from the U′ model and metadata only. It listed the gene-less reactions and the alternative-gene rules of reactions essential in at least one mapped condition where U′ grows. Each came with its equation, annotations and current rule. Every gene came with its description from the Fitness Browser gene table, its RefSeq definition and its genomic neighbourhood. One fresh subagent per organism wrote its decisions in one attempt. It was instructed to read only its packet and the written procedure and not to use the web. The procedure was not revised after development. All decisions were locked with the inputs before any outcome was downloaded.

**Scoring.**

- **Simulation.** Carbon uptake was set to −10 mmol/gDW/h and growth counted above 0.001 per hour. Medium components that the model could not take up received an uptake-only transporter (medium completion), and the energy-generating-cycle gate applied. The solver was GLPK through COBRApy 0.32.1 (Python 3.11).
- **Calls.** A gene was predicted important when its deletion stops growth in a condition where the wild-type model grows. It was observed important when its fitness is −2 or lower, with no t-score filter. Only genes with fitness values are scored.
- **Metric.** MCC \[14\] on all gene × condition pairs where both arms grow. The primary gene set was the union of both arms' genes: a gene absent from a model counts as predicted unimportant there.
- **Uncertainty and reading.** Per-organism intervals: 1,000 gene-bootstrap resamples. Means: 10,000 organism-bootstrap resamples, with t-based intervals as a sensitivity analysis. Organisms changing by less than 0.001 counted as unchanged. Sign tests were exact and two-sided. The pre-declared reading is given in Results.

**Locking, time order and verification.** Each lock is a manifest of file hashes with a content fingerprint, committed to git and recorded in a private project document whose server timestamp follows the commit by less than half a minute. Outcome files were downloaded only after the corresponding lock, and their hashes and download times were recorded. All locks, downloads and runs took place on 3 October 2026. The method was locked first, and each panel's inputs were locked before that panel's fitness tables were first downloaded; each arm was then run once. Supplementary Table S1 lists the times.

All four manifests verify at their own commits. At the current head the replication plan's manifest differs in one file, `filter_report.json`, which gained panel B entries before the panel B inputs were locked; its earlier entries are unchanged.

On 3 October a separate agent, using its own code, recomputed the primary and secondary results from the run matrices, with exact agreement. It checked the fitness values of every arm against the downloaded tables (376,620 arm × gene × condition cells on panel A and 881,911 on panel B), and confirmed the time order from git history, the manifests, the download records and the external timestamps. A second agent checked this paper's derived tables, figures and counts against the result files on 5 October.

**Hand-curated reference models (added after the study).** iML1515 \[17\], from the validation repository of Bernstein et al. \[11\], and iJN1463 from BiGG (BiGG's version of the iJN1462 model of Nogales et al. \[18\]) were scored with the study's media and carbon-source tables, parameters, medium completion and solver, without gap-filling. Their gene identifiers are locus tags that match the Fitness Browser's. iML1515 was scored on the *E. coli* BW25113 screens; the strain's deletions (araBAD, rhaBAD, lacZ) were not applied, because the protocol has no strain adjustment, and the three of these genes that are scored show no fitness defect. The runs used Python 3.13 with the study's versions of COBRApy and optlang (the study ran under Python 3.11). The *P. putida* models were compared on the genes both cover and the conditions in which both grow, as in the study's common-gene analysis, and on the union of genes. Coverage is the number of mapped conditions in which a model grows. Replicate agreement compares the calls (fitness ≤ −2) of pairs of replicate experiments of the same compound and medium, on the same genes and conditions; single replicates give a lower bound on the reliability of the averaged calls the benchmark scores.

## Data and code availability

All code, inputs, lock manifests, run matrices and curation decisions are in the project repository (github.com/Unsaif/metabolic\_modelling\_advancements, branch `claude/opus-continuation`). It should be archived with a DOI before submission.

| What | Where in the repository |
| --- | --- |
| Method, replication plan and adjudication procedure (with the verbatim prompt) | `docs/studies/transfer-method-v1.md`, `transfer-v1-replication-plan.md`, `transfer-adjudication-procedure-v1.md`, `transfer-adjudication-prompt-v1.txt` |
| Lock manifests | `results/study_freezes/transfer_v1_method.json`, `transfer_v1_inputs_panel_A.json`, `transfer_v1_replication_plan.json`, `transfer_v1_inputs_panel_B.json` |
| Every curation decision | `data/studies/transfer_v1/decisions/<organism>.json` |
| Run matrices, paired comparisons and aggregates | `results/transfer_v1/development/`, `evaluation_panel_A/`, `evaluation_panel_B/`, `pooled_panels_A_B_*.json` |
| Fitness download records (hashes and times) | `data/fitness_browser_panel/` |
| Hand-curated reference models and comparison | `scripts/score_reference_model.py`, `compare_reference_models.py`; `results/transfer_v1/reference_models/` (iML1515 itself from github.com/dbernste/E\_coli\_GEM\_validation, not redistributed) |
| Figures and Table 1 | `scripts/plot_transfer_timeline.py`, `plot_transfer_paired.py`, `plot_transfer_attribution.py`, `plot_curation_precision.py`, `make_transfer_table1.py` |

Fitness data come from the Fitness Browser (fit.genomics.lbl.gov) and draft models from the EMBL GEMs collection.

## References

A reviewing agent checked entries 1–14 against publisher or repository records on 5 October 2026 (Chicco and Bernstein from search records only). Entries 15 and 16 were added on its advice, and 17 and 18 with the hand-curated comparison.

1. Machado D, Andrejev S, Tramontano M, Patil KR. Fast automated reconstruction of genome-scale metabolic models for microbial species and communities. *Nucleic Acids Res.* 2018;46(15):7542–7553. doi:10.1093/nar/gky537
2. Zimmermann J, Kaleta C, Waschina S. gapseq: informed prediction of bacterial metabolic pathways and reconstruction of accurate metabolic models. *Genome Biol.* 2021;22:81. doi:10.1186/s13059-021-02295-1
3. Henry CS, DeJongh M, Best AA, Frybarger PM, Linsay B, Stevens RL. High-throughput generation, optimization and analysis of genome-scale metabolic models. *Nat Biotechnol.* 2010;28:977–982. doi:10.1038/nbt.1672
4. Mendoza SN, Olivier BG, Molenaar D, Teusink B. A systematic assessment of current genome-scale metabolic reconstruction tools. *Genome Biol.* 2019;20:158. doi:10.1186/s13059-019-1769-1
5. Hsieh YE, Tandon K, Verbruggen H, Nikoloski Z. Comparative analysis of metabolic models of microbial communities reconstructed from automated tools and consensus approaches. *npj Syst Biol Appl.* 2024;10:54. doi:10.1038/s41540-024-00384-y
6. Thiele I, Palsson BØ. A protocol for generating a high-quality genome-scale metabolic reconstruction. *Nat Protoc.* 2010;5:93–121. doi:10.1038/nprot.2009.203
7. Heinken A, et al. Genome-scale metabolic reconstruction of 7,302 human microorganisms for personalized medicine. *Nat Biotechnol.* 2023;41:1320–1331. doi:10.1038/s41587-022-01628-0
8. Lieven C, et al. MEMOTE for standardized genome-scale metabolic model testing. *Nat Biotechnol.* 2020;38:272–276. doi:10.1038/s41587-020-0446-y
9. Wetmore KM, et al. Rapid quantification of mutant fitness in diverse bacteria by sequencing randomly bar-coded transposons. *mBio.* 2015;6:e00306-15. doi:10.1128/mBio.00306-15
10. Price MN, et al. Mutant phenotypes for thousands of bacterial genes of unknown function. *Nature.* 2018;557:503–509. doi:10.1038/s41586-018-0124-0
11. Bernstein DB, Akkas B, Price MN, Arkin AP. Evaluating *E. coli* genome-scale metabolic model accuracy with high-throughput mutant fitness data. *Mol Syst Biol.* 2023;19:e11566. doi:10.15252/msb.202311566
12. Luo J, Wang H, Moyer D, Guo Z, Robinson JL, Gustafsson J, Anton M, Chen Y, Kerkhoven EJ, Nielsen J, Li F. Reconstruction of human metabolic models with large language models. *Proc Natl Acad Sci USA.* 2026;123(15):e2516511123. doi:10.1073/pnas.2516511123
13. King ZA, et al. BiGG Models: a platform for integrating, standardizing and sharing genome-scale models. *Nucleic Acids Res.* 2016;44:D515–D522. doi:10.1093/nar/gkv1049
14. Chicco D, Jurman G. The advantages of the Matthews correlation coefficient (MCC) over F1 score and accuracy in binary classification evaluation. *BMC Genomics.* 2020;21:6. doi:10.1186/s12864-019-6413-7
15. Chen C, Liao C, Liu YY. Teasing out missing reactions in genome-scale metabolic networks through hypergraph learning. *Nat Commun.* 2023;14:2375. doi:10.1038/s41467-023-38110-7
16. Collins MD, Jones D. Distribution of isoprenoid quinone structural types in bacteria and their taxonomic implication. *Microbiol Rev.* 1981;45(2):316–354. doi:10.1128/mr.45.2.316-354.1981
17. Monk JM, Lloyd CJ, Brunk E, et al. iML1515, a knowledgebase that computes *Escherichia coli* traits. *Nat Biotechnol.* 2017;35:904–908. doi:10.1038/nbt.3956
18. Nogales J, Mueller J, Gudmundsson S, et al. High-quality genome-scale metabolic modelling of *Pseudomonas putida* highlights its broad metabolic capabilities. *Environ Microbiol.* 2020;22(1):255–269. doi:10.1111/1462-2920.14843

## Supplementary information

**Supplementary Table S1. Time order of the transfer study, 3 October 2026 (UTC).** The external record is the server timestamp of the private project document made for each lock. Times are truncated to the second.

| Step | Event | Time | External record |
| --- | --- | --- | --- |
| 1 | Method locked | 09:40:45 | 09:41:01 |
| 2 | Panel A inputs locked (media, carbon sources, base models, curation decisions) | 10:11:43 | 10:11:50 |
| 3 | Panel A fitness tables first downloaded; every arm then run once | 10:12:18 to 10:12:40 | — |
| 4 | H3 arm explored on three organisms already scored (one development, two panel A) | 12:28:54 to 12:29:25 | — |
| 5 | Replication plan for panel B locked | 12:31:09 | 12:31:17 |
| 6 | Panel B inputs locked | 13:06:21 | 13:06:45 |
| 7 | Panel B fitness tables first downloaded; every arm then run once | 13:07:00 to 13:07:19 | — |
