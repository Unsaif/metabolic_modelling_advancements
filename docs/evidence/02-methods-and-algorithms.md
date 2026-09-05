# Evidence brief 02 — Methods and algorithms (state as of September 2026)

Compiled 4 Sept 2026 by Claude research agents for the "Metabolic modelling improvements" project. Every claim carries a source URL; items marked *unverified* or "details unverified" could not be confirmed in-session (several publisher hosts blocked full-text fetches) and should be checked before being cited. This brief feeds the Landscape & Roadmap document (00-landscape-and-roadmap.md).

## (A) Key developments 2023–2026

**Sampling / degeneracy**
- **dingo** (2024): Multiphase Monte Carlo + billiard-walk-with-rounding; samples Recon3D "in less than a day" where PolyRound+hopsy CDHR "did not converge" after 10 days. https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae037/7633919
- **LooplessFluxSampler** (2024): ADSB sampler for the *loopless* space; ~4× faster than ll-ACHRB on small models, ~1000× on iMM904 (MATLAB/COBRA). https://link.springer.com/article/10.1186/s12859-023-05616-2
- Constrained **Riemannian HMC** used for community flux sampling; sampling vs. FBA shifted predicted antagonistic interactions 74%→61% and cooperative 30%→44% (anaerobic). https://link.springer.com/article/10.1186/s12859-024-05655-3
- COBRApy 0.32 added a **CHRR** sampler (uniformity guarantee); 0.29 added a HiGHS+OSQP hybrid LP/QP solver. https://opencobra.github.io/cobrapy/releases/
- Review of sampling + context-specific models (Trends Biotechnol. 2025; details unverified). https://pubmed.ncbi.nlm.nih.gov/40781001/

**Uncertainty / Bayesian / ensembles**
- **BayFlux** (2023): MCMC over genome-scale flux posteriors given 13C data; matches 13CFLUX2 point estimates on *E. coli*; genome-scale models gave narrower posteriors than core models. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011111
- **Trans-dimensional Bayesian model averaging for 13C-MFA** (arXiv May 2026): RJ-MCMC + diffusive nested sampling over "billions of model variants"; synthetic validation only. https://arxiv.org/abs/2605.25079
- **GEMsembler** (bioRxiv 2025): ensembles across CarveMe/gapseq/etc.; *E. coli* essentiality AUCPR 0.556→0.711 (iML1515 gold standard 0.767). https://www.biorxiv.org/content/10.1101/2025.04.09.647174v1.full

**Omics integration**
- **Localgini thresholding** (2025): benchmarked with FASTCORE/MBA/iMAT/INIT/GIMME/mCADRE on NCI-60 + HPA over Recon2.2/Recon3D/Human1; thresholding choice reduced inter-method variance more than the extraction method itself. https://www.nature.com/articles/s41540-025-00617-8
- **RNA-seq normalisation benchmark** (2024): between-sample normalisation (RLE/TMM/GeTMM) gave ~0.80 accuracy (AD) vs. within-sample (TPM/FPKM), at the cost of missed true positives. https://www.nature.com/articles/s41540-024-00448-z
- Expression-weighted FBA in *Arabidopsis* cut error vs. 13C-MFA from 169–180% to 10–13%. https://academic.oup.com/bioinformatics/article/39/5/btad186/7114031
- 2026 review (Mol. Omics) reorganises methods as switches/valves/allocation/thermodynamic anchors and states head-to-head comparative studies "remain sparse". https://academic.oup.com/molecular-omics/article/22/2/aaiag005/8469221
- Older benchmarks still cited: Opdam et al. 2017 https://www.cell.com/fulltext/S2405-4712(17)30010-8; "Guidelines for extracting… context-specific models" (Metab. Eng.; details unverified) https://www.sciencedirect.com/science/article/abs/pii/S1096717622001525

**Enzyme/proteome constraints**
- **GECKO 3.0** (Nat. Protoc. 2024): DLKcat integration, "light" ecModels, proteomics constraints; MATLAB+RAVEN+Gurobi; latest release v3.2.5 (Mar 2026). https://www.nature.com/articles/s41596-023-00931-7 ; https://github.com/SysBioChalmers/GECKO
- **ECMpy 2.0** (Python; DLKcat/AutoPACMEN/BRENDA kcats), release Dec 2023, 19 stars. https://github.com/tibbdc/ECMpy
- **kcat predictors**: TurNuP (2023) test R²=0.44, <40% identity R²=0.33 vs. DLKcat ≈0; ~18% better proteome-allocation predictions than DLKcat in yeast ecModels. https://www.nature.com/articles/s41467-023-39840-4 — CatPred (2025) kcat R²=0.61 held-out / 0.39 OOD, with ensemble uncertainty. https://www.nature.com/articles/s41467-025-57215-9 — **WILDkCAT** (2026) pipeline: experimental kcats cover only 29.4% of iML1515 enzyme–reaction pairs; gaps filled with CataPro. https://academic.oup.com/bioinformatics/article/42/8/btag510/8732673
- **ME-models**: coralME built 495 gut ME-models (Cell Systems, Nov 2025) https://today.ucsd.edu/story/uc-san-diego-researchers-develop-new-tool-to-predict-how-bacteria-influence-health ; *P. putida* ME-model: pathway-level r=0.71/0.76 vs. RNA-seq/Ribo-seq (M-model 0.50/0.49), solved with Quad MINOS at 1e-16 tolerance, ~5 min vs. 72 ms. https://www.nature.com/articles/s41540-025-00521-1

**Kinetic / hybrid**
- **RENAISSANCE** (Nat. Catal. 2024): NES-trained generators of kinetic parameter sets (502 params, 123 reactions) — ~92% valid models; steady-state data only. https://www.nature.com/articles/s41929-024-01220-6
- Review of genome-scale kinetic modelling (ACS Synth. Biol. 2025; details unverified). https://pubs.acs.org/doi/10.1021/acssynbio.4c00868

**Communities / dFBA**
- **PhyloCOBRA** (2025): merges phylogenetically close GEMs before MICOM/OptCom; 186 metagenomes; ~50% runtime cut; better growth-vs-replication-rate correlation at order level. https://academic.oup.com/bioinformatics/article/41/7/btaf328/8211154
- Multi-objective host–microbiota coupling (iScience 2024): epithelial cell + 5-organism ecosystem. https://www.cell.com/iscience/fulltext/S2589-0042(24)01317-8
- **PyCoMo** (2024) https://academic.oup.com/bioinformatics/article/40/4/btae153/7635576 ; **BN-BacArena** (2024; 19 species, 44 foods train/11 test) https://academic.oup.com/bioinformatics/article/40/5/btae266/7660539 ; **MetaBiome** agent-based+GEM (2025) https://journals.asm.org/doi/full/10.1128/msystems.01652-24
- Reconstruction-tool choice dominates predicted exchanges: "exchanged metabolites was more influenced by the reconstruction approach rather than the… community". https://www.nature.com/articles/s41540-024-00384-y
- Scale resource: 247,092 human-microbe GEMs (Cell Systems 2025). https://www.cell.com/cell-systems/fulltext/S2405-4712(25)00029-8

**Whole-body / multiscale**
- Infant WBMs (Cell Metab. 2024): 360 organ-resolved, sex-specific models (days 0–180), growth matched WHO curves, 10,000 personalised newborns, IEM biomarker trajectories. https://research.universityofgalway.ie/en/publications/personalized-metabolic-whole-body-models-for-newborns-and-infants/
- **Human2** (PNAS 2026): LLM + GitHub-Actions-assisted curation; demographic-specific organ models; "enzyme-constrained dynamic" whole-body simulation of feeding/fasting. https://research.chalmers.se/en/publication/551634
- Harvey/Harvetta in use (81,094/83,521 reactions; CPLEX; deterministic pFBA, no experimental validation). https://www.frontiersin.org/journals/physiology/articles/10.3389/fphys.2025.1586750/full
- COBRA Toolbox v3.7 "Persephone Pipeline Publication Release" (quadratic-LP interface renewal, WBM support). https://github.com/opencobra/cobratoolbox/releases
- Numerics baseline: ME-models span ~6 orders of magnitude; DQQ (double→quad MINOS) reaches 1e-20 tolerances in 4.5 h on 70k-variable *E. coli* ME. https://www.nature.com/articles/srep40863

**Solvers**
- Machado 2024 benchmark: Gurobi fastest, CPLEX second, HiGHS/SCIP competitive on LP; CPLEX ~100× slowdown for >4-member communities; GLPK/COIN failed Recon3D MILP within a week. https://pmc.ncbi.nlm.nih.gov/articles/PMC10878033/
- No peer-reviewed application of PDLP/cuPDLP GPU LP to GEMs was found (unverified absence).

## (B) Used in practice vs. published-and-abandoned

- **Actively maintained**: COBRApy (0.32.1; HiGHS default-capable) https://github.com/opencobra/cobrapy ; MICOM 0.39.1 released 11 Aug 2026 https://pypi.org/project/micom/ ; GECKO v3.2.5 (Mar 2026); COBRA Toolbox v3.7; COMETS (2025 citation, 28 stars) https://github.com/segrelab/comets .
- **Dormant/legacy**: MASSpy last release 0.1.7, Dec 2022, Python ≤3.9, 27 stars https://pypi.org/project/masspy/ ; SKiMpy pins `cobra<=0.24.0` (current 0.32) https://github.com/EPFL-LCSB/skimpy ; AutoPACMEN last PyPI 0.6.1 https://pypi.org/project/autopacmen-Paulocracy/0.6.1/ ; ECMpy 2.0 no release since Dec 2023, 19 stars; COBRAme superseded by coralME.
- **Omics methods**: the 2026 review calls GIMME/iMAT "canonical" and GECKO "the most widely adopted" enzyme-constraint implementation; PROM/TRFBA/MADE appear only in reviews, not in 2023–26 benchmarks found. https://academic.oup.com/molecular-omics/article/22/2/aaiag005/8469221
- **Samplers**: ACHR/OptGP (COBRApy defaults) and CHRR (COBRA Toolbox, now COBRApy) dominate; dingo/ADSB are 2024 entrants with limited uptake evidence (unverified).
- **Loop removal**: CycleFreeFlux (2015) and ll-FBA remain the cited tools; no 2023–26 successor with adoption found. https://academic.oup.com/bioinformatics/article/31/13/2159/195895

## (C) Limitations, critiques, failure modes

- Extraction-method outputs depend more on thresholding/normalisation than on the algorithm (Localgini 2025; normalisation benchmark 2024); between-sample normalisation trades recall for precision.
- Reconstruction-tool bias propagates into community predictions (exchanged metabolites tool-dependent, https://www.nature.com/articles/s41540-024-00384-y).
- FBA optimisation vs. sampling yields qualitatively different interaction calls (74%→61% antagonism), i.e., community conclusions are objective-function artefacts.
- kcat ML: DLKcat R²≈0 for <40% identity enzymes (TurNuP paper); CatPred lists database noise, substrate-mapping errors, and that 3D features add little; only ~29% of iML1515 pairs have experimental kcats (WILDkCAT). Thus ecModels rest mostly on predicted or fitted parameters.
- ME/WBM numerics: quad precision needed (1e-16 tolerance; 5 min per solve for one ME-model; coralME's 495 models imply heavy compute); WBM studies are deterministic single-subject runs with "no experimental validation" (Alessi 2025).
- Validation statistics: GECKO 3 protocol concedes GEMs "are unable to correctly predict" many phenotypes; FlowGAT beats FBA on essential genes but underperforms on non-essential due to class imbalance — typical essentiality benchmarks are imbalanced and organism/condition-specific. Kaste & Shachar-Hill 2024 argue for formal model-validation/selection statistics in FBA/MFA (details unverified). https://aiche.onlinelibrary.wiley.com/doi/full/10.1002/btpr.3413
- Solver sensitivity: CPLEX 100× slowdown for community models; open-source MILP failures on Recon3D.

## (D) ML/LLM/AI work and evaluation

- **AMN** (2023): FBA embedded in trainable network; Q²=0.77–0.81 on *E. coli* media/knockouts. https://www.nature.com/articles/s41467-023-40380-0 → **dAMN** (2026): dynamic hybrid, median R²≈0.98 forecast / 0.97 novel media (*E. coli*, 280 media), 0.90–0.94 *P. putida*; code https://github.com/brsynth/dAMN-main-release
- **FlowGAT** (2024): GAT on mass-flow graphs; PRAUC ~0.85 (glucose), 255 labelled nodes, 11 carbon sources. https://www.nature.com/articles/s41540-024-00348-2
- **SIMBA-GNN** (2025): 2,850 GEM-simulated pairwise interactions as edge features; Spearman 0.85 vs. <0.6 baselines on 186 subjects. https://www.nature.com/articles/s41540-025-00631-w
- **FBA-informed surrogates for cybergenetic control** (2024): *in silico* only. https://arxiv.org/abs/2401.00670
- **RENAISSANCE**: generative kinetic parameters (above). Related 2026 work: EINN (title-only, unverified) https://www.sciencedirect.com/science/article/pii/S2405805X26001845
- **LLMs**: ChatGEM (PNNL, 2026): multi-agent RAG + COBRApy; RAG raised task score 2.63→4.20; one experimentally validated *P. putida* succinate prediction; backbones/code not in abstract. https://www.biorxiv.org/content/10.64898/2026.07.20.739662v1 — "Comprehensive evaluation of LLM capabilities for… GEMs in metabolic engineering" (2026; content unverified) https://www.biorxiv.org/content/10.64898/2026.06.03.730004v1 — Human2 used LLM-assisted curation (validation metrics not in abstract).

## (E) Opportunities for frontier-model intelligence

1. **A living, executable benchmark for omics-integration methods** (all extraction methods × thresholds × normalisations on Human1/Recon3D with 13C-MFA/CRISPR ground truth). Rationale: 2026 review says comparisons are sparse; 2025 work shows preprocessing dominates. *Success*: public leaderboard reproducing Localgini/Opdam results and ranking ≥10 methods with CIs.
2. **Numerically robust WBM/ME solving in Python**: implement DQQ-style double→quad refinement, scaling diagnostics and solution certification on HiGHS/Gurobi. *Success*: Harvey/Harvetta + coralME models solve to ≤1e-9 residuals without CPLEX, with runtime ≤2× current.
3. **Certified sampling at genome scale**: port ADSB/billiard-walk into COBRApy with convergence diagnostics (PSRF, ESS) and loopless option. *Success*: Recon3D loopless samples in <1 day, diagnostics passing, replacing ACHR/OptGP defaults.
4. **Uncertainty-aware ecModels**: propagate CatPred/TurNuP predictive variances (not point kcats) into ecModels via sampling/robust LP. *Success*: proteome-allocation prediction with calibrated intervals covering ≥90% of measured proteomics in yeast/E. coli.
5. **Auto-derivation of consistent community formulations**: formal proof/verification of MICOM/SteadyCom/cFBA objective equivalences and degeneracy, plus PhyloCOBRA-style reduction. *Success*: growth-vs-replication-rate correlation ≥ MICOM on the 186-sample set at ≤50% cost.
6. **Validation-statistics standard**: derive MCC/PRAUC with class-imbalance corrections, cross-condition holdouts and 13C-MFA agreement tests; ship as a MEMOTE-like check. *Success*: adopted in ≥3 tool papers within a year.
7. **LLM agent for model debugging** (unbounded loops, blocked reactions, infeasibility diagnosis) with objective test suites, extending ChatGEM's approach with grounded, reproducible evaluations. *Success*: ≥90% correct diagnoses on a curated infeasibility/loop corpus.

## (F) Most important sources

1. dingo (2024) https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae037/7633919
2. LooplessFluxSampler (2024) https://link.springer.com/article/10.1186/s12859-023-05616-2
3. BayFlux (2023) https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011111
4. Localgini benchmark (2025) https://www.nature.com/articles/s41540-025-00617-8
5. RNA-seq normalisation benchmark (2024) https://www.nature.com/articles/s41540-024-00448-z
6. Multi-omics integration review (2026) https://academic.oup.com/molecular-omics/article/22/2/aaiag005/8469221
7. GECKO 3.0 (2024) https://www.nature.com/articles/s41596-023-00931-7
8. TurNuP (2023) https://www.nature.com/articles/s41467-023-39840-4
9. CatPred (2025) https://www.nature.com/articles/s41467-025-57215-9
10. WILDkCAT (2026) https://academic.oup.com/bioinformatics/article/42/8/btag510/8732673
11. RENAISSANCE (2024) https://www.nature.com/articles/s41929-024-01220-6
12. Solver benchmark (2024) https://pmc.ncbi.nlm.nih.gov/articles/PMC10878033/
13. Quad-precision ME solving (2017) https://www.nature.com/articles/srep40863
14. Infant WBMs (2024) https://www.cell.com/cell-metabolism/fulltext/S1550-4131(24)00182-7
15. dAMN (2026) https://academic.oup.com/bioinformatics/article/42/5/btag230/8674704
