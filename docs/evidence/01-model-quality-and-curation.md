# Evidence brief 01 — Model quality, reconstruction and curation (state as of September 2026)

Compiled 4 Sept 2026 by Claude research agents for the "Metabolic modelling improvements" project. Every claim carries a source URL; items marked *unverified* could not be confirmed in-session and should be checked before being cited in any publication. This brief feeds the Landscape & Roadmap document (00-landscape-and-roadmap.md).

## (A) Key developments 2023–2026

- **AGORA2** (Jan 2023, Nat Biotech): 7,302 strain reconstructions via KBase drafts + DEMETER refinement; validated vs 3 datasets (metabolite uptake/secretion accuracy 0.81–0.84, enzyme activity 0.72, 253 drug–microbe pairs 0.81); avg MEMOTE-style score 73%; literature review of 732 papers. https://www.nature.com/articles/s41587-022-01628-0 . Post-publication patch: v2.01 replaced futile-cycle-enabling reactions in 151 reconstructions; an `sbml_files_fixed/` directory was added July 2024. https://www.vmh.life/files/reconstructions/AGORA2/version2.01/README.txt
- **APOLLO** (Cell Systems 2025; preprint Oct 2023): 247,092 MAG-based reconstructions (KBase + DEMETER). Quality at scale: mass/charge balance, stoichiometric & flux consistency (63.2% ± 6.2% flux-consistent reactions after refinement vs 60.9% for drafts), ATP/biomass sanity on aerobic/anaerobic media; agreement with experimental data 81%/73% (refined) vs 11%/10% (drafts). Limitations: smaller models than AGORA2, ~2,000 drug-metabolism reactions absent, MAG incompleteness. https://www.biorxiv.org/content/10.1101/2023.10.02.560573v1.full ; https://www.cell.com/cell-systems/fulltext/S2405-4712(25)00029-8
- **Bactabolize** (eLife, Oct 2023): reference-pangenome approach; *K. pneumoniae* KPPR1: 85.8% substrate-growth accuracy (190 substrates), 80.6% gene-essentiality accuracy (1,220 genes); ~98 s vs 5.46 h for gapseq; ModelSEED/KBase captured fewer genes. https://elifesciences.org/articles/87406
- **ModelSEED v2** (preprint Oct 2023): core-metabolism-first reconstruction with ATP-yield safety threshold, ML template assignment, 57% less gap-filling than v1; cyanobacteria and endosymbiont biomass remain weak. https://www.biorxiv.org/content/10.1101/2023.10.04.556561v1.full
- **Hsieh et al. 2024** (npj Syst Biol Appl): CarveMe vs gapseq vs KBase on identical genomes → reaction Jaccard only 0.23–0.24, metabolite 0.37; predicted exchanged metabolites depend more on the tool than on the community. https://www.nature.com/articles/s41540-024-00384-y
- **pan-Draft** (Genome Biol, Oct 2024): species-representative gapseq models from multiple MAGs; outperforms single-MAG models at 50–90% completeness; needs ≥15 genomes/species. https://link.springer.com/article/10.1186/s13059-024-03425-1
- **MACAW** (Genome Biol, Mar 2025): dead-end/dilution/duplicate/loop tests flag 40% of Human-GEM 1.15, 48% of yeast-GEM 9.0, 52% of iML1515 reactions; authors stress flags ≠ error rates. https://link.springer.com/article/10.1186/s13059-025-03533-6
- **Human2 / Human-GEM 2.0** (PNAS, Apr 2026): GPT-4-assisted GPR audit plus GitHub Action CI (see D). https://www.pnas.org/doi/10.1073/pnas.2516511123
- **gapseq 2.0** (Jan 2026) re-implemented alignment (one large MSA; BLAST/DIAMOND/MMseqs2); v2.1.0 (May 2026) updated to UniProt 2026_01. https://github.com/jotech/gapseq/releases
- **CarveMe-GutMicrobes** (bioRxiv, Jun 2026): gut-specific curated universal DB with taxonomic restriction; validated on gene essentiality and metabolite production (numbers not retrievable). https://www.biorxiv.org/content/10.64898/2026.06.26.734454v1.full
- **MetaNetX 4.5** (NAR, Jan 2026): 1.50 M metabolites, 83,796 reactions; added VMH namespace; dual-level protonation normalisation. https://academic.oup.com/nar/article/54/D1/D617/8343507
- **EHMN 2026** (Metabolites, Mar 2026): 22,642-reaction human network, 37 energy-generating cycles removed; 61% of species lack formulas so full ΔG parameterisation was impossible. https://www.mdpi.com/2218-1989/16/4/236
- **Thermodynamics tooling**: Thermo-Flux (bioRxiv Nov 2025) converted 87/107 BiGG models to thermodynamic-stoichiometric form, failing on 20 (multi-compartment transport, undefined structures) https://www.biorxiv.org/content/10.1101/2025.11.20.689566v1 ; ThermOptCOBRA (iScience 2025, MATLAB) https://github.com/NiravBhattLab/ThermOptCOBRA
- **BioModTool** (Bioinf Adv, Feb 2025): structured biomass objective functions from composition data (cannot infer composition). https://academic.oup.com/bioinformaticsadvances/article/5/1/vbaf036/8029660
- **DNNGIOR** (bioRxiv 2023): neural-network gap-filling; F1 = 0.85 for reactions in >30% of training genomes; 14× better than unweighted gap-filling on drafts. https://api.biorxiv.org/details/biorxiv/10.1101/2023.07.10.548314

## (B) Tools / resources and maintenance status

| Tool | Status (verified) |
|---|---|
| CarveMe | v1.6.6, 12 Sep 2025 (SCIP solver since 1.6.0) https://pypi.org/project/carveme/ |
| gapseq | v2.1.0, 30 May 2026; 216 stars https://github.com/jotech/gapseq |
| ModelSEED/KBase | MS2 is KBase's recommended path; v2 paper still preprint (journal status unverified) https://docs.kbase.us/workflows/metabolic-models ; DB: 33,978 cpds/36,645 rxns https://github.com/ModelSEED/ModelSEEDDatabase |
| Pathway Tools/MetaFlux | v29.0 (Apr 2025), 29.5 (Dec 2025), 30.0 listed; no MetaFlux changes since v23.5 (2019) https://bioinformatics.ai.sri.com/ptools/release-notes.html |
| Bactabolize | v1.0.5, 23 May 2025 https://github.com/kelwyres/Bactabolize |
| Reconstructor (KEGG-based) | v1.2.0, Oct 2025; moved to csbl/reconstructor https://github.com/emmamglass/reconstructor |
| merlin | v4 (NAR 2022); no 2023–26 release found (unverified) https://academic.oup.com/nar/article/50/11/6052/6606174 |
| AuReMe | no 2023–26 update found (unverified) |
| MEMOTE | last release 0.17.0, 8 Jan 2024; 134 open issues; recent commit activity unverified https://github.com/opencobra/memote/blob/develop/HISTORY.rst |
| MACAW | v1.0.0, Nov 2024 https://github.com/Devlin-Moyer/macaw |
| Human-GEM | 12,931 rxns / 8,461 mets / 2,848 genes; v2.0.0 tag (year rendered inconsistently across fetches — verify) https://github.com/SysBioChalmers/Human-GEM/releases/tag/v2.0.0 |
| Recon3D | 2018; no update found; superseded in practice by Human-GEM/VMH https://www.nature.com/articles/nbt.4072 |
| AGORA2 | v2.01 (Jan 2023) + SBML fixes Jul 2024; GitHub repo has 3 commits https://github.com/VirtualMetabolicHuman/AGORA2 |
| BiGG | legacy site frozen since Oct 2019; new bigg.bio lists 4,881 models (contents unverified) http://bigg.ucsd.edu/ ; https://bigg.bio/ |
| Rhea | release 142 (2 Sep 2026), 18,611 reactions; UniProt's annotation standard https://www.rhea-db.org/ |
| KEGG | FTP/bulk requires paid subscription; commercial licence required https://www.kegg.jp/kegg/legal.html |
| eQuilibrator 3.0 | ~1 M compounds via MetaNetX; group-contribution limits https://academic.oup.com/nar/article/50/D1/D603/6432058 |
| FROG (BioModels) | reproducibility standard (FVA, reaction/gene deletion, objective); 4 implementations https://www.ebi.ac.uk/biomodels/curation/fbc |

## (C) Limitations, critiques, failure modes

- **Tool disagreement dominates biology**: same genome → 23–24% reaction overlap between CarveMe/gapseq/KBase; community-exchange predictions are tool artefacts. https://www.nature.com/articles/s41540-024-00384-y
- **MEMOTE score ≠ predictive quality**: in *L. kefiri*, CarveMe had the highest MEMOTE score but a wrong fermentation profile; Bactabolize/PanGEM models did not grow on defined medium. https://www.mdpi.com/2218-1989/15/12/767 . MACAW authors warn flag counts overstate errors; resolution stays manual. https://link.springer.com/article/10.1186/s13059-025-03533-6
- **Draft quality is poor without curation**: APOLLO drafts matched experiments only 10–11% vs 73–81% after DEMETER. https://www.biorxiv.org/content/10.1101/2023.10.02.560573v1.full
- **Energy-generating cycles persist**: AGORA2 needed a post-release patch for 151 models; EHMN removed 37 EGCs; original EGC study https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005494
- **Annotation/GPR uncertainty**: 62% (16,277/26,246) of Human-GEM gene–reaction pairs were "inconclusive" even with GPT-4 + UniProt/HPA; 203 genes removed as non-metabolic. https://www.pnas.org/doi/10.1073/pnas.2516511123
- **Biomass & environment uncertainty**: swapping biomass between species changes up to 30% of predicted essential reactions; only 5/21 plant GEMs had measured biomass. https://link.springer.com/article/10.1186/s13059-021-02289-z ; sensitivity study https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/bit.70181 (content not fetched)
- **Namespace fragmentation**: MetaNetX must reconcile ≥13 databases with protonation/tautomer conflicts; KEGG is licence-restricted; legacy BiGG frozen; mergem needed for cross-namespace merging https://academic.oup.com/nargab/article/6/1/lqae010/7597502
- **Thermodynamic coverage gaps**: 20/107 BiGG models fail automatic conversion; 39–61% of human metabolites lack formulas.
- **Curation bottleneck**: "manual curation … remains labor-intensive" (Zuniga et al. 2023) https://www.preprints.org/manuscript/202311.0461 ; AGORA2 required 732 papers read manually.
- **No modern multi-tool benchmark**: the only broad systematic assessment remains Mendoza et al. 2019; no 2023–26 successor covering Biolog + essentiality across many organisms/tools was found. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1769-1

## (D) LLM / agentic work

- **Human2 (PNAS 2026)** — GPT-4 API with prompts tuned on 3,008 manually reviewed entries audited all 26,246 gene–reaction pairs (UniProt/HPA descriptions vs reaction role): 7,774 consistent, 2,195 inconsistent, 16,277 inconclusive. Manual review: 1,985 true / 210 false flags → 1,135 GPRs updated, 203 genes removed ("10% GPR consistency gain"). CI: testYamlConversion, 57 metabolic tasks, sanityCheck, memoteTest, macawTests block merges. Human2 MEMOTE 81%; IEM simulation accuracy 65.4% (ec version); improved essentiality MCC. Only a GPT-4 model was used; no precision/recall on the inconclusive majority. https://www.pnas.org/doi/10.1073/pnas.2516511123
- **Tewari et al., eLife (Jul 2026)** — frontier LLMs asked to generate networks de novo: core *E. coli* reactions recall 63.6–90.9%, but most generated networks gave infeasible flux/no growth (missing PGM, LDH); signalling-network perturbation accuracy 6–33% vs 78–95% for curated. https://elifesciences.org/articles/109709
- **Yeoh et al., bioRxiv (Jun 2026)** — GPT-4, Gemini, Claude, DeepSeek-R1 on GEM interpretation/coding: conversational prompting failed to find injected stoichiometric errors; tool-assisted (COBRApy) procedures with domain constraints succeeded; errors from identifier assumptions and context limits. https://api.biorxiv.org/details/biorxiv/10.64898/2026.06.03.730004
- **ChatGEM (PNNL, bioRxiv Jul 2026)** — multi-agent (ADEPT) + RAG over COBRApy; RAG raised code-generation score 2.63→4.20; experimentally validated succinate strains in *P. putida*. Simulation-focused, not curation. https://www.biorxiv.org/content/10.64898/2026.07.20.739662v1
- **MechAInistic (bioRxiv May–Jul 2026)** — Reviewer/Architect agents with executable COBRApy verification; drug-target case studies; "auditable chain" to literature. https://www.biorxiv.org/content/10.64898/2026.05.11.723319v2
- **Li et al. (bioRxiv Sep 2024)** — LLM mining of >29,000 metabolic-engineering entries (1,210 products, 751 organisms). https://api.biorxiv.org/details/biorxiv/10.1101/2024.09.09.612023
- **Gaps**: no published LLM gap-filling method, no LLM-driven de novo microbial reconstruction validated on Biolog/Tn-seq, no LLM evaluation of reaction mass/charge balancing at scale (searches returned none; DNNGIOR is the only ML gap-filler found).

## (E) Opportunities for frontier-model intelligence

1. **Evidence-graded GPR auditing** — Human2 left 62% of pairs inconclusive with GPT-4; retrieval over UniProt/Rhea/literature plus reasoning could resolve most. Evidence: precision/recall vs held-out manual decisions, inconclusive fraction <20%.
2. **Cross-namespace reaction reconciliation** — 23% tool overlap is largely identifier/representation mismatch; chemistry-aware matching (InChIKey, protonation, Rhea master reactions). Evidence: Jaccard on same-genome models rises and merged models pass MetaNetX/Rhea consistency checks.
3. **Flag-to-fix agents for MACAW/EGC/duplicate outputs** — propose specific directionality/transport/GPR fixes with citations. Evidence: curator acceptance rate, zero EGCs post-fix, essentiality MCC not degraded.
4. **Literature-to-evidence extraction with provenance** — replicate DEMETER's 732-paper review automatically for thousands of species. Evidence: recall of AGORA2/APOLLO phenotype tables from raw papers; new species reaching ≥0.8 experimental agreement without manual refinement.
5. **Closed-loop reconstruction benchmark** — agent reconstructs, tests against Biolog/Tn-seq, iterates; Yeoh & Tewari show tool-grounded loops work where prompting fails. Evidence: a public Mendoza-2019-style benchmark with held-out organisms.
6. **Structure/formula/charge assignment for orphan metabolites** — unblocks thermodynamics (Thermo-Flux failures, EHMN 61% formula gap). Evidence: formula coverage >95%, ΔG coverage, fewer thermodynamically infeasible loops.
7. **Organism-specific biomass/media from literature** — feed BioModTool-style templates. Evidence: predicted yields/essentiality shift within measured ranges.
8. **Agent-maintained CI for community GEMs** — generalise Human-GEM's Action stack (memote+macaw+tasks) to AGORA/APOLLO-scale collections. Evidence: adoption count, reduced post-release patches like AGORA2 v2.01.

## (F) Most important sources

1. Human2/LLM reconstruction, PNAS 2026 — https://www.pnas.org/doi/10.1073/pnas.2516511123
2. Benchmarking LLM-generated biochemical networks, eLife 2026 — https://elifesciences.org/articles/109709
3. LLM capabilities for GEM analysis, bioRxiv 2026 — https://www.biorxiv.org/content/10.64898/2026.06.03.730004v1
4. APOLLO, Cell Systems 2025 — https://www.cell.com/cell-systems/fulltext/S2405-4712(25)00029-8
5. AGORA2, Nat Biotech 2023 — https://www.nature.com/articles/s41587-022-01628-0
6. MACAW, Genome Biol 2025 — https://link.springer.com/article/10.1186/s13059-025-03533-6
7. Hsieh et al. tool comparison, npj SBA 2024 — https://www.nature.com/articles/s41540-024-00384-y
8. Bactabolize, eLife 2023 — https://elifesciences.org/articles/87406
9. ModelSEED v2, bioRxiv 2023 — https://www.biorxiv.org/content/10.1101/2023.10.04.556561v1.full
10. gapseq releases (v2.x, 2026) — https://github.com/jotech/gapseq/releases
11. MetaNetX 2026, NAR — https://academic.oup.com/nar/article/54/D1/D617/8343507
12. Bernstein et al. uncertainty review, Genome Biol 2021 — https://link.springer.com/article/10.1186/s13059-021-02289-z
13. Thermo-Flux, bioRxiv 2025 — https://www.biorxiv.org/content/10.1101/2025.11.20.689566v1
14. Sazonov et al. tool comparison, Metabolites 2025 — https://www.mdpi.com/2218-1989/15/12/767
15. Mendoza et al. systematic assessment, Genome Biol 2019 — https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1769-1

**Caveats**: Recon3D-specific 2025–26 updates, AuReMe/merlin status, and identifiers.org sustainability were not checked (marked unverified above); the CarveMe-GutMicrobes and Nègre et al. papers were only accessible as abstracts.
