# Evidence brief 06 — LLMs, agents and foundation models × constraint-based modelling (2023–2026)

Compiled 4 Sept 2026 by Claude research agents for the "Metabolic modelling improvements" project. Legend: **[PR]** peer-reviewed · **[PP]** preprint · **[GH]** GitHub/blog only · "unverified" = could not confirm (PNAS/PMC/bioRxiv full texts were blocked; abstracts and APIs were used). This brief feeds the Landscape & Roadmap document (00-landscape-and-roadmap.md).

## (A) Catalogue of LLM/agent work on GEMs

| Name | What it does | Evaluation | Status / URL |
|---|---|---|---|
| **Human2** (Luo, Wang, Moyer … Kerkhoven, Nielsen, F. Li; PNAS 123(15), Apr 2026) | Consensus human GEM that "leverages LLMs and GitHub Action checks to streamline automated, efficient, and collaborative curation"; demographic tissue models; enzyme-constrained dynamic whole-body model | See brief 01 for figures obtained from a separate fetch (62% of GPR pairs inconclusive after GPT-4 audit) | [PR] https://www.pnas.org/doi/10.1073/pnas.2516511123 ; Chalmers record https://research.chalmers.se/en/publication/551634 |
| Human-GEM release history | v1.14.0 (Feb 2021): "remove 149 non-metabolic genes, that are discovered using the GPT-3 language model"; v2.0.0: automated MACAW error summary on every PR | Release notes only | [GH] https://github.com/SysBioChalmers/Human-GEM/releases |
| **MechAInistic** (Loecker … Helikar, U. Nebraska; bioRxiv v1 May 2026, v4 Jul 2026; arXiv 2607.18249) | Architect–Reviewer multi-agent system turning NL questions into executable COBRApy workflows (pathway comparison, perturbation, drug targets) on paired cell states | Two case studies (RA naive B cells; MS Th17) proposing devimistat/CPI-613 and ivosidenib; claims plain LLMs give "plausible-sounding but unverifiable narratives"; no quantitative benchmark in abstract | [PP] https://www.biorxiv.org/content/10.64898/2026.05.11.723319v2 ; https://arxiv.org/abs/2607.18249 ; workflows https://github.com/drahmedabdeenhamed/mechanistic-workflows |
| **LLM capability evaluation for GEMs** (Yeoh, Patro, Wong, Poh; NUS; bioRxiv Jun 2026) | Systematic test of GPT-4, Gemini, Claude, DeepSeek-R1 on GEM interpretation + coding with a multi-LLM scoring rubric | DeepSeek-R1 best on theory, Gemini on coding; failures: context-window limits, "incorrect identifier assumptions", organism-specific reasoning errors; deliberately injected stoichiometric errors were poorly detected conversationally, better with domain constraints + tools | [PP] https://www.biorxiv.org/content/10.64898/2026.06.03.730004v1 |
| **D2Cell** (F. Li lab, Tsinghua; bioRxiv Sep 2024 → Trends Biotechnol 2026) | LLM mining of metabolic-engineering literature + hybrid DL-GEM target predictor + RAG chatbot | Fine-tuned Qwen1.5-14B: 85% precision / 86% F1; 74–79% accuracy on 100 test papers; 29,006 entries, 751 organisms; target prediction acc: E. coli 0.914 (0.804 OOD), yeast 0.842 (0.770), C. glutamicum 0.736 (0.671) | [PR] https://www.cell.com/trends/biotechnology/abstract/S0167-7799(26)00135-6 ; [PP] https://www.biorxiv.org/content/10.1101/2024.09.09.612023v1.full.pdf ; code https://github.com/LiLabTsinghua/D2Cell |
| **BioPromptX** (F. Li lab) | LLM extraction of gene function and enzyme kinetic parameters (Llama/GPT/Qwen/DeepSeek), RL-optimized prompts, validated vs UniProt/BRENDA | Paper unverified | [GH] https://github.com/LiLabTsinghua/BioPromptX |
| **MetaKnogic-Alpha** (OHSU; bioRxiv Feb 2026) | Hyper-relational KB from >100k articles; reasoning agent; every finding verified against a curated reaction network | 0.98 mechanistic accuracy when evidence exists | [PP] https://www.biorxiv.org/content/10.64898/2026.02.05.704050v1.full |
| **ModelSEEDagent** (J. Faria, Argonne) | NL agent: 24 tools (12 COBRApy, 3 ModelSEED build/gapfill, media tools, ID resolution) | Self-reported "100% success across 4 model types" (e_coli_core, iML1515…) | [GH] https://github.com/jplfaria/ModelSEEDagent |
| **gem-flux-mcp** (Faria) | GEM MCP server: build_media/build_model/gapfill_model/run_fba + ModelSEED DB lookups; 785 tests | No paper | [GH] https://github.com/jplfaria/gem-flux-mcp |
| **Biomni** (Stanford; bioRxiv May 2025) | General bio-agent; `perform_flux_balance_analysis` (COBRApy, SBML/JSON, constraints, top-10 fluxes) and mass-action `simulate_metabolic_network_perturbation` | Biomni-Eval1: 433 instances/10 tasks, none GEM-specific; **BiomniBench** (May 2026): 100 tasks, agents weak at "method selection, biological interpretation, and scientific reasoning" | [PP] https://www.biorxiv.org/content/10.1101/2025.05.30.656746v1 ; tools https://raw.githubusercontent.com/snap-stanford/Biomni/main/biomni/tool/systems_biology.py ; https://www.biorxiv.org/content/10.64898/2026.05.12.724604v1.full |
| **Comp2GPR** (Castillo, VTT; bioRxiv Jun 2026) | Complex Portal + sequence evidence + "AI-assisted filtering" → Boolean GPRs; tested on Yeast9 | Improved gene-essentiality prediction | [PP] https://www.biorxiv.org/content/10.64898/2026.06.24.734174v1.full.pdf |
| **gapsmith-db** (KAUST BORG) | Licence-clean pathway DB; LLM proposals (RAG over Europe PMC) treated as untrusted and passed through 9 symbolic checks (atom/charge balance, EC validity, ΔG, ATP-cycle, pathway flux) | v0.1 fixtures only | [GH] https://github.com/bio-ontology-research-group/gapsmith-db |
| **KinModGPT** (Maeda & Kurata, IJMS 2023) | GPT→Antimony→SBML kinetic models | Direct GPT→SBML: "all generated SBML models were invalid"; via Antimony all 4 tests passed | [PR] https://www.mdpi.com/1422-0067/24/8/7296 |
| GPT metabolic-network extraction (Phongwattana & Chan, 2023) | GPT mines metabolites/enzymes/relations from abstracts → NetworkX graphs; expert-checked | Qualitative | [PP] https://www.biorxiv.org/content/10.1101/2023.06.27.546560v1.full |
| **MetaboliteChat** (NYU, Nov 2025) / **MetaboT** (arXiv 2510.01724) | Multimodal metabolite LLM / metabolomics-KG multi-agent | details unverified | [PP] https://www.biorxiv.org/content/10.1101/2025.11.07.687008v1.full ; https://arxiv.org/abs/2510.01724 |
| **Self-driven biological discovery** (King lab, Chalmers; bioRxiv 2025) | LLM agents + lab automation + logic scaffold; yeast findings (glutamate/spermine synergy) | Wet-lab validated; no FBA | [PP] https://www.biorxiv.org/content/10.1101/2025.06.24.661378v3.full |
| **ChatGEM** (PNNL; bioRxiv Jul 2026) | Multi-agent (ADEPT) + RAG over COBRApy for simulation/strain design | RAG raised code-generation score 2.63→4.20; experimentally validated succinate strains in *P. putida* | [PP] https://www.biorxiv.org/content/10.64898/2026.07.20.739662v2 |

Not found (after ~40 searches): any verified LLM work on biomass composition or media formulation for GEMs; any LLM/AI publication from the Thiele lab as of the lab's 2024–25 list (https://www.digitalmetabolictwin.org/msp-publications) — note the ERC PoC ChatRD award (Jun 2026) in brief 04; no LLM/agent repos in the openCOBRA org (https://api.github.com/search/repositories?q=org:opencobra).

## (B) Enzyme / transporter / kinetics models

**EC prediction.** CLEAN (Science 2023): New-392 precision 0.596, recall 0.479, F1 0.497 (https://github.com/tttianhao/CLEAN). EC-Bench (Davoudi, Henry et al., Bioinformatics Advances 2026; 10 tools, leakage-controlled splits): at EC level 4 BLASTp F1 0.502 vs EnzBert 0.475; CLEAN most robust OOD (https://academic.oup.com/bioinformaticsadvances/article/6/1/vbag004/8417619). HIT-EC (Nat Commun Jan 2026): micro-F1 0.93 vs CLEAN 0.88, macro-F1 0.84 vs 0.80, evidential uncertainty (https://www.nature.com/articles/s41467-026-68727-3). PLM+MLP benchmark (bioRxiv Apr 2026): 97.0% EC4 accuracy; +31.8 pp over BLAST for *Giardia* (https://www.biorxiv.org/content/10.64898/2026.03.31.715487v1.full). **LLMs alone are poor at EC**: EC-Reason-Bench (arXiv Jul 2026) — closed-book uniformly low, open-book rises sharply, best LLM ≈ neighbour voting; reasoning acts as "arbiter of conflicting neighbors" (https://arxiv.org/abs/2607.26397). Enzyme↔reaction retrieval: ReactZyme (NeurIPS 2024, https://arxiv.org/abs/2408.13659); EnzymeCAGE (~1M pairs; de-orphaning, pathway reconstruction; https://www.biorxiv.org/content/10.1101/2024.12.15.628585v1.full).

**Kinetics.** DLKcat (Nat Catal 2022): r = 0.71 test; 343 fungal ecGEMs; over-predicts activity on random substrates (https://www.nature.com/articles/s41929-022-00798-z). TurNuP (Nat Commun 2023): R² 0.44; DLKcat R² 0.02 for <40% identity; ecGEM proteome MSE −18% (https://www.nature.com/articles/s41467-023-39840-4). CatPred (Nat Commun 2025): kcat R² 0.607 held-out / 0.390 OOD, Km 0.648/0.536, with uncertainty (https://www.nature.com/articles/s41467-025-57215-9). Zheng 2026 (Research Square): all 5 SOTA models "failed to predict absolute kcat"; ranking collapses OOD; 98 xylanase mutants r < 0.3 (https://www.researchsquare.com/article/rs-9154046/v1). GotEnzymes2: 59.6M predictions, 7.3M enzymes, 10,765 species (https://github.com/LiLabTsinghua/GotEnzymes2).

**Transporters.** SPOT (PLOS Biol 2024): 92.4% accuracy, MCC 0.80, AUC 0.96; ~84% below 40% identity; not yet integrated into reconstruction pipelines (https://journals.plos.org/plosbiology/article?id=10.1371%2Fjournal.pbio.3002807).

**ML gap-filling / surrogates.** DNNGIOR (iScience Dec 2024; NN-guided gapfilling for ModelSEED/BiGG; https://github.com/MGXlab/DNNGIOR — numbers unverified). AMN (Nat Commun 2023): differentiable FBA; Q² 0.81 on 17,400 knockout growth rates where FBA ≈ 0 (https://www.nature.com/articles/s41467-023-40380-0). dAMN (Bioinformatics 2026): median R² ≈ 0.98 on E. coli growth curves (https://academic.oup.com/bioinformatics/article/42/5/btag230/8674704). FlowGAT (npj Syst Biol Appl 2024): PRAUC 0.85–0.90 for essentiality (https://www.nature.com/articles/s41540-024-00348-2).

## (C) Failure modes and missing benchmarks

- **Identifier/stoichiometry errors**: LLMs assume wrong identifiers and miss injected stoichiometric mistakes in chat mode (NUS eval, above); direct SBML generation was 100% invalid (KinModGPT).
- **Unverifiable narratives** from plain LLMs vs auditable tool chains (MechAInistic).
- **Chemistry reasoning**: BioMol-LLM-Bench (Bioinformatics 2026; 13 models, 26 tasks): "none … achieves meaningful performance on challenging regression tasks"; chemically inconsistent reasoning; tools improve numerical stability (https://academic.oup.com/bioinformatics/article/42/8/btag550/8740988). ChemBench (Nat Chem 2025): overconfidence, 22% on NMR-signal counting, knowledge deficits (https://www.nature.com/articles/s41557-025-01815-x).
- **Agentic analysis**: BixBench — GPT-4o/Claude 3.5 17% open-answer, MCQ ≈ random (https://arxiv.org/abs/2503.00096); BiomniBench — weak method selection/interpretation.
- **Knowledge vs retrieval**: EC-Reason-Bench shows LLM parametric knowledge of EC is insufficient.
- **Kinetics OOD collapse** (TurNuP, CatPred, Zheng 2026).
- **Missing**: no public benchmark for LLM-driven reconstruction, gap-filling candidate ranking, GPR inference, ID mapping, or model repair. Deterministic checkers exist and should gate LLM proposals: MACAW (Genome Biology 2025: dead-end, dilution, diphosphate, duplicate, loop tests; https://github.com/Devlin-Moyer/macaw ; DOI 10.1186/s13059-025-03533-6), memote (https://github.com/opencobra/memote), gapsmith's nine checks.

## (D) Adjacent AI-for-science frameworks

- **AI co-scientist** (Google, Feb 2025): 6 agents + tournament Elo; AML repurposing, liver-fibrosis targets (p < 0.01 in organoids), cf-PICI rediscovery (https://research.google/blog/accelerating-scientific-breakthroughs-with-an-ai-co-scientist/).
- **Robin** (FutureHouse; Nature May 2026): ripasudil for dry AMD via iterative hypothesis→experiment→RNA-seq loop (https://www.futurehouse.org/research-announcements/demonstrating-end-to-end-scientific-discovery-with-robin-a-multi-agent-system).
- **Kosmos** (Edison, Nov 2025): 1,500 papers + 42,000 lines of code per run; 79.4% statement accuracy by experts; reproduced a metabolomics finding (nucleotide metabolism in hypothermic brain) (https://edisonscientific.com/news/announcing-kosmos).
- **Coscientist** (Nature 2023) and **ChemCrow** (Nat Mach Intell 2024; 18 tools) — tool-grounded chemistry agents (https://www.nature.com/articles/s42256-024-00832-8.pdf).
- **BioChatter** (Nat Biotech 2025): KG-grounded LLM platform with benchmarking (https://www.nature.com/articles/s41587-024-02534-3).
- Lessons for GEMs: (i) ground every claim in an executable tool call (COBRApy, MACAW); (ii) reviewer/critic agents and tournament ranking; (iii) expert rubrics + process-level grading (BiomniBench); (iv) closed loops with wet-lab validation (Robin, King lab); (v) persistent structured world-models across long tasks (Kosmos).

## (E) Prioritised opportunities for a frontier model now

1. **GEM-Bench** (highest priority; unlocks everything else). Tasks: reaction balancing/charge, cross-DB ID mapping (BiGG↔MetaNetX↔ModelSEED), GPR inference vs Yeast9/iML1515, held-out-reaction gap-fill ranking, MACAW/memote defect repair, COBRApy coding tasks; automatic grading. *Risk*: leakage from training data; mitigate with post-cutoff curations. *Success*: public leaderboard; frontier model+tools ≥ 90% on balancing/ID tasks, reported error rates per task.
2. **Propose-then-verify curation agent.** LLM proposes reactions/GPR/localisation with citations; deterministic gate (mass/charge, ΔG, MACAW, memote, growth-phenotype regression) accepts; opens PRs to Human-GEM/Yeast-GEM/AGORA2/APOLLO (gapsmith-db and Human-GEM CI are templates). *Risk*: subtle ID hallucinations passing checks. *Success*: ≥100 maintainer-merged fixes with >50% PR acceptance; zero new MACAW failures.
3. **Literature-evidence extraction at scale** (D2Cell/BioPromptX/MetaKnogic style) for orphan reactions, GPRs, kcat/Km, media. *Success*: ≥90% precision on an expert-annotated held-out set; ≥20% reduction of unannotated reactions in a target GEM.
4. **Annotation ensemble in draft reconstruction**: CLEAN/HIT-EC/EnzymeCAGE + SPOT + LLM as evidence-arbiter (per EC-Reason-Bench) inside CarveMe/ModelSEED/gapseq. *Risk*: OOD errors compounding. *Success*: improved Biolog growth/no-growth MCC on ≥50 strains vs baseline pipelines.
5. **Audited NL interface to COBRA** (MCP + Architect/Reviewer). *Success*: ≥90% expert-graded correctness on 100 standard analyses (FBA/FVA/knockouts/context-specific models) with reproducible notebooks; adoption by openCOBRA.
6. **Kinetic parameterisation with uncertainty** (CatPred/GotEnzymes2 + extracted measurements) for ecGEMs. *Success*: proteome-allocation MSE gains ≥ TurNuP's 18% on new species.
7. **Mechanistically grounded hypothesis generation** (MechAInistic/Kosmos pattern) with pre-registered predictions tested in wet lab. *Success*: ≥3 validated novel predictions (essential genes, drug targets, cross-feeding).
8. **Error-detection red-teaming**: measure and fix LLM sensitivity to injected stoichiometric/ID errors (NUS finding). *Success*: detection rate from chat-level baseline to >95% with tool-augmented prompts.
9. **Media/transporter/biomass formulation** — the empty niche: SPOT + literature extraction → exchange reactions and defined media. *Success*: AUC > 0.85 on MediaDive-style growth data.

## (F) Key sources

1. Human2, PNAS 2026 — https://www.pnas.org/doi/10.1073/pnas.2516511123
2. Human-GEM releases (GPT-3 gene pruning, MACAW CI) — https://github.com/SysBioChalmers/Human-GEM/releases
3. MechAInistic, bioRxiv/arXiv 2026 — https://arxiv.org/abs/2607.18249
4. LLM evaluation for GEMs, bioRxiv 2026 — https://www.biorxiv.org/content/10.64898/2026.06.03.730004v1
5. D2Cell, Trends Biotechnol 2026 — https://www.cell.com/trends/biotechnology/abstract/S0167-7799(26)00135-6
6. MetaKnogic-Alpha, bioRxiv 2026 — https://www.biorxiv.org/content/10.64898/2026.02.05.704050v1.full
7. gem-flux-mcp / ModelSEEDagent — https://github.com/jplfaria/gem-flux-mcp
8. Biomni, bioRxiv 2025 — https://www.biorxiv.org/content/10.1101/2025.05.30.656746v1
9. BiomniBench, 2026 — https://www.biorxiv.org/content/10.64898/2026.05.12.724604v1.full
10. BixBench, 2025 — https://arxiv.org/abs/2503.00096
11. EC-Bench, Bioinf Adv 2026 — https://academic.oup.com/bioinformaticsadvances/article/6/1/vbag004/8417619
12. EC-Reason-Bench, 2026 — https://arxiv.org/abs/2607.26397
13. HIT-EC, Nat Commun 2026 — https://www.nature.com/articles/s41467-026-68727-3
14. CatPred, Nat Commun 2025 — https://www.nature.com/articles/s41467-025-57215-9
15. TurNuP, Nat Commun 2023 — https://www.nature.com/articles/s41467-023-39840-4
16. kcat "what works…", 2026 — https://www.researchsquare.com/article/rs-9154046/v1
17. SPOT, PLOS Biol 2024 — https://journals.plos.org/plosbiology/article?id=10.1371%2Fjournal.pbio.3002807
18. AMN, Nat Commun 2023 — https://www.nature.com/articles/s41467-023-40380-0
19. MACAW, Genome Biol 2025 — https://github.com/Devlin-Moyer/macaw
20. AI co-scientist, 2025 — https://research.google/blog/accelerating-scientific-breakthroughs-with-an-ai-co-scientist/ ; Robin, Nature 2026 — https://www.futurehouse.org/research-announcements/demonstrating-end-to-end-scientific-discovery-with-robin-a-multi-agent-system ; Kosmos — https://edisonscientific.com/news/announcing-kosmos

**Caveats**: PNAS/PMC/bioRxiv full texts were blocked for this brief, so Human2's LLM-usage details, MechAInistic's model versions, and the NUS benchmark's per-model scores are reported only at abstract level (marked unverified). KEGG/BiGG/MetaCyc MCP servers and LLM-based media/biomass work were not found.
