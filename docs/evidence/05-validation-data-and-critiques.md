# Evidence brief 05 — Validation data, benchmarks and field-level critiques (state as of September 2026)

Compiled 4 Sept 2026 by Claude research agents for the "Metabolic modelling improvements" project. ~20 searches and ~60 primary-source fetches. Items marked *unverified* could not be confirmed in-session (PMC, PLoS, bioRxiv and DepMap throttled or captcha-blocked) and should be checked before being cited. This brief feeds the Landscape & Roadmap document (00-landscape-and-roadmap.md).

## (A) Ground-truth datasets

| Dataset | Measures | Organisms | Size | Machine-readable? | URL |
|---|---|---|---|---|---|
| Keio collection (Baba 2006) | Single-gene deletion viability, LB + MOPS-glucose | *E. coli* K-12 | 3,985 mutants of 4,288 targets; 303 essential | Tables via GenoBase | [msb4100050](https://link.springer.com/article/10.1038/msb4100050) |
| Fitness Browser (RB-TnSeq) | Conditional gene fitness across carbon/N sources, stresses | 50 bacteria | 7,593 genome-wide experiments; 2018 paper: 32 bacteria, 11,779 genes given phenotypes | Yes: per-organism TSV + sqlite (~5 GB); DOE-ENIGMA funded | [orgAll](https://fit.genomics.lbl.gov/cgi-bin/orgAll.cgi), [Price 2018](https://www.nature.com/articles/s41586-018-0124-0) |
| DEG 15 | Essential genes (prokaryotes + eukaryotes) | Multi-kingdom | Counts **unverified** (prokaryote records "more than doubled" vs DEG 10) | Download | [DEG 15](https://academic.oup.com/nar/article/49/D1/D677/5937083) |
| Nichols 2011 chemical genomics | Mutant × condition growth landscape | *E. coli* | ~3,979 mutants × 324 conditions (**unverified** numbers) | Supplementary tables | [Cell 2011](https://doi.org/10.1016/j.cell.2010.11.052) |
| DepMap / Broad–Sanger CRISPR | Cancer-cell-line gene dependency | Human | >900 lines in unified 2020 release; current larger (**unverified**) | CSV | [depmap.org/broad-sanger](https://depmap.org/broad-sanger/) |
| Biolog PM (per-study) | Substrate use (C/N/P/S) | Varies; e.g. *K. pneumoniae* KPPR1 190 substrates | Per-study | Supplementary tables, no central repository | [Bactabolize](https://pmc.ncbi.nlm.nih.gov/articles/PMC10564454/) |
| BacDive 2025 | Strain phenotypes, API tests, enzyme, growth conditions | Prokaryotes | 97,334 strains (2024.1; >100k Dec 2025); 1.59 M API data points for 24,112 strains | Yes: CSV/JSON/REST/SPARQL; CC BY 4.0 | [NAR 2025](https://academic.oup.com/nar/article/53/D1/D748/7848838), [news](https://bacdive.dsmz.de/news/40) |
| ProTraits (as used by gapseq) | Carbon-source utilisation | 526 organisms | 1,795 high-confidence phenotypes / 48 C-sources | Yes | [gapseq 2021](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02295-1) |
| KOMODO | Organism–medium pairs | 18,049 strains | 3,335 media; 20,824 pairings; 1,324 components | Web/GROWREC | [Oberhardt 2015](https://www.nature.com/articles/ncomms9493) |
| CeCaFDB | 13C-MFA central-carbon flux distributions | 36 organisms | 581 flux distributions, 118 refs; last major update 2015 | Excel download | [NAR 2015](https://academic.oup.com/nar/article/43/D1/D549/2438317) |
| NCI-60 exometabolomics + Achilles (benchmark bundle) | Uptake/secretion, growth, CERES essentiality | 22–60 human lines | Bundled in Jamialahmadi 2019 | Yes (supplement) | [PLoS CB 2019](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006936) |
| AGORA2 validation sets (NJC19, BacDive, Madin; 253 drug–microbe pairs) | Uptake/secretion, drug transformation | 7,302 strains modelled | 253 drug–microbe pairs | Supplement | [Heinken 2023](https://www.nature.com/articles/s41587-022-01628-0) |
| PaxDb v6.0 | Protein abundance | ~390 species | 1,639 datasets; 1.38 M proteins; LLM-assisted curation | TSV | [NAR 2026](https://academic.oup.com/nar/article/54/D1/D427/8313448) |
| BRENDA 2026.1 | Km, kcat, etc. | All | 8,831 EC classes; 186,457 Km; 94,073 kcat; CC BY 4.0 | Yes | [statistics](https://www.brenda-enzymes.org/statistics.php) |
| GotEnzymes2 (2025) | Predicted kcat/Km/thermal | 10,765 species | 59.6 M predicted entries | Yes | [abstract](https://tagteam.harvard.edu/hub_feeds/4353/feed_items/17114468) |
| SABIO-RK | Curated reaction kinetics | All | Size **unverified** | Yes (web services) | [bio.tools](https://bio.tools/sabio-rk) |
| eQuilibrator 3.0 | ΔG′ via component contribution | — | ~1 M compounds (MetaNetX); Python API; MIT/CC | Yes | [PubMed](https://pubmed.ncbi.nlm.nih.gov/34850162/) |
| HMDB 5.0 | Human metabolites + concentrations | Human | 217,920 compounds; >19,715 concentration values added in 5.0; CC BY-NC | Yes (XML/JSON/CSV) | [NAR 2022](https://academic.oup.com/nar/article/50/D1/D622/6431815) |
| MetaboLights | Raw metabolomics studies | 6,815 organism/parts | 8,544 studies, 1,358 public (Sept 2023) | ISA-Tab | [NAR 2024](https://academic.oup.com/nar/article/52/D1/D640/7424432) |

Gaps: no central Biolog/growth-rate compendium; no maintained fluxomics compendium after 2015; exometabolomics scattered in supplements.

## (B) Benchmark efforts and reported accuracies

- **MEMOTE (2020)** tests format, annotation, stoichiometric consistency, biomass, cycles; experimental growth/essentiality tests exist but need user data. It is a *structural* QC score, not predictive accuracy ([Lieven 2020](https://www.nature.com/articles/s41587-020-0446-y), [docs](https://memote.readthedocs.io/en/latest/)). Yeast9 reports "+27% MEMOTE score" with only "moderately improved" essentiality ([Yeast9](https://link.springer.com/article/10.1038/s44320-024-00060-7)).
- **Tool head-to-heads**: Mendoza 2019 (7 tools, 3 bacteria) — max reaction coverage vs curated models only 72% (*L. plantarum*) / 58% (*B. pertussis*); "none of the tools outperforms the others" ([Genome Biol](https://link.springer.com/article/10.1186/s13059-019-1769-1)). gapseq 2021 — carbon-source TP 47% vs 31% (ModelSEED) vs 24% (CarveMe); SCFA products 91% vs 12% vs 0%; essentiality on 5 species, curated models generally better ([gapseq](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02295-1)). Bactabolize 2023 — *K. pneumoniae* substrate accuracy 85.8% (Bactabolize) vs 77.9% (CarveMe) vs 77.0% (gapseq); TraDIS essentiality 80.6–84.0% across tools ([Bactabolize](https://pmc.ncbi.nlm.nih.gov/articles/PMC10564454/)). AGORA2 — drug-transformation accuracy 0.81; uptake/secretion 0.81–0.84; enzyme activity 0.72 ([AGORA2](https://www.nature.com/articles/s41587-022-01628-0)).
- **Fitness-data benchmark (best current practice)**: Bernstein et al. 2023 evaluated four *E. coli* GEM generations against RB-TnSeq over 25 carbon sources; recommend AUC-PR; overall accuracy ~93.8% at fitness<−2 but AUC-PR much lower; the published version reports iML1515 AUC-PR ≈ 0.81 rising to ≈ 0.85 after fixes (a separate fetch of the published article confirmed these values; the preprint's lower figures are superseded); 21 false-negative cofactor-biosynthesis genes, 8 isoenzyme GPR errors fixed ([MSB 2023](https://link.springer.com/article/10.15252/msb.202311566), [preprint](https://www.biorxiv.org/content/10.1101/2023.01.05.522875v1.full)).
- **Growth rates**: Yeast9 R² = 0.842 (aerobic/anaerobic); Pearson 0.66 (single KO) / 0.78 (double KO); synthetic lethality ~80% accuracy but precision 10.0%, recall 13.1% ([Yeast9](https://link.springer.com/article/10.1038/s44320-024-00060-7)).
- **Fluxes vs 13C-MFA**: Schuetz 2007 — 99 objective/constraint combos, 6 conditions; "no single objective describes the flux states under all conditions" ([MSB](https://link.springer.com/article/10.1038/msb4100162)). Plant leaf: pFBA weighted error 94–180% vs 9–22% with expression constraints ([Bioinformatics 2023](https://academic.oup.com/bioinformatics/article/39/5/btad186/7114031)). Kaste & Shachar-Hill 2024 argue validation/model-selection practice is immature and propose χ²-complementary frameworks ([Biotech Prog](https://api.crossref.org/works/10.1002/btpr.3413)).
- **Human biomarkers**: Shlomi 2009 precision 0.73–0.76, recall 0.40–0.56 ([MSB](https://link.springer.com/article/10.1038/msb.2009.22)); whole-body models 85.3%/84.9% correct on 252 biomarkers over 57 IEMs vs Recon3D 50.2% ([Thiele 2020](https://link.springer.com/article/10.15252/msb.20198982)).
- **Context-specific extraction**: 9 methods; only 5 gave significant uptake/secretion correlations; TRFBA best for growth ([Jamialahmadi 2019](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006936)).
- **Solvers**: 6 solvers, 7 models — identical objectives ([mSystems 2024](https://journals.asm.org/doi/10.1128/msystems.00833-23)); earlier discrepancies traced to SBML parsing, 9/88 models needed encoding fixes ([Ebrahim 2015](https://link.springer.com/article/10.15252/msb.20156157)).
- **Community challenges**: no DREAM challenge on flux/GEMs ([DREAM collection](https://collections.plos.org/collection/dream/)); no 2024–26 standard benchmark found — only calls for "a robust validation framework" ([Quinn-Bohmann 2025](https://www.nature.com/articles/s41564-025-01972-2)).

## (C) Critiques and failure modes

**Well-evidenced**
- Biomass/objective: swapping biomass reactions changes predicted essential reactions by up to 30% ([Bernstein 2021](https://link.springer.com/article/10.1186/s13059-021-02289-z)); no universal objective ([Schuetz 2007](https://link.springer.com/article/10.1038/msb4100162)); biomass coefficients "considerably affect the attainable fluxes" ([Simensen 2022](https://doaj.org/article/e952e80f4d5e4313a67323368887da63)).
- Alternate optima/degeneracy: foundational ([Mahadevan & Schilling 2003](https://api.crossref.org/works/10.1016/j.ymben.2003.09.002)); Bernstein 2021 lists optimal-solution degeneracy as core uncertainty.
- Missing regulation: false positives where machinery is "present but not active" ([Bactabolize](https://pmc.ncbi.nlm.nih.gov/articles/PMC10564454/)); isoenzyme GPRs without regulation create "overly promiscuous" networks ([Bernstein 2023](https://link.springer.com/article/10.15252/msb.202311566)).
- Annotation error propagation: ~35% of *E. coli* genes lack experimental function evidence ([Fang 2020](https://www.nature.com/articles/s41579-020-00440-4)); database misannotation documented by [Schnoes 2009](https://api.crossref.org/works/10.1371/journal.pcbi.1000605) (range 5–63% **unverified**); ~50% of reactions automated tools miss lack gene associations in curated models ([Mendoza 2019](https://link.springer.com/article/10.1186/s13059-019-1769-1)).
- Evaluation fragility: apparent accuracy depends on assumed medium (cofactor leakage) and metric; naive metrics suggested accuracy *declining* with model size until corrected ([Bernstein 2023](https://www.biorxiv.org/content/10.1101/2023.01.05.522875v1.full)); "the small amount of experimental data available is likely to introduce biases" ([Carter 2023](https://academic.oup.com/bib/article/25/1/bbad439/7457942)).
- Reproducibility: encoding/parsing errors ([Ebrahim 2015](https://link.springer.com/article/10.15252/msb.20156157)); in adjacent ODE modelling only 51% of 455 models reproduced directly ([Tiwari 2021](https://link.springer.com/article/10.15252/msb.20209982)).

**Opinion / weakly evidenced**
- "Predicts everything and nothing": no citable source found for the phrase; substantively covered by degeneracy above.
- Overfitting curated models to the same essentiality/Biolog data used for validation: no direct study found — a real structural risk but **unverified** empirically.
- Clinical/industrial disconnect: industry adoption "in its infancy", blocked by skills, workflows, software ([Richelle 2020](https://www.nature.com/articles/s41540-020-0127-y)); gut community modelling lacks validation frameworks ([Quinn-Bohmann 2025](https://www.nature.com/articles/s41564-025-01972-2)).

## (D) Field-level signals

- Scale: 6,239 GEMs by Feb 2019, only 183 manually reconstructed; CarveMe produced 91.8% of bacterial GEMs ([Gu 2019](https://link.springer.com/article/10.1186/s13059-019-1730-3)); BiGG holds 108 curated GEMs ([Carter 2023](https://academic.oup.com/bib/article/25/1/bbad439/7457942)); AGORA2 7,302 strains; ~247k human-microbe GEMs (2025). Automated volume vastly outpaces validation.
- Training: systems-biology master's programmes struggle to balance biology, programming and maths; "learning programming languages ... is challenging" for biologists ([Sauter 2025](https://erc.bioscientifica.com/view/journals/erc/32/6/ERC-25-0024.xml)); industry needs "broad interdisciplinarity of competences" ([Richelle 2020](https://www.nature.com/articles/s41540-020-0127-y)); teaching modules exist ([Kaste 2023](https://iubmb.onlinelibrary.wiley.com/doi/10.1002/bmb.21777)).
- Funding/infrastructure: key ground truth sits in infrastructure grants (Fitness Browser via DOE ENIGMA; BRENDA/BacDive as ELIXIR/DSMZ resources). Per-year GEM paper counts and conference attendance: **not verified** this session.
- No community challenge ever run ([DREAM](https://collections.plos.org/collection/dream/)).

## (E) Opportunities for frontier-model intelligence

1. **GEM-Bench harness** — one command evaluates any SBML model against Fitness Browser (7,593 experiments), BacDive, ProTraits, Keio, DepMap, CeCaFDB using AUC-PR and explicit media definitions. *Success*: ≥50 organisms, leaderboard reproduced by two independent groups within ±0.02 AUC-PR.
2. **Automated error triage** — agents explain each false prediction and propose literature-backed fixes in Bernstein-2023's three classes (medium cofactors, GPR/isoenzymes, reversibility). *Success*: ≥50% of FN/FP resolved with citations, gains confirmed on held-out conditions.
3. **Leakage-controlled validation** — track provenance of every curation edit; report train/test splits. *Success*: benchmark cards state which ground truth touched curation; held-out AUC-PR published alongside in-sample.
4. **Literature-to-phenotype extraction** — grow ProTraits/BacDive-style tables (gapseq hand-collected only 24 organisms' fermentation products) and IEM biomarker tables. *Success*: ≥10⁴ new records with provenance, ≥95% precision on expert audit.
5. **Fluxomics compendium refresh** — CeCaFDB frozen since 2015; extract 13C-MFA fluxes into a standard schema, run FBA-vs-MFA leaderboard per objective. *Success*: ≥1,000 flux distributions; per-objective error reported.
6. **Agentic GEM task benchmark** — current bio benchmarks (BixBench 17% ([arXiv](https://arxiv.org/abs/2503.00096)), ScienceAgentBench 42.2% ([arXiv](https://arxiv.org/abs/2410.05080)), LAB-Bench ([arXiv](https://arxiv.org/abs/2407.10362)), MetaBench ([arXiv](https://arxiv.org/html/2510.14944v1))) contain no GEM tasks. Seeds exist: LLMs recover 64–91% of *E. coli* core reactions but many networks infeasible ([Tewari 2026](https://elifesciences.org/articles/109709)); conversational prompting "failed to identify injected errors" ([Yeoh 2026](https://api.biorxiv.org/details/biorxiv/10.64898/2026.06.03.730004)). *Success*: injected-error detection ≥90%; de novo reconstructions match curated coverage (>72%) and Biolog accuracy ≥85%.
7. **Reproducibility audit** — re-run every BiGG/BioModels GEM through SBML-fbc validation and two solvers. *Success*: 100% reproduce published growth within 1e-6.
8. **Calibrated ensembles** — Bernstein-2021 ensemble reconstruction with probabilities scored against fitness data. *Success*: expected calibration error <0.05.

## (F) Key sources

1. Evaluating E. coli GEM accuracy with high-throughput mutant fitness data (2023) — https://link.springer.com/article/10.15252/msb.202311566
2. Addressing uncertainty in GEM reconstruction and analysis (2021) — https://link.springer.com/article/10.1186/s13059-021-02289-z
3. MEMOTE (2020) — https://www.nature.com/articles/s41587-020-0446-y
4. Systematic assessment of GEM reconstruction tools (2019) — https://link.springer.com/article/10.1186/s13059-019-1769-1
5. gapseq (2021) — https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02295-1
6. Bactabolize (2023) — https://pmc.ncbi.nlm.nih.gov/articles/PMC10564454/
7. AGORA2 (2023) — https://www.nature.com/articles/s41587-022-01628-0
8. Systematic evaluation of objective functions (2007) — https://link.springer.com/article/10.1038/msb4100162
9. Model validation and selection in MFA/FBA (2024) — https://api.crossref.org/works/10.1002/btpr.3413
10. Mutant phenotypes for thousands of bacterial genes (2018) — https://www.nature.com/articles/s41586-018-0124-0
11. BacDive 2025 — https://academic.oup.com/nar/article/53/D1/D748/7848838
12. Yeast9 (2024) — https://link.springer.com/article/10.1038/s44320-024-00060-7
13. Benchmarking biochemical networks generated by LLMs (2026) — https://elifesciences.org/articles/109709
14. LLM capabilities for GEM analysis (2026) — https://api.biorxiv.org/details/biorxiv/10.64898/2026.06.03.730004
15. Community-scale gut models perspective (2025) — https://www.nature.com/articles/s41564-025-01972-2
