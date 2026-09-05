# Evidence brief 03 — Tooling, infrastructure, standards and reproducibility (state as of September 2026)

Compiled 4 Sept 2026 by Claude research agents for the "Metabolic modelling improvements" project. Facts were verified by direct fetches of PyPI JSON, GitHub pages and Semantic Scholar where publisher pages were blocked. GitHub release pages do not render years; years come from PyPI unless flagged "inferred". Items marked *unverified* should be checked before being cited. This brief feeds the Landscape & Roadmap document (00-landscape-and-roadmap.md).

## (A) Key facts and developments 2023–2026

**COBRApy**
- Latest release **0.32.1 (2026-08-11)**; 0.32.0 (2026-08-10) added a `chrr` uniform sampler via optional `hopsy` and Python 3.14 support; 0.31.1 (2026-03-26); 0.30.0 (2025-10-14, dropped Py3.8); 0.29.1 (2024-09-19); 0.29.0 (2023-10-31, hybrid HiGHS/OSQP solver) — https://pypi.org/project/cobra/ , https://opencobra.github.io/cobrapy/releases/0.32.0/
- Cadence collapsed: 22 releases in 2017 → 1 (2024), 1 (2025), 3 (2026) (PyPI release JSON). A 13-month gap separated 0.29.1 and 0.30.0.
- Every release since 0.26 is cut by one person (Christian Diener); PyPI maintainer-of-record is Moritz Beber — https://github.com/opencobra/cobrapy/releases , https://pypi.org/project/cobra/
- Repo: 581 stars, 238 forks, 77 open issues, 20 open PRs — https://github.com/opencobra/cobrapy . Human PRs sit unmerged for years: FBC-v3 metadata/key-value-pair support (#1440, opened Apr 2025, reviewed positively Oct 2025, still open) https://github.com/opencobra/cobrapy/pull/1440 ; fast gap-filling (#1450, Jun 2025, awaiting rebase Mar 2026) https://github.com/opencobra/cobrapy/pull/1450 ; loopless "potentials" constraints (#1472, Apr 2026); PRs from 2022 still open — https://github.com/opencobra/cobrapy/pulls
- Hard dependency on `optlang~=1.9` and `python-libsbml~=5.19`; optlang 1.9.1 (2026-06-02), 1.9.0 added Py3.14, 1.8.0 (2023) introduced the matrix interface and HiGHS+OSQP hybrid — https://pypi.org/project/optlang/ , https://github.com/opencobra/optlang/blob/1.9.1/CHANGELOG.rst
- Defaults: solver preference Gurobi→CPLEX→GLPK; bounds ±1000 — https://cobrapy.readthedocs.io/en/latest/configuration.html
- No governance/roadmap document found; user-defined-constraints API discussion (2022) unresolved — https://github.com/opencobra/cobrapy/discussions/1240

**COBRA Toolbox (MATLAB)**
- Releases v3.5 (2 Jul; Persephone, nutrition toolbox, Fast-SL, CHRR barrier rounding, entropic FBA, 22 new contributors), v3.6 (20 Oct; MicroMap), v3.7 "Persephone Pipeline Publication Release" (8 May; conserved-moiety framework "Rahou et al., JTB 2025", renewed QP interface, 13 contributors) — https://github.com/opencobra/cobratoolbox/releases , https://github.com/opencobra/cobratoolbox/releases/tag/v3.7 . Years are not rendered; v3.7 inferred 2026 (unverified). Releases since v3.5 are made by @farid-zare, earlier ones by @rmtfleming.
- Repo: 286 stars, 336 forks, 65 open issues, **0 open PRs / 2,352 closed**, 10,365 commits, GPLv3 — https://github.com/opencobra/cobratoolbox , https://github.com/opencobra/cobratoolbox/pulls
- Development plan targets "v4.0", lead developer Ronan Fleming (Galway); no mention of Python — https://opencobra.github.io/cobratoolbox/stable/plan.html
- MATLAB-only capabilities: 35+ analysis sub-modules incl. whole-body/PSCM (`optimizeWBModel` with scaling/numeric-emphasis flags for MOSEK/CPLEX, `checkIEM_WBM`), Persephone (metagenomics→personalised host–microbiome WBMs via SeqC/MARS/mgPipe), DEMETER, rBioNet, XomicsToModel, thermoKernel, MetaboAnnotator, nutrition algorithm, EFMviz/Paint4Net; ~115 tutorials — https://opencobra.github.io/cobratoolbox/stable/modules/index.html , https://opencobra.github.io/cobratoolbox/stable/modules/analysis/wholeBody/PSCMToolbox/index.html , https://opencobra.github.io/cobratoolbox/stable/modules/analysis/persephone/index.html , https://opencobra.github.io/cobratoolbox/stable/tutorials/
- MATLAB→Python: no published migration plan found. The 2017 parity issue (#533) was closed under a "1.0 milestone" without a comparison table — https://github.com/opencobra/cobrapy/issues/533

**Ecosystem momentum (PyPI dates)**
- Gaining: straindesign 1.19.1 (2026-07-26), MICOM 0.39.1 (2026-08-11, requires cobra≥0.32, highspy), CNApy 1.2.8 (2026-04-02), coralME 1.2.3.1 (2026-04-10), COMETSpy 0.6.3 (2026-02-18), MEWpy 1.0.0 (2026-01-04), cobramod 1.3.1 (2026-03-11), hopsy 1.7.0, PolyRound 0.5.0 (2026-08-24), gapseq v2.1.0 (2026-05-30, R) — https://pypi.org/project/straindesign/ , https://pypi.org/project/micom/ , https://github.com/jotech/gapseq
- Stalling: MEMOTE 0.17.0 (2024-01-08; 134 open issues) https://github.com/opencobra/memote ; Escher 1.8.1 (2024-10-25; canonical repo moved to opencobra/escher with 6 stars) https://github.com/zakandrewking/escher ; reframed 1.6.0 (2025-04-28); CarveMe 1.6.6 (2025-09-12); ModelSEEDpy 0.4.2 (2025-03-06).
- Dead: cameo 0.13.6 (2021-11-09), pytfa 0.9.4 (2022-06-27), dfba 0.1.8 (2020), geckopy 2.0.2 (2021), mackinac (2017) — https://pypi.org/project/cameo/ , https://pypi.org/project/pytfa/
- R: sybil last published 2021-05-31 and described as "discontinued" by its successor; gapseq now depends on `cobrar` (GitHub-only, libsbml-backed, GLPK/CPLEX) — https://rdrr.io/cran/sybil/ , https://github.com/Waschina/cobrar , https://gapseq.readthedocs.io/en/latest/install.html
- Julia: COBREXA 2 (Bioinformatics, Feb 2025; ConstraintTrees-based, 3 citations on S2); repo has 19 stars (legacy repo 46) — https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaf056/8005852 , https://github.com/COBREXA/COBREXA.jl

**Solvers**
- Machado's benchmark (mSystems 2023/24): Gurobi fastest; SCIP/HiGHS/GLPK similar on LP; CPLEX 100× slowdown on >4-member communities; GLPK/CBC failed largest MILP within a week; all solvers agreed on objectives — https://journals.asm.org/doi/10.1128/msystems.00833-23
- HiGHS 1.15.1 (2026-07-02, MIT); 1.10.0 added GPU cuPDLP-C; release notes warn PDLP "cannot be used" where a basic solution is required — https://github.com/ERGO-Code/HiGHS/releases , https://pypi.org/project/highspy/
- cuPDLPx (Lu, Peng, Yang 2025) benchmarks at 1e-4/1e-8 KKT tolerance; no GEM-specific benchmark found — https://arxiv.org/html/2507.14051
- Gurobi academic: named-user 1-year; WLS 2 concurrent sessions, 90-day renewable, permanent internet connection required (page updated 2026-03-24) — https://support.gurobi.com/hc/en-us/articles/34672988479633-What-are-the-restrictions-on-using-an-academic-WLS-license ; gurobipy 13.0.3 (2026-08-25), cplex 22.2.0.1 (2026-07-15) on PyPI.

**Standards & repositories**
- SBML FBC v3 (quadratic objectives, user-defined constraints, key-value pairs, double charges) proposed 2019; implemented in libSBML 5.20.x and cbmpy 0.8.8 ("alpha", Nov 2023); still not in a COBRApy release; sbml.org's fbc page lists only v1/v2 — https://sourceforge.net/p/sbml/mailman/message/36712110/ , https://github.com/sbmlteam/libsbml/releases , https://github.com/SystemsBioinformatics/cbmpy/releases , https://sbml.org/documents/specifications/level-3/version-1/fbc/
- BiGG: legacy site states data "not updated since October 2019"; SBRG GitHub v1.6 added 23 models; new **bigg.bio** (DTU BRIGHT + UCSD) reports 4,881 models, 7,322 reactions, ChEBI/InChI/RHEA-grounded IDs, non-commercial licence — http://bigg.ucsd.edu/ , https://github.com/SBRG/bigg_models/releases , https://bigg.bio/about
- MetaNetX MNXref 4.5 (Nov 2025): 1.50M compounds, 83.8K reactions, VMH newly integrated, SPARQL endpoint with LLM "ExpasyGPT" examples; authors admit cross-refs "may be contradictory, outdated or missing" — https://academic.oup.com/nar/article/54/D1/D617/8343507
- ModelSEED database has no GitHub releases (70 stars) — https://github.com/ModelSEED/ModelSEEDDatabase/releases
- BacDive 2024.1: 97,334 strains, links to MediaDive; no model-media conversion — https://academic.oup.com/nar/article/53/D1/D748/7848838
- VMH homepage returned no parsable content (unverified status); AGORA2/APOLLO (247,092 reconstructions) hosted via GitHub/Cell Systems — https://github.com/VirtualMetabolicHuman/AGORA2 , https://www.cell.com/cell-systems/fulltext/S2405-4712(25)00029-8

## (B) Tool table

| Tool | Lang | Latest version/date | Maintainer activity | Unique capabilities | Notable gaps |
|---|---|---|---|---|---|
| COBRApy | Python | 0.32.1, 2026-08-11 | 1 release engineer; 20 open PRs, 77 issues | optlang abstraction, hybrid HiGHS/OSQP, chrr sampling | No FBC v3, no thermo/WBM/community tooling, PR backlog |
| COBRA Toolbox | MATLAB | v3.7, 8 May (year unverified) | Fleming lab; 0 open PRs | WBM/PSCM, Persephone, DEMETER, rBioNet, XomicsToModel, thermoKernel | MATLAB licence, no Python path |
| optlang | Python | 1.9.1, 2026-06-02 | Sonnenschein/Diener | solver-agnostic LP/QP/MILP | sympy overhead |
| straindesign | Python | 1.19.1, 2026-07-26 | Klamt lab, active | MCS, OptKnock, GPR compression, SCIP | small community (52★) |
| MICOM | Python | 0.39.1, 2026-08-11 | Diener lab, active | community models, cooperative tradeoff | gut-centric |
| MEMOTE | Python | 0.17.0, 2024-01-08 | stalled; 134 issues | standard QC tests | no FBC3, aging deps |
| Escher | JS/Py | 1.8.1, 2024-10-25 | low | pathway maps | repo split, low activity |
| COBREXA.jl | Julia | v2 (paper Feb 2025) | LCSB/HHU | HPC, ConstraintTrees composability | 19★, tiny community |
| cobrar (sybil successor) | R | GitHub only | Waschina | gapseq backend | not on CRAN |
| gapseq | R | v2.1.0, 2026-05-30 | Kiel, active | pathway-aware reconstruction | R/GLPK stack |
| cameo / pytfa | Python | 2021 / 2022 | none | strain design / TFA | abandoned |
| cbmpy | Python | 0.8.8.1, Nov 2023 | Olivier | only FBC v3 writer | low adoption |

## (C) Reproducibility & interoperability failure modes

1. Same model, same solver (Gurobi), different tool: 355 "active" reactions in COBRApy vs 1,143 in the Toolbox; attributed to solver defaults and alternative optima; unresolved — https://groups.google.com/g/cobra-toolbox/c/7ThHMssO_s4
2. Loopless methods non-deterministic: `loopless_solution` gave fluxes ranging 0–5 across runs (#1262, fixed via PR #1341); loopless FVA changed several-fold between 0.9.0 and 0.13.3 (#738, closed without explanation) — https://github.com/opencobra/cobrapy/issues/1262 , https://github.com/opencobra/cobrapy/issues/738
3. Solver-state leakage: FVA on a copied model differs by 1e-3 with CPLEX but not GLPK (#924, open) — https://github.com/opencobra/cobrapy/issues/924
4. FROG study of 65 BioModels submissions: ~40% reproduced unaided, 28% needed technical fixes, 32% needed authors; causes were metadata gaps, invalid SBML, solver numerical precision, default bounds — https://livermetabolism.com/paper/Raman2024_FROG.pdf , https://www.ebi.ac.uk/biomodels/curation/fbc
5. FROG tooling itself is fragile: BioModels lists runFROG as "not actively updated" and Fluxer "unmaintained"; fbc-curation 0.3.2 (2025-10-21) wraps cobrapy+cameo (cameo is dead) — https://pypi.org/project/fbc-curation/
6. Standards drift: FBC v3 encodable only via cbmpy/libSBML; COBRApy has no FBC v3 support, so v3 user-defined constraints are lost on read; PR #1440 addresses only the v3 metadata parts (annotations, key–value pairs, float charges), not user-defined constraints — https://github.com/opencobra/cobrapy/pull/1440
7. Namespace fragmentation: bigg.ucsd.edu frozen since 2019 while bigg.bio re-grounds IDs; MetaNetX concedes contradictory cross-refs — https://bigg.bio/about , https://academic.oup.com/nar/article/54/D1/D617/8343507
8. The SBML Test Suite (v3.5.0) has no advertised fbc semantic cases — https://github.com/sbmlteam/sbml-test-suite
9. Classic reference (not re-fetched this session): Ebrahim et al. 2015, "Do genome-scale models need exact solvers or clearer standards?", https://doi.org/10.15252/msb.20156157

Post-review additions (4 Sept 2026, from the independent review of the roadmap draft): (a) counts of "active" reactions (item 1) are not a reproducible quantity under non-unique optima, which is why FROG scores objective, FVA and deletions — item 1 is a definitional issue, not an open bug; (b) FBCModelTests.jl (LCSB, Kratochvíl) is an actively maintained FROG+MEMOTE implementation with CI, so FROG tooling is not uniformly unmaintained — https://github.com/LCSB-BioCore/FBCModelTests.jl ; (c) the MARS component of Persephone is Python (github.com/ThieleLab/mars-pipeline; authors include Hulshof, Nap, Martinelli, Widder), so Persephone is not entirely MATLAB-only; (d) MetaboAnnotator (Preciat et al., Bioinformatics 2022; Thiele/Fleming) already performs metabolite formula, charge and structure assignment in the Toolbox and should be carried forward rather than rebuilt; (e) in Machado's benchmark the MILP failures were GLPK/COIN on the single-organism Recon3D problem — SCIP and HiGHS did not fail.

## (D) Agentic / LLM interfaces

- **ChatGEM** (Chowdhury et al., PNNL, bioRxiv 2026-07-21): multi-agent ADEPT framework + RAG over COBRApy; RAG raised mean score 2.63→4.20 on three task tiers; applied to ecGEMs of *P. putida* — https://www.biorxiv.org/content/10.64898/2026.07.20.739662v2
- **Yeoh et al.** (NUS, bioRxiv 2026-06-08): benchmark of GPT-4/Gemini/Claude/DeepSeek-R1 on domain knowledge, flux prediction, pathway construction, optimisation; rubric-based LLM-judge ensemble, 9 metrics ×1–5; failure modes: context limits, wrong identifier assumptions, strain-dependent errors; "blind search" fails to localise injected stoichiometric faults unless prompted to call COBRApy mass-balance tools — https://www.biorxiv.org/content/10.64898/2026.06.03.730004v1
- **MechAInistic** (Helikar lab, bioRxiv 2026-05-08): Reviewer-supervised Architect agent over COBRApy; auditable evidence chains; immune-cell RA/MS targets — https://www.biorxiv.org/content/10.64898/2026.05.11.723319v1.full.pdf
- **Human2** (Li/Nielsen, PNAS 2026-04-08): LLM+expert reconstruction pipeline producing tissue/WBM/enzyme-constrained human models — https://www.pnas.org/doi/10.1073/pnas.2516511123
- **MCP servers**: gem-mcp (GSoC 2026/NRNB; 7 stateful tools, server-side model sessions, unit+protocol tests, 0 stars) https://github.com/praneeshlabs/gem-mcp ; FluxForge (18 tools, Recon3D+GLPK) https://glama.ai/mcp/servers/ali-kishk/FluxForge_MCP ; metabo-idmapper https://glama.ai/mcp/servers/Wooyoung-kim91/metabo-idmapper ; gem-flux-mcp (Faria, Argonne; see brief 06)
- **Skills**: K-Dense "cobrapy" Claude skill (2 installs; no evaluation) https://mcp.directory/skills/cobrapy ; Biomni has no dedicated GEM tools beyond a generic FBA call https://github.com/snap-stanford/Biomni
- Evaluation is uniformly weak: no shared task suite, LLM-as-judge dominant, no ground-truth flux datasets.

## (E) Opportunities for frontier-model intelligence

1. **Automated MATLAB→Python porting with differential tests.** Port PSCM/Persephone/DEMETER/thermoKernel to a COBRApy extension, using the Toolbox as an oracle. Success: bit-for-bit (within 1e-6) agreement on Harvey/Harvetta FBA and ≥3 Persephone tutorials.
2. **Maintainer-multiplier for COBRApy.** Triage 77 issues/20 PRs, rebase #1440/#1450/#1472, land FBC v3 read/write. Success: FBC v3 round-trip in a 0.33 release; median PR age <90 days.
3. **A cross-tool FROG-style conformance suite** (COBRApy/Toolbox/COBREXA/cobrar × GLPK/HiGHS/Gurobi/CPLEX) run in CI, with agents diagnosing divergences (tolerance, loop removal, default bounds). Success: published divergence matrix for ≥100 BiGG/BioModels models; ≥5 root-caused fixes upstream.
4. **Identifier-reconciliation agent** over MetaNetX 4.5 + bigg.bio + ModelSEED using structure (InChI) not strings. Success: ≥95% metabolite mapping accuracy on a curated hold-out; contradictions flagged with provenance.
5. **Solver-configuration advisor** for WBM/community scale: auto-scaling, tolerance selection, fallback from PDLP-GPU to simplex crossover. Success: reproduce optimizeWBModel objectives on open solvers within 1e-6 at ≤2× Gurobi time.
6. **Grounded GEM agent benchmark**: executable tasks (knockouts, media, gap-fill) with numeric ground truth from FROG reports, replacing LLM-judge scoring. Success: adopted by ≥2 agent papers; reports hallucinated-identifier rate.
7. **Media-to-constraints translator** from MediaDive/BacDive recipes to exchange bounds with uncertainty. Success: growth/no-growth prediction on BacDive phenotypes with reported AUC.

## (F) Most important sources

1. COBRApy PyPI release history (2026) — https://pypi.org/project/cobra/
2. COBRApy open PRs incl. FBC3 metadata #1440 (2025–26) — https://github.com/opencobra/cobrapy/pulls
3. COBRA Toolbox releases v3.5–v3.7 — https://github.com/opencobra/cobratoolbox/releases
4. COBRA Toolbox development plan — https://opencobra.github.io/cobratoolbox/stable/plan.html
5. Machado, solver benchmark, mSystems (2023/24) — https://journals.asm.org/doi/10.1128/msystems.00833-23
6. Raman et al., FROG analysis (2024) — https://livermetabolism.com/paper/Raman2024_FROG.pdf
7. BioModels FBC curation standard — https://www.ebi.ac.uk/biomodels/curation/fbc
8. MetaNetX 2026, NAR — https://academic.oup.com/nar/article/54/D1/D617/8343507
9. bigg.bio About (BiGG relaunch) — https://bigg.bio/about
10. HiGHS releases (GPU cuPDLP-C, PDLP caveats) — https://github.com/ERGO-Code/HiGHS/releases
11. Gurobi academic WLS restrictions (2026) — https://support.gurobi.com/hc/en-us/articles/34672988479633-What-are-the-restrictions-on-using-an-academic-WLS-license
12. ChatGEM (2026) — https://www.biorxiv.org/content/10.64898/2026.07.20.739662v2
13. Yeoh et al., LLM evaluation for GEMs (2026) — https://www.biorxiv.org/content/10.64898/2026.06.03.730004v1
14. Luo et al., Human2 via LLMs, PNAS (2026) — https://www.pnas.org/doi/10.1073/pnas.2516511123
15. COBREXA 2, Bioinformatics (2025) — https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaf056/8005852
