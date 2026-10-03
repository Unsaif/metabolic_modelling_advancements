# Benchmark paper: skeleton and evidence map (draft for Tim)

3 October 2026, Claude (Opus 5.5). A starting structure for the first publication (roadmap decision D12), built only from results already in this repository. Every number links to its source. Framing, authorship and venue are Tim's decisions.

## Working titles
- Do corrections to draft metabolic models transfer to new organisms? A pre-registered, replicated test on twelve bacteria
- Fitness-blind curation of automatically reconstructed metabolic models, tested prospectively against genome-wide mutant fitness

## One-paragraph abstract (numbers as of 3 October 2026)
Automatically reconstructed genome-scale metabolic models are the only models available for most bacteria. Their gene-level predictions are rarely tested on organisms that played no part in their curation. We scored historical CarveMe drafts (EMBL GEMs) against RB-TnSeq gene fitness from the Fitness Browser. On four development organisms we derived five annotation-only correction rules and a written procedure for blind AI curation of gene–reaction rules. We froze both, with file hashes and external timestamps, and then scored them once on six new organisms. We then replicated the test unchanged on six more. The automatic rules raised MCC by +0.036 in development, +0.037 on the first panel and +0.036 on the second. Nine of ten evaluable new organisms improved and none worsened. Blind AI curation added +0.028 pooled over both panels (eight of ten up). Nearly all of that gain came from removing non-catalytic gene alternatives and joining enzyme subunits; assigning genes to gene-less reactions gave no net gain. Two of twelve organisms could not be evaluated because of rules fixed before testing. The evaluation code, manifests and decisions are public.

## Claims and their evidence

| # | Claim | Evidence (numbers) | Source |
|---|---|---|---|
| C1 | Historical CarveMe drafts agree only moderately with mutant fitness on carbon sources; none of the development drafts grows on its experimental medium as shipped | Baseline union-gene MCC on the ten evaluable new organisms 0.33–0.52, except *R. palustris* (0.14) and *M. tuberculosis* (about 0); development drafts need a minimal gap-fill to grow | `results/transfer_v1/evaluation_panel_{A,B}/paired.json`; Sprint 3 note |
| C2 | Five annotation-only rules improve agreement on organisms never used to develop them, by the same amount as in development | +0.036 / +0.037 / +0.036 (development / panel A / panel B); pooled +0.037 [0.025, 0.048], 9 up, 0 down, 1 unchanged; sign test p = 0.004 | `docs/studies/transfer-v1-replication-results.md` |
| C3 | Blind AI curation of gene rules from annotation alone adds a further, smaller gain; its useful part is narrow | Pooled +0.028 [0.009, 0.046], 8 up, 2 down; R1/R2 +0.026 [0.014, 0.038]; R6 +0.003 [−0.008, 0.013]; R6 errors concentrate on ion transporters | same; `transfer-v1-results.md` |
| C4 | A self-custodied but auditable prospective protocol is feasible | Two freezes per panel (git, manifests, Claude-project server timestamps); first outcome access after both; independent recomputation reproduced every number exactly | `results/study_freezes/`; verification sections of both results documents |
| C5 | The pipeline's own fixed rules set the limits | 2 of 12 organisms not evaluable (gap-fill energy-cycle cuts; reference rule); 2 more misrepresented (organic carbon in a "no carbon" medium; a phototroph run in the dark) | failure records in `results/transfer_v1/evaluation_panel_{A,B}/` |

## Figures and tables
- **Figure 1. Study design and time order** (to make): development organisms → frozen method → panel A inputs → freeze → outcomes → score; same for panel B.
- **Figure 2. Per-organism paired differences** (exists): `docs/studies/transfer-v1-replication-results.png`.
- **Figure 3. Attribution of the gain** (to make): reaction patches, ATP synthase rule, normalisation, menaquinone rule, R1/R2 curation, R6 curation; per panel.
- **Figure 4. Precision of curation decisions** (to make): new "important gene" calls confirmed by the data, by decision kind and organism.
- **Table 1.** Organisms, panels, media, conditions mapped and grown, gap-fill reactions, status.
- **Table 2.** Primary results per organism and pooled (the replication results' tables).
- **Supplement.** Adjudication procedure and prompt (verbatim); all decisions with their evidence; freeze manifests; media and carbon-source mapping tables.

## Methods outline
1. Data: Fitness Browser RB-TnSeq (Price et al. 2018 and later releases), carbon-source experiments; FEBA medium recipes.
2. Models: EMBL GEMs (Machado et al. 2018) pinned at commit 260d0f1; CarveMe bacterial universe.
3. Mapping: gene maps by protein sequence; media and carbon sources by written rules.
4. Base model: minimal gap-fill on one reference condition chosen by rule, with an energy-generating-cycle gate.
5. Corrections: five automatic transforms; blind adjudication by fresh AI subagents limited to an annotation packet.
6. Scoring: knockout growth versus fitness ≤ −2; MCC on matched genes and conditions; union-gene paired differences; gene and organism bootstraps; sign test; pre-declared reading.
7. Pre-registration and custody: manifests, external timestamps, quarantine of outcome files, independent recomputation.

## Limits to state plainly
- Self-custodied evaluation; separation by recorded order only.
- One assay platform and one historical draft collection; current CarveMe output is not tested.
- Ten evaluable new organisms; small effects (≈0.03 MCC).
- The AI curator's training may include published phenotypes for some organisms; development showed one procedure example leaking into one development organism.
- RB-TnSeq fitness in pooled competition is not single-mutant growth.

## Before submission
1. Decide whether to add a robustness analysis with current CarveMe drafts of the same twelve organisms. These organisms are now exposed, so this would be robustness, not a new independent test.
2. Reconcile the Bernstein et al. 2023 reference numbers (roadmap Q7) if that benchmark is cited.
3. Check data and model licences for redistribution (Fitness Browser, EMBL GEMs, BiGG/CarveMe universe).
4. Agree an AI-contribution statement and authorship with the lab.
5. Correct the hard-coded provenance strings in run cards (cosmetic; see the errata in both results documents).
