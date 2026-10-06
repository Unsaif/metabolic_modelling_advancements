# Metabolic modelling advancements

Working repository for the "Metabolic modelling improvements" project: bringing frontier-model
effort to constraint-based metabolic modelling, with verification first. Tim Hulshof (Thiele lab,
University of Galway), with AI-assisted research and independent review.

**October 2026:** first drafts of two papers. The [benchmark paper](docs/paper/paper1-draft.md) reports a pre-specified, replicated test of draft-model corrections on ten held-out bacteria. The [second paper](docs/paper/paper2-draft.md) (open whole-body IEM simulation and the robustness of its predictions) reports the [IEM v0.3](docs/studies/wbm-iem-v0.3-results.md) and [v0.4](docs/studies/wbm-iem-v0.4-results.md) studies: the Python port matches the COBRA Toolbox bound for bound on Harvey and Harvetta and makes the same call wherever MATLAB returns an optimum; 217 and 218 of 251 directions are correct; and comparing flux ranges does not fix values whose maximum is the same in health and disease.

Start with the **[6 September scientific audit](docs/reviews/2026-09-06-scientific-audit.md)** and its **[reproducible artifacts](results/audit_2026_09_06/)**. The benchmark improvements are retrospective development results; independent validation remains to be done. The IEM port required corrections and its legacy results must not be used as validation.

The follow-up development sprint completes both existing patch configurations across all four draft organisms. Read the **[new sprint results](docs/sprints/2026-09-06-development-sprint.md)** together with two mechanism investigations: the **[historical ATP-synthase scoring issue](docs/evidence/07-atp-synthase-mechanism-audit.md)** and the **[closed quinone pool hidden by a biomass deletion](docs/evidence/08-quinone-biomass-audit.md)**. A higher development score does not resolve these biological questions. The **[evaluation protocol](docs/studies/evaluation-protocol-v1.md)** records exposure and the procedure for a future independent test; its file-freeze utility is available now.

The **[quinone repair sprint](docs/sprints/2026-09-06-quinone-repair.md)** supplies a source-derived five-reaction candidate that restores net synthesis and permits growth with a quinone requirement. That requirement also creates 65 new disagreements with measured gene fitness, so the candidate is not accepted as a biological correction. The sprint includes explicit gene-evidence grades, a coverage audit, saved candidate models, independent numerical/metric verification and reusable [cofactor checks](docs/studies/cofactor-audit-method.md).

The **[PpnP recycling sprint](docs/sprints/2026-09-06-ppnp-repair.md)** tests a source-supported explanation for 34 of those disagreements. Two provisional nucleoside recycling reactions remove the ribokinase dependency while preserving quinone demand; all 34 changed predictions agree with the exposed development measurements. The other 31 quinone-related disagreements remain. Controls recover the predicted dependency when the new route is disabled. The enzyme assignment remains inferred from homology, and neither the new gene nor a key upstream gene has an exported fitness row. A separate [numerical audit](docs/studies/growth-rank-numerics.md) shows why small changes in raw ranking scores require caution even when all binary predictions reproduce.

The **[precursor-supply diagnostic](docs/sprints/2026-09-06-precursor-supply.md)** investigates those remaining 31 disagreements without changing the model or its score. The PP_5317 mutant reaches 95% of parental growth in all 31 conditions with hypothetical external 4-hydroxybenzoate supply, at the exact amount required by the model's quinone demand. Transport and downstream synthesis controls support that mechanism. The source audit retains the existing UbiC gene assignment; the experimental audit confirms that all 67 contributing gene-fitness measurements remain above the defect threshold. Pooled experimental fitness and isolated deletion growth measure different things. Actual precursor sharing remains unproven, and all 31 benchmark errors stay reported.

The **[medium-preparation sprint](docs/sprints/2026-09-06-medium-preparation.md)** fixes a demonstrated state-handling problem: switching media can leave withdrawn supplements available when compartment labels are mixed. The new opt-in API uses explicit compartment roles and returns an isolated model with a boundary report. All 156 fresh conditions preserve the inherited equations and gene rules; all five synthetic withdrawal controls become independent of prior use. The 113 draft conditions have independently checked numerical solutions. A declared 129-solve follow-up of the curated reference provides accepted cross-solver evidence for 28 conditions, retains eight numerical certificate gaps, and records seven unanimously reported infeasible conditions. Every result, including rejections, is independently checked. The suite passes 500 tests plus six subtests. This strengthens the software checks without changing a biological model or benchmark score.

**First independent test (3 October 2026).** The **[transfer study](docs/sprints/2026-10-03-transfer-study.md)** froze a fitness-blind correction method on the four development organisms, then scored it once on six new organisms whose fitness data were downloaded only after two recorded freezes. The automatic, annotation-only rules transferred: +0.037 MCC on average over the five evaluable organisms (95% interval +0.017 to +0.054; four up, none down), matching the development estimate. Blind AI curation of gene rules was inconclusive: +0.027 on average (−0.007 to +0.060), clearly positive in three organisms and negative in the two whose experiments the draft models cannot represent (a medium with hidden organic carbon; a phototroph). See the **[method](docs/studies/transfer-method-v1.md)** and **[results](docs/studies/transfer-v1-results.md)**; an independent agent reproduced every number. The evaluation is self-custodied, and the panel is small.

**Replication (3 October 2026).** The same frozen method was then tested once on the six reserved organisms. The automatic rules replicated: +0.036 MCC, five of five evaluable organisms up, against +0.037 on the first panel and +0.036 in development. Blind AI curation improved all five (+0.029). Pooled over both panels (ten evaluable organisms), the automatic rules give +0.037 (nine up, none down) and blind curation +0.028 (eight up, two down); both pass the pre-declared test. See the **[replication results](docs/studies/transfer-v1-replication-results.md)**.

Then read **`docs/roadmap/00-landscape-and-roadmap.md`** — the standing charter: where the field
stands in 2026, the diagnosis, the workstreams, the decision log, and the current sprint. The
`docs/evidence/` briefs hold the sourced findings behind it; `docs/sprints/` holds the results of
each sprint.

## Layout

```
docs/          roadmap, evidence briefs, sprint notes
gembench/      GEM-Bench: benchmark harness (media, datasets, metrics, benchmark cards, FROG fixtures,
               whole-body-model loading/solving, IEM biomarker protocol)
scripts/       runnable entry points (see gembench-README.md)
results/       benchmark cards and result tables produced so far
data/
  fitness_browser/   RB-TnSeq gene fitness for 9 organisms (Fitness Browser; see PROVENANCE.md)
  genpept/           NCBI GenPept-derived maps from EMBL GEM gene ids (RefSeq accessions) to Fitness Browser locus tags
  reference/         curated mapping tables: Fitness Browser media and carbon sources -> BiGG identifiers
  iem/               inborn-error-of-metabolism biomarker ground truth v0.1 (linked to Orphanet/OMIM/HPO/HMDB)
models/gapfilled/  the four EMBL draft models minimally gap-filled for their experimental media (added reactions in *_gapfill.json)
models/bigg/       curated BiGG source model iJN1463 (UC Regents research/non-profit terms; see models/bigg/NOTICE)
tools/         small helpers (browser-download collectors, GenPept parser)
external/      (git-ignored) clones of public repositories used as data/model sources — see gembench-README.md
```

## Reproducing

`gembench-README.md` lists the Python dependencies, the external repositories to clone into
`external/`, and the commands for each benchmark. Reference data that is too large or not
redistributable (HPO, Orphanet XML, Human-GEM, Recon3D, Harvey/Harvetta) is fetched by the
commands documented there rather than stored here.

## Historical status (5 September 2026; superseded where noted by the audit)

- E. coli carbon-source fitness benchmark (Bernstein et al. 2023 protocol) reproduced with bootstrap CIs.
- yeast-GEM 9.1.1 essential-gene test reproduced exactly.
- Multi-organism carbon-source fitness benchmark (Sprint 3): EMBL draft models of B. thetaiotaomicron, P. putida,
  S. oneidensis and S. meliloti scored against Fitness Browser fitness, as shipped and after a one-condition minimal
  gap-fill, with iML1515 as control — `results/carbon_fitness_multi/`, summary in `summary.tsv`.
- Harvey/Harvetta 1.03d solve with HiGHS (interior point, ~20 s) and GLPK in Python, coupling constraints included, solutions certified.
- IEM biomarker protocol and the Toolbox's diet and physiological constraints ported to Python. They are verified against a MATLAB run of `runIEM_HH.m`: identical bounds, and identical calls wherever MATLAB returned an optimum. On Harvey 1.03d, 217 of 251 biomarker directions are correct; see docs/studies/wbm-iem-v0.3-results.md.
- IEM ground truth v0.1: 279 tuples from the lab's benchmark, linked outward; HPO cross-check done.
- Fitness Browser data for 9 organisms downloaded (5 Sept 2026).

## Reproducing the audit

The audit was tested with Python 3.14 on macOS arm64; `requirements-audit.txt` records the exact installed packages. It is an audit environment snapshot, not a claim of compatibility with every platform or every optional workstream.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-audit.txt
.venv/bin/python -m pytest -q
.venv/bin/python scripts/verify_carbon_fitness_multi.py
.venv/bin/python scripts/audit_saved_results.py
.venv/bin/python scripts/survey_embl_atp_synthase.py --summarize results/embl_atp_synthase_survey_sample500.tsv
```

`audit_saved_results.py` recomputes existing arrays into a separate audit directory. `--correct-leakage` corrects historical card metadata only; numerical outputs are preserved. New simulations should use a fresh `--output-dir`, as in `results/audit_2026_09_06/carbon_fitness_multi/`. The corrected IEM runner uses protocol v0.2 and separate outputs; it requires external whole-body model inputs and a fresh run.
