# Hand-curated reference models outside BiGG: plan (Paper 1, post hoc)

Written 6 October 2026, before any of the models below was scored. Committed before the scoring runs, so the
commit time shows the order.

## Why

Paper 1 compares the transfer study's draft arms (B0, U′, M) with a hand-curated model in one development organism
(*P. putida*, iJN1463 from BiGG). That comparison found the drafts about level with the curated model on shared genes
and behind on coverage. One organism cannot carry that finding. This extension scores the published curated models
of the other development organisms with the same protocol.

## Which models (rule fixed beforehand)

The curated model named for each organism in `data/fitness_browser/PROVENANCE.md`, committed on 5 September 2026
(6cabc1e) before any curated model was scored:

| Organism | Model | Source used | sha256 |
|---|---|---|---|
| *S. oneidensis* MR-1 (MR1) | iSO783 (Pinchuk et al. 2010) | BioModels MODEL1507180036, `MODEL1507180036_url.xml` | 2d4d60b9522bfe9d329ca614949e3af016ecf3a635d6df667ae8f9a1075d56ee |
| *S. meliloti* 1021 (Smeli) | iGD1575 (diCenzo et al. 2016) | Nature Communications Supplementary Data 6 (`41467_2016_BFncomms12219_MOESM594_ESM.zip`, SBML `ncomms12219-s7.xml`) | zip 1f01e9c7b0d5f03d690b53e1b3e2bb028656cd0e11daf94a6fe78381660fcb1b; xml 49a0b6e1d7a675f3fe6bc245772557d5fc4d2d58c579455832a2a9926d021894 |
| *B. thetaiotaomicron* (Btheta) | iAH991 (Heinken et al. 2013) | No public SBML was found. Either the lab's original file, or a rebuild from the published Supplementary Table S10a/S10b that passes the checks below | — |

Newer curated models exist and are **not** scored, so the choice cannot follow the results: iMR1_799 (Ong et al.
2014) and iLJ1162 (Luo et al. 2022) for *S. oneidensis*, iGD1348 (diCenzo et al. 2020) for *S. meliloti*. The paper
names them. iMR1_799 was compared with an earlier fitness dataset for MR-1 (Deutschbauer et al. 2011), and iGD1348
used Tn-seq data, so they are less independent of the screens than the models scored here.

## How each model is prepared (relabelling only)

`scripts/translate_curated_model.py` gives each model the BiGG identifiers the study's media and carbon-source tables
use.

- **Metabolites.** Identifiers that already are BiGG identifiers are kept. Older BiGG identifiers are mapped through
  the ModelSEED database's BiGG1 → compound → BiGG aliases. ModelSEED compounds are mapped through their BiGG
  aliases (ModelSEED database commit 194ac8a). Explicit decisions are listed with their evidence in
  `data/reference/namespace/curated_model_id_overrides.tsv`. Two metabolites that would get the same identifier keep
  their own.
- **Exchanges.** Exchange reactions are renamed `EX_<id>_e`. KBase-style two-step exchanges in iGD1575
  (`X_e0 <=> X_b`, `X_b <=>`) are collapsed to `X_e <=>` with the published bounds; this is the same constraint.
- **Gene rules.** One unparseable rule in iGD1575 (a missing closing parenthesis) is restored to the authors' rule,
  documented in `data/reference/namespace/curated_model_gpr_repairs.tsv`.
- **Check.** The script stops unless the optimum is unchanged, both with the published bounds and with every
  exchange open.
- **Nothing else changes.** Reactions, stoichiometry, bounds, gene rules and the objective are as published. There
  is no gap-filling and no curation.
- **Gene identifiers.** These are mapped to Fitness Browser sysNames by fixed rules: iSO783 by identity; iGD1575
  `smc…`/`sma…`/`smb…` → `SMc…`/`SMa…`/`SM_b…`; iAH991 `BT_0554` → `BT0554`.

## Scoring

`scripts/score_reference_model.py` uses the transfer study's fixed protocol:

- the organism's media and carbon-source tables;
- carbon uptake −10, growth threshold 1e-3, fitness threshold −2;
- medium completion and GLPK;
- the energy-cycle gate.

A model that fails the gate is reported and not scored. Conditions whose carbon source has no exchange reaction in a
model count as no growth, as for the drafts. Results are never overwritten.

## Analyses (`scripts/compare_reference_models.py`)

The analyses below are the same as for iJN1463, per organism:

- coverage (mappable conditions with growth);
- own MCC;
- four-way union MCC;
- curated − arm MCC differences on the union and on common genes, with 1,000 gene-bootstrap 95% intervals;
- error overlap of U′ and the curated model on common genes;
- the curated model's calls on genes only it contains;
- no-growth conditions that lack an exchange reaction.

How results will be described:

- A difference whose interval excludes zero is reported as a difference.
- Otherwise the result is reported as no clear difference, with the interval.
- Organisms are listed side by side and not pooled into one test.

Anything beyond this plan is labelled exploratory.

## iAH991 rebuild checks (only if the original file is not available)

A model rebuilt from the published tables is scored only if all of these hold:

1. It has the published counts: 1,488 reactions, 1,152 metabolites and 991 genes.
2. It reproduces the growth rates in Supplementary Table S3 under the conditions of Table S8.
3. It reproduces the sole carbon sources in Supplementary Table S1.
4. It reproduces the published gene-essentiality predictions where the supplement lists them.

The rebuild is done without access to fitness data. Every manual intervention is listed. If the lab's original file
arrives, the original is scored, and the rebuild is compared with it.
