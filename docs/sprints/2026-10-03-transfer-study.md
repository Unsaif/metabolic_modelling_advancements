# Sprint 5 (3 October 2026): the first independent test

Claude (Opus 5.5), continuing from Codex's audit branch (`codex/scientific-audit`, 2ae5273) on branch `claude/opus-continuation`.

## Why this sprint

The [scientific audit](../reviews/2026-09-06-scientific-audit.md) found that every improvement so far was measured on the four organisms whose fitness data had guided it, so none could be called validated. This sprint set up and ran the first test on organisms that played no part in development, with the method fixed in writing before their outcomes were downloaded.

## What was done, in order

1. **Panel declared before any download** (commit aa25b6d): eligibility criteria E1–E6 for Fitness Browser organisms with an exact-strain EMBL GEMs draft. Twelve were eligible; a seeded random split gave panel A (Cola, Dino, Dyella79, MycoTube, PS, RPal_CGA009) and a reserved panel B (Caulo, Cup4G11, Marino, Miya, PV4, SB2B) (commit 8d0d52b). Only gene and experiment metadata were downloaded.
2. **Gene maps by protein sequence** for panel A (commit 9a702ab), because the newer drafts name genes by WP_ accessions.
3. **Development arms on the four development organisms**: the gap-filled base (B0), the universal rules (reaction patches, conjunction-rule normalisation, glycolaldehyde diffusion, ATP synthase if annotated, menaquinone rule) in combinations, and blind AI adjudication of gene rules by fresh subagents working only from annotation packets (procedure v1). Development result: U′ (= UNQ) +0.036 MCC over B0 (4/4 up); M (U′ + blind adjudication) +0.028 over U′ (3/4 up).
4. **Rule-based media and reference-condition mapping** written and checked against the development tables (the FEBA recipe mapper reproduces all six hand-curated development media; the reference rule reproduces all four hand-chosen references).
5. **Method freeze** (commit 81883ca; manifest e43a3c2, 190 files, fingerprint `7ffb4213…`; also recorded in the Claude project): [docs/studies/transfer-method-v1.md](../studies/transfer-method-v1.md). It fixes hypotheses, arms, the primary metric (paired MCC difference on the union of both arms' genes), the organism-level aggregate and the reading rule.
6. **Panel A inputs prepared under the frozen rules** (commit 2364a2d): media mapped (R2A undefined), 16 new carbon-source names mapped, reference conditions chosen, base models built (PS failed: every gap-fill solution created an energy-generating cycle), blind adjudication of the five others (265 candidates: 25 R6, 2 R2, 76 R1, 162 abstain). **Inputs freeze** (manifest d15bbed, fingerprint `ddf371e2…`, also recorded in the project).
7. **First download of panel A fitness tables** (commit 3b6096a), then every arm run once.
8. **Independent verification** by a separate agent with its own code: all numbers reproduce exactly, both freezes verify, and the time order is consistent.

## Result

Full write-up: [docs/studies/transfer-v1-results.md](../studies/transfer-v1-results.md).

- **H1, automatic rules: supported.** +0.037 MCC on average over five evaluable organisms (organism bootstrap 95% interval +0.017 to +0.054); four up, none down, one unchanged (*M. tuberculosis*, where the model has no signal at all).
- **H2, blind AI curation: inconclusive.** +0.027 on average (−0.007 to +0.060): +0.04 to +0.07 in Cola, Dino and Dyella79, slightly negative in *M. tuberculosis* and −0.031 in *R. palustris*.
- Two organisms expose pipeline limits rather than rule failures: the *M. tuberculosis* "no carbon" medium contains organic carbon (asparagine, citrate, ethanol), so the model grows everywhere; *R. palustris* is a phototroph that the model can only run on trace nitrate.
- Removals and joins of gene alternatives (R1/R2) helped again (+0.020); gene assignments to gene-less reactions (R6) were mixed again, with losses from metal-ion transporter assignments in both phases.

## Correction to the development sprint record

`gembench.patches.apply_gpr_patches` applies a gene-rule patch whose genes are absent from the draft only for class R4. Four of the six accepted R6 assignments in `data/reference/gpr_patches_v0.4.json` name genes the drafts do not contain (MR1 ACGAMK → SO3507; Btheta ASPO2y → BT3184 and OCBT_2 → BT3717; Smeli DXPS → SMc00972), so they were skipped with only a log message in the cycle-6/7 arms of the [development sprint](2026-09-06-development-sprint.md). Only DHORD6 and PDX5PO2 (Btheta) were applied. The development numbers stand as computed; their description as including the cycle-6/7 R6 assignments is wrong for those four. The transfer study's `gembench.transfer.apply_decisions` adds absent genes explicitly and is not affected.

## Whole-body IEM protocol v0.2 (background)

The corrected protocol (per-IEM state isolation, shipped bounds only) is running on Harvey 1.03d: 24 of 57 IEMs finished, 92 of 105 expected biomarker directions correct so far. These numbers use shipped bounds without the physiological and diet constraints, so they are not comparable with the published 85 percent. The run resumes from `results/wbm_iem/Harvey_1_03d_iem_results_v0.2.json` with `python3 scripts/run_wbm_iem.py Harvey_1_03d`.

## Next

1. Method v2, developed on development organisms and panel A (now exposed), tested once on panel B: trace counter-ions limited in the media rule; "no carbon" media with organic carbon flagged; phototroph handling; more gap-fill cut rounds; no R6 assignments to metal-ion transporters.
2. Ask an external custodian (for example a lab member) to hold panel B's outcomes, so that v2 is not self-custodied.
3. Finish the IEM v0.2 run; then port the physiological and diet constraints.

## How to continue

```
python3 scripts/run_transfer_study.py run --org <org> --arms B0,UNQ,M      # results/transfer_v1/<role>/<org>/<arm>/
python3 scripts/transfer_compare.py --role evaluation_panel_A --orgs Cola,Dino,Dyella79,MycoTube,PS,RPal_CGA009 --pairs B0:UNQ,UNQ:M --out paired.json
python3 scripts/transfer_aggregate.py paired.json --metric union
python3 scripts/freeze_study.py verify --root . --manifest results/study_freezes/transfer_v1_method.json
```
