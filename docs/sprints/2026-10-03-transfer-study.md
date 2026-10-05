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

## Replication on panel B (same day, after the panel A results)

Tim confirmed that no external custodian is available, so the replication is self-custodied in the same way. Order: replication plan frozen (manifest 19dc96c, fingerprint `09312d9d…`, also in the project) → panel B inputs prepared under the frozen v1 rules and frozen (manifest d790b27, fingerprint `922b12b8…`, also in the project) → first download of panel B fitness tables → every arm run once → independent verification by a separate agent (exact agreement). Full write-up: [docs/studies/transfer-v1-replication-results.md](../studies/transfer-v1-replication-results.md).

- **H1 replicated.** Automatic rules on panel B: +0.036 MCC [0.024, 0.051], five of five organisms up. Development, panel A and panel B estimates agree: +0.036, +0.037, +0.036.
- **H2 supported on panel B**: +0.029 [0.011, 0.047], five of five up.
- **Pooled over both panels (ten evaluable organisms):** H1 +0.037 [0.025, 0.048], nine up and none down. H2 +0.028 [0.009, 0.046], eight up and two down. Both are supported by the pre-declared reading.
- Within curation, removals and joins (R1/R2) carry the gain (pooled +0.026); gene assignments to gene-less reactions (R6) do not (pooled +0.003).
- **H3 (new, secondary):** dropping R6 assignments to inorganic-ion transport reactions. It pointed the right way (+0.003; Caulo +0.005, PV4 +0.012) but is inconclusive, because only two organisms had such assignments.
- Miya (*Desulfovibrio*) was not evaluable: the frozen reference rule found no candidate. Its fitness table was never downloaded.

## Correction to the development sprint record

`gembench.patches.apply_gpr_patches` applies a gene-rule patch whose genes are absent from the draft only for class R4. Four of the six accepted R6 assignments in `data/reference/gpr_patches_v0.4.json` name genes the drafts do not contain (MR1 ACGAMK → SO3507; Btheta ASPO2y → BT3184 and OCBT_2 → BT3717; Smeli DXPS → SMc00972), so they were skipped with only a log message in the cycle-6/7 arms of the [development sprint](2026-09-06-development-sprint.md). Only DHORD6 and PDX5PO2 (Btheta) were applied. The development numbers stand as computed; their description as including the cycle-6/7 R6 assignments is wrong for those four. The transfer study's `gembench.transfer.apply_decisions` adds absent genes explicitly and is not affected.

## Whole-body IEM protocol v0.2 (complete)

The corrected protocol (per-IEM state isolation, shipped bounds only) finished on Harvey 1.03d on 3 October. The run was interrupted twice when the workspace was reclaimed; it resumed from its results file each time, and the last 22 IEMs ran as two disjoint shards that were merged afterwards (same run fingerprint `8f9b6aa0…`). Results: `results/wbm_iem/Harvey_1_03d_iem_results_v0.2.json` and `Harvey_1_03d_iem_summary_v0.2.json`.

- All 57 IEMs completed; 251 of 252 expected biomarkers were scored (one exchange is absent from the model).
- **220 of 251 biomarker directions are correct (87.6%)**, and 37 of 57 IEMs have every biomarker right.
- Of the 31 errors, 18 predict no change where a change is expected and 13 predict the opposite direction. The reversals include citrulline in three proximal urea-cycle disorders (CPS1, NAGS, OTC), where the model raises blood citrulline instead of lowering it.
- These numbers use the shipped model bounds without the physiological constraints and the average European diet of the published protocol. They are **not comparable with the published 85 percent**; the next step is to port those constraints and rerun.
- **Correction (4 October 2026):** the previous point is wrong. Porting the constraints showed that Harvey 1.03d is shipped with them already applied: its `SetupInfo` records the EU average diet, and the port reproduces the stored bounds once two older parameter values are used (GFR 129.75 rather than 90 ml/min; CSF export from 0.35 rather than 0.52 ml/min). v0.2 therefore ran with the constraints as released with the model. v0.2 also deviated from `runIEM_HH.m` in one step: it set ub = 100 on all 261 bile-duct exits instead of the Toolbox's 28. See [the v0.3 plan](../studies/wbm-iem-v0.3-plan.md).

## IEM v0.3 and the MATLAB reference run (4–5 October)

Results: [wbm-iem-v0.3-results.md](../studies/wbm-iem-v0.3-results.md).
- **Port confirmed against MATLAB.** The Python port of the Toolbox's diet and physiological setup equals the COBRA Toolbox's own setup, bit for bit, for all 81,094 reactions. The MATLAB run used the Toolbox on Tim's Mac with Gurobi.
- **Accuracy.** With the protocol as the Toolbox runs it today, 217 of 251 biomarker directions are correct (86.5%).
- **MATLAB comparison.** MATLAB makes the same call wherever it returned an optimum. It scores 204 of 252 (81.0%) by `runIEM_HH`'s formula, because 13 solves without an optimum count as "no change".
- **What moved the calls.** The bile-duct correction changes 5 calls (net −3); the newer constraint parameters change none.

## Next

1. Write up the benchmark result (draft models, the transfer study and its replication) as the core of the first paper.
2. A method v2 (trace counter-ions limited in the media rule; base-media organic carbon flagged; phototroph handling; more gap-fill cut rounds; a reference rule that skips media without named carbon sources; no R6 assignments to ion transporters) needs new data: every Fitness Browser organism with an EMBL GEMs draft has now been used.
3. Port the physiological and diet constraints to the IEM protocol and rerun v0.2 with them.

## How to continue

```
python3 scripts/run_transfer_study.py run --org <org> --arms B0,UNQ,M      # results/transfer_v1/<role>/<org>/<arm>/
python3 scripts/transfer_compare.py --role evaluation_panel_A --orgs Cola,Dino,Dyella79,MycoTube,PS,RPal_CGA009 --pairs B0:UNQ,UNQ:M --out paired.json
python3 scripts/transfer_aggregate.py paired.json --metric union
python3 scripts/freeze_study.py verify --root . --manifest results/study_freezes/transfer_v1_method.json
```
