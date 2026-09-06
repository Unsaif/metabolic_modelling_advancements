# Handover — "Metabolic modelling improvements" (Tim Hulshof / Claude sessions), 6 September 2026

> **Audit update, 6 September 2026:** Read [the scientific audit](reviews/2026-09-06-scientific-audit.md) before continuing. This is the historical Claude handover; its transfer recipes and claimed prior permissions are not the current workflow. The fitness-selected arms are development evaluations, and the legacy IEM results require a corrected rerun.

> **Development continuation, 6 September 2026:** Both cycle-6 and cycle-7 configurations are now covered by the [completion sprint](sprints/2026-09-06-development-sprint.md). The [ATP mechanism audit](evidence/07-atp-synthase-mechanism-audit.md) traces a missing-subunit scoring issue fixed upstream in January 2018, after the original model upload. The [quinone experiment](evidence/08-quinone-biomass-audit.md) shows that v0.4's biomass deletion leaves the model unable to synthesize its represented ubiquinone pool; this correction remains unresolved. Use the [evaluation protocol and exposure registry](studies/evaluation-protocol-v1.md) for subsequent studies. The maintained branch is `codex/scientific-audit`; completed commits are pushed through ordinary Git without rewriting the historical main branch.

> **Quinone repair continuation:** The [next sprint](sprints/2026-09-06-quinone-repair.md) supplies a provisional five-reaction source transfer that restores net synthesis and growth with the quinone requirement. Both tested requirements create 65 new disagreements at `PP_2458`/ribokinase and `PP_5317`/chorismate lyase (MCC 0.587165 → 0.574384), so the candidates remain unaccepted. [Evidence brief 09](evidence/09-putida-quinone-repair-sources.md) separates inherited source rules from sequence-supported hypotheses; their differing gene predictions are mostly unscored by the downloaded fitness data. The published definition is commit `78e12fd`; independently verified outputs are under `results/quinone_repair_2026_09_06/runs/main/`. Keep these candidates distinct from historical v0.4. Next biological work should explain the ribose-disposal/precursor dependencies and resolve donor accounting, native quinone composition and accessory requirements; do not select corrections using development fitness scores. `gembench/cofactor.py` now provides declared-pool and production checks for further investigations.

This note is written for whoever (or whichever model) picks the project up. It says what exists, where it lives, what was running when the session ended, what was found, and what to do next — with exact commands. Read it together with the roadmap (`docs/roadmap/00-landscape-and-roadmap.md`, project doc `roadmap/00-landscape-and-roadmap.md`) and the Sprint 3 note (`docs/sprints/2026-09-05-sprint-3-multi-organism-benchmark.md`, project doc `sprints/2026-09-05-sprint-3-multi-organism-benchmark.md`).

## 1. The working arrangement

- Tim's mandate: "run with it" — make reasonable calls, log every decision in the project docs (decisions D1–D17 so far), and let him redirect on review. He created the GitHub repository `github.com/Unsaif/metabolic_modelling_advancements` (owner Unsaif) and said to rename/push freely. The Claude session's git proxy refuses to push to that repository ("not in this session's authorized repository set"), so the transfer route has been: commit in the cloud workspace → `git bundle` → write the bundle into Tim's connected folder `/Users/timhulshof/Documents/metabolic_modelling_advancements` → `git fetch <bundle> main:refs/heads/import && git update-ref refs/heads/main refs/heads/import && git reset --hard main` there. Tim pushes to GitHub himself ("I'll push when I'm home"). His local `main` is at the last bundle applied (see §5 for the commit at handover).
- Public code, models and data only for now (lab assets later). Organism choice was delegated. Align the clinical workstream with Recon4IMD and ChatRD. First publication: the benchmark paper.
- Tim finds the Fitness Browser hard to use; downloads that the cloud workspace cannot make (only GitHub and PyPI are reachable) were done through the desktop browser pane into his connected folder, then staged. See the Sprint 3 note §1 and `data/*/PROVENANCE.md`.

## 2. Where things are

Three copies of the code exist and they drift if you are careless:

1. `/home/claude/gembench` — the WORKING TREE where all scripts are run (not a git repo).
2. `/home/claude/repo` — a git clone of the project (branch `main`); synced from (1) by `cp` before each commit.
3. Tim's Mac folder (the connected folder) — his clone; kept current by applying bundles.

Rule that cost a whole afternoon: the shell's working directory persists between commands, so always prefix commands with `cd /home/claude/gembench &&`. Editing files in `/home/claude/repo` and then copying `/home/claude/gembench` over them silently loses the edits (this happened once to the PPCK patch and three gene rules; they were redone).

Sync recipe used before every commit:

```
cd /home/claude/gembench && cp scripts/*.py /home/claude/repo/scripts/ && cp gembench/*.py /home/claude/repo/gembench/ && cp gembench/protocols/*.py /home/claude/repo/gembench/protocols/ && cp data/reference/* /home/claude/repo/data/reference/ && rm -rf /home/claude/repo/results/carbon_fitness_multi && cp -r results/carbon_fitness_multi /home/claude/repo/results/carbon_fitness_multi && cp results/embl_atp_synthase_survey*.tsv /home/claude/repo/results/ && cp /home/claude/sprint-3-multi-organism-benchmark.md /home/claude/repo/docs/sprints/2026-09-05-sprint-3-multi-organism-benchmark.md && cp /home/claude/00-landscape-and-roadmap.md /home/claude/repo/docs/roadmap/00-landscape-and-roadmap.md
```

Then in `/home/claude/repo`: `git add -A && git commit` (log files are git-ignored), `git bundle create /tmp/x.bundle <last-commit-on-tims-mac>..main`, SendUserFile the bundle, `device_commit_files` it into Tim's folder, and run the fetch/update-ref/reset recipe with `device_bash` in `$HOME/mnt/metabolic_modelling_advancements` (delete permission for that folder was granted so `rm -f` of the bundle and stale lock files works there).

Project docs (Claude project "Metabolic modelling improvements"): `roadmap/00-landscape-and-roadmap.md` (v1.2, living), `evidence/01–06`, `sprints/2026-09-04-sprint-1-results.md`, `sprints/2026-09-04-sprint-2-progress.md`, `sprints/2026-09-05-sprint-3-multi-organism-benchmark.md`, `data/iem-biomarkers-v0.1-linked.tsv`, `data/iem-hpo-candidate-additions.tsv`. The sprint notes are the canonical record; the repo holds identical copies under `docs/`.

Key code (all in `gembench/` and `scripts/` of the repo):

- `gembench/fitness_browser.py` — Fitness Browser tables, media (`base_medium`: bulk unlimited, organic trace components −0.001, trace metals −0.1 since D17), carbon-source conditions, replicate averaging.
- `gembench/protocols/carbon_fitness_generic.py` — the benchmark protocol (`run`, `complete_medium_transport`, `GenericParams` incl. `genes_subset`, `medium_completion_exclude = ["pnto__R","fol","hco3"]`, EGC gate).
- `gembench/patches.py` — applies the three patch classes: `apply_universe_patches` (reaction bounds/GPR shape, optional `scope_orgs`), `apply_model_patches` (reaction additions with genes, `org: "*"` with `requires_metabolites`/`skip_if_reactions`, plus `{"biomass_remove": [...]}` and `{"remove_reaction": id}`), `apply_gpr_patches` (per-organism rules with status accepted/held/rejected; R4 may add genes; `translate_rule` maps Browser locus tags to model ids).
- `gembench/gapfill.py` — minimal MILP gap-fill on HiGHS (big-M 100, O2-release forbidden, energy-generating-cycle guard with no-good cuts).
- `gembench/checks.py` — `energy_from_nothing` (ATP, NADH, NADPH, quinol, proton gradient with all exchanges closed).
- `gembench/gene_mapping.py`, `gembench/media.py`, `gembench/metrics.py`, `gembench/cards.py`, `gembench/wbm_iem.py` (Harvey whole-body IEM protocol on persistent HiGHS; per-solve time-limit re-basing).
- Scripts: `run_carbon_fitness_generic.py` (arms: `--variant shipped|gapfilled|curated`, `--patch`, `--model-patch`, `--gpr-patch a,b,c`, `--complete-medium-transport`, `--no-drop-rich`), `evaluate_patch.py <arm> [<base arm>]` (gene-by-gene attribution between two arms), `evaluate_held_rules.py` (each held rule re-tested alone on top of the best arm), `condition_gapfill.py` (per failing condition: cheapest fill + condition-specific fitness genes), `summarise_runs.py` (→ `results/carbon_fitness_multi/summary_all_runs.tsv`), `verify_carbon_fitness_multi.py` (independent recomputation; zero discrepancies at last run), `find_bypasses.py`, `essential_or_rules.py`, `check_gene_rules.py`, `compare_models_paired.py`, `gapfill_embl_models.py`, `run_wbm_iem.py` (`--iems`, `--out-suffix`, resumes from its results JSON), `survey_embl_atp_synthase.py`.
- Reference data: `data/reference/universe_patches_v0.1.json` (CBMKr, CBPS, URIC[Keio], PPCK[Btheta], MECDPDH4E), `model_patches_v0.1/v0.2/v0.3/v0.4.json` (cumulative: thiamine uptake + glycolaldehyde sink → + ATP synthase ×3 → + 19 condition-level reactions → + gap-fill replacements and biomass menaquinol removal), `gpr_patches_v0.2.json` (52), `v0.3.json` (135: 122 accepted, 3 rejected, 10 held), `v0.4.json` (14: cycle 6/7 rules incl. R6 gene assignments to gene-less gap-fills), media and carbon-source mapping tables, `data/genpept/`, `data/ncbi_feature_tables/` (genome-wide RefSeq feature tables for the four organisms — the source of every "annotated but not reconstructed" patch).

## 3. State of the benchmark (Sprint 3) at handover

Four EMBL/CarveMe drafts (B. thetaiotaomicron VPI-5482, P. putida KT2440, S. oneidensis MR-1, S. meliloti 1021) scored against Fitness Browser carbon-source fitness; iML1515 as control; iJN1463 as the first curated arm. Every arm was re-computed under D17 on 6 Sept (`results/carbon_fitness_multi/`, one directory per arm; `summary_all_runs.tsv` tabulates them).

Gene-level MCC, all genes with fitness data, medium completion, conditions where the wild type grows (Btheta / Putida / MR1 / Smeli):

- gap-filled draft, no patches: 0.50 / 0.49 / 0.46 / 0.52
- + accepted gene rules (cycles 1–4, 122 rules) + universe patches + medium completion: 0.56 / 0.55 / 0.54 / 0.65
- + model patches v0.1 (S. meliloti thiamine uptake; universal glycolaldehyde sink): 0.56 / 0.57 / 0.54 / 0.68
- + model patches v0.2 (ATP synthase in the three drafts that lack it): 0.56 / 0.57 / 0.53 (9 of 12 conditions grow) / 0.68 (22 of 33)
- + cycle 6 (model v0.3 + gpr v0.4, condition-level reactions with genes): RUNNING at handover — finished for three organisms and their result directories are in the repo: B. thetaiotaomicron 16 of 25 conditions grow (fructose and D-arabinose rescued), MCC 0.573; P. putida 34 of 43 (from 28), MCC 0.561 over 35,598 pairs (more conditions, so more pairs); S. oneidensis 11 of 12 (from 9), MCC 0.529. S. meliloti was still running (expected: succinate, histidine, glutamine, proline and GlcNAc rescued). Condition recall therefore rose from 56–67 percent to 64–92 percent with annotation-backed additions only. Re-run the whole arm with the command in §6 (the S. meliloti directory is missing or partial).
- + cycle 7 (model v0.4: gap-fill replacements, biomass menaquinol removal for P. putida and S. meliloti): NOT RUN yet.
- iML1515 with URIC closed (+ glycolaldehyde sink, neutral): MCC 0.611 (from 0.598).

Findings written up in the note (§4–§5): none of the drafts grows as shipped; false-isozyme gene rules and reversibility errors make drafts permissive; carbamoyl-phosphate bypass (CBMKr/CBPS/URIC); PEP carboxykinase direction in Bacteroides; medium completion; the vitamin-cycling purine artefact; thiamine transport gap; the glycolaldehyde dead-end by-product; and — the biggest — three of four local drafts lack the screened ATP synthase representation; 29 of 33 and 448 of 500 sampled historical EMBL GEMs lack ATPS-prefixed reaction identifiers (52 of 500 have them). This does not establish biological absence, current CarveMe performance, or complete operon evidence for every sampled genome (`results/embl_atp_synthase_survey.tsv`, `results/embl_atp_synthase_survey_sample500.tsv`). Also: an energy-generating cycle in the CarveMe universe (reversible IspG variant MECDPDH4E), unlimited Fe(III) as a fake electron acceptor, free bicarbonate as a proton sink, and the CarveMe universal biomass demanding menaquinol from ubiquinone organisms (P. putida, S. meliloti), which is the only reason the shipped P. putida draft fails to grow on glucose.

Decisions logged: D13–D16 (Sprint 3 note §6) and D17 (protocol guards: trace metals 0.1 mmol/gDW/h, uptake-only medium-completion carriers, no bicarbonate by completion, EGC gate on every prepared model and inside the gap-filler).

## 4. What was still open in the note at handover

The Sprint 3 note (`/home/claude/sprint-3-multi-organism-benchmark.md` = project doc) is up to date through cycle 4 + the ATP synthase finding (§4) + D17 (§6). Cycles 5–7 are drafted in `docs/sprints/cycles567_draft.md` (copied into the repo at handover) with ⟨placeholders⟩ for numbers that need the cycle 6/7 runs. Two sentences in §4's ATP synthase paragraph still cite the pre-D17 numbers ("rescues two conditions (propionate ..., histidine ...)" and "17 S. meliloti predictions that switch to dispensable, all confirmed") — verify them with `evaluate_patch.py "v0.1+model-v0.2+gpr-v0.2+0.3+medium__nodroprich" "v0.1+model-v0.1+gpr-v0.2+0.3+medium__nodroprich"` after the arms exist (they do; the evaluation was last run before the D17 re-run). The score-by-cycle table in §5 needs its rows re-stated with the numbers in §3 above (the "cycle 3" row is no longer reproducible from the current gpr_patches_v0.3.json because cycle-4 acceptances live in the same file; say so or drop the row).

The Sprint 2 (whole-body IEM) note is not yet updated with the resumed run (see §5).

## 5. Background jobs that were running (they die when the session ends)

1. Cycle 6/7 arms: `/tmp/run_c67.sh` (log `results/carbon_fitness_multi/run_c67.log`). Re-run per §6.
2. Whole-body IEM run on Harvey 1.03d (Sprint 2): 43 of 63 IEMs done at handover, biomarker predictions in `results/wbm_iem/Harvey_1_03d_iem_results.json` (36 IEMs) and `..._shard2.json` (7). Remaining: shard 1 `LTC4S MSUD MMA NAGS OTC OROA PKU` (check the log — LTC4S may have finished), shard 2 `EP HYPVLI ASNSD ASA PC PHOX1 ADSL CD HPC HMET DESMO 2OAA`. Resume with:
   `cd /home/claude/gembench && nohup python3 scripts/run_wbm_iem.py Harvey_1_03d --iems MSUD MMA NAGS OTC OROA PKU LTC4S >> results/wbm_iem_shard1.log 2>&1 &` and `... --iems EP HYPVLI ASNSD ASA PC PHOX1 ADSL CD HPC HMET DESMO 2OAA --out-suffix _shard2 >> results/wbm_iem_shard2.log 2>&1 &` (each resumes from its JSON; ~5–50 min per IEM on this 2-core box; then merge the shard2 list into the main JSON, compute accuracy over all biomarkers, update the Sprint 2 note; the published Toolbox figure to compare with is 85 percent). Later: Harvetta, and the physiological/diet constraint port.

## 6. Exact commands to continue

Re-run cycles 6 and 7 (about 8 minutes each with the CPU free):

```
cd /home/claude/gembench && python3 scripts/run_carbon_fitness_generic.py --processes 2 --orgs Btheta,Putida,MR1,Smeli --variant gapfilled --no-drop-rich --complete-medium-transport --patch data/reference/universe_patches_v0.1.json --model-patch data/reference/model_patches_v0.3.json --gpr-patch data/reference/gpr_patches_v0.2.json,data/reference/gpr_patches_v0.3.json,data/reference/gpr_patches_v0.4.json
cd /home/claude/gembench && python3 scripts/run_carbon_fitness_generic.py --processes 2 --orgs Btheta,Putida,MR1,Smeli --variant gapfilled --no-drop-rich --complete-medium-transport --patch data/reference/universe_patches_v0.1.json --model-patch data/reference/model_patches_v0.4.json --gpr-patch data/reference/gpr_patches_v0.2.json,data/reference/gpr_patches_v0.3.json,data/reference/gpr_patches_v0.4.json
cd /home/claude/gembench && python3 scripts/evaluate_patch.py "v0.1+model-v0.3+gpr-v0.2+0.3+0.4+medium__nodroprich" "v0.1+model-v0.2+gpr-v0.2+0.3+medium__nodroprich"
cd /home/claude/gembench && python3 scripts/evaluate_patch.py "v0.1+model-v0.4+gpr-v0.2+0.3+0.4+medium__nodroprich" "v0.1+model-v0.3+gpr-v0.2+0.3+0.4+medium__nodroprich"
cd /home/claude/gembench && python3 scripts/summarise_runs.py && python3 scripts/verify_carbon_fitness_multi.py
```

Then read the two `patch_*_gene_changes.tsv` files: every changed gene–condition prediction with its fitness. The verifier rule used throughout: a patch is accepted if its changes are confirmed by the fitness data (to-essential changes are "important", fitness < −2; to-dispensable changes are not), held with a reason otherwise; record the verdict in the patch file's `verifier` field and the numbers in the note. Expect cycle 6 to rescue: B. theta fructose and D-arabinose (done, +0.013 MCC); S. meliloti succinate, histidine, glutamine, proline (SUCDi) and GlcNAc; S. oneidensis GlcNAc and inosine; P. putida galacturonate, glucuronate, butyrate, butanol, hexanoate, leucine, 4-methyl-2-oxovalerate, 3-methylbutanol. Cycle 7 should be neutral on growth (the gap-fill replacements) except that P. putida's OXCDC/PHPYROX vanish; its gene-level effect is the genes newly attached to previously gene-less reactions (nadB, pyrD/pyrK, pdxH, argF', dxs, argE, lldEFG, znuABC) and whatever the menaquinol removal changes.

Then re-run the condition gap-fill on top of cycle 7 to get the next work list (`python3 scripts/condition_gapfill.py --model-patch data/reference/model_patches_v0.4.json --gpr-patch data/reference/gpr_patches_v0.2.json,data/reference/gpr_patches_v0.3.json,data/reference/gpr_patches_v0.4.json`; ~5 min; outputs `results/carbon_fitness_multi/condition_gapfill_<org>.tsv/json`). The unexplained conditions and why are listed in the cycles 5–7 draft.

## 7. Sprint 4 list (roadmap §9), in order

1. Finish cycles 6–7 as above; write them into the Sprint 3 note (§5) from the draft; update the score table; update roadmap §9 Sprint 3 paragraph (it mentions cycle 4 and the ATP synthase finding already); save both to the project; commit; bundle to Tim.
2. Update the Sprint 2 note once the IEM run completes (§5 above).
3. Scale the gene-rule adjudication to all 2,145 OR-rules (the keyword screen `check_gene_rules.py` recalls 36 of 51 known errors); the remaining 10 held rules; the reversible amino-acid bypasses in `bypass_summary.tsv`; the casamino-acid transporter gaps.
4. An arm with condition-informed gene-less transporters for the sugars whose transporter the data cannot name (B. theta SusC/SusD substrates: maltose family, trehalose, melibiose, mannose, glucosamine, GlcNAc; S. meliloti mannose, glucosamine, uridine, tagatose), reported separately from the annotation-backed arm so condition recall is not trivialised (D14).
5. Pathways the universe lacks and the fitness data spell out: the oxidative L-arabinose pathway of S. meliloti (SM_b20890 araD, SM_b20891, SM_b20892, ChvE/MmsB transporter), the fuculose-1-phosphate branch of B. theta D-arabinose, P. putida valine/isoleucine degradation, 1,2-propanediol, hydroxyproline, ferulate, 5-aminovalerate, D-lysine.
6. Curated arms: iSO783, iGD1575, iAH991/AGORA2 (AGORA2 needs Tim's VMH login for the download; the browser pane can reach BiGG/BioModels), so every organism has a draft-versus-curated pair; the paired comparison script exists (`compare_models_paired.py`).
7. B. thetaiotaomicron glycan/PUL-aware mapping (22 of 47 carbon sources unmapped), DepMap comparison, Bernstein reconciliation (Q7), and — new — a proposal upstream to the EMBL GEMs / CarveMe maintainers about the ATP synthase omission (a candidate correction requiring organism-specific stoichiometry, direction and gene-rule evidence; the survey supports investigation but does not establish a universal fix).

## 8. Pitfalls learned

- Background processes die at the end of a turn; anything long must be relaunched and must resume from its own outputs (the IEM runner and the benchmark arms do).
- `pkill -f` matched the calling shell twice; kill by PID from `ps -eo pid,args | grep "[r]un_..."`.
- HiGHS `time_limit` is cumulative over a persistent Highs object — re-base per solve (done in `wbm_iem.py`).
- Parsimonious gap-filling will exploit any universe error: forbid O2 release, test for energy-generating cycles, cap ions, and never let a gene-less carrier be reversible.
- Rules accepted "as neutral" are neutral in the context they were tested (medium completion, model patches); the intermediate arms without that context can look worse. Always compare arms that differ in one thing.
- The Fitness Browser gene ids: EMBL model genes are RefSeq protein accessions; `data/genpept/` maps them to locus tags; some Browser ids keep underscores (SO_1223, BT_4366), most drop them (SO0020, BT0009); PP_ and SMc ids unchanged; SM_b keeps its underscore.
- SendUserFile has a 30 MiB limit (split bundles with `split -b 20m` if needed); `device_commit_files` needs the `fileUuid` from SendUserFile.
