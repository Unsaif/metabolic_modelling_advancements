<!-- Saved verbatim by the main session from the verifying agent's final message (6 October 2026). The main session's actions on each discrepancy are in RESPONSE.md. -->

# Independent verification of IEM study v0.4 and the Paper 2 draft

6 October 2026. I checked every number with my own scripts, which are in `results/wbm_iem/independent_verification_v0.4/`. They import nothing from `gembench` or `scripts/`. The only project code I ran was `scripts/freeze_study.py verify`. I ran no LP solves and changed no other file or git state; `git status` shows only the new folder. Some web pages could not be read (rate limits, captcha, blocked hosts); they are listed in section (c).

**Bottom line**
- **Data and provenance hold.** Almost every count, value, hash, time and table entry reproduces exactly.
- **One major factual error.** The published Thiele et al. 2020 accuracies are swapped between Harvey and Harvetta. The error is in the frozen plan, in the results note and three times in Paper 2.
- **One false statement.** Paper 2 says every run finished on 6 October; two finished on 5 October.
- **Wording that goes beyond the data.** Several sentences overstate: the port reproducing "every stored bound", caps as "physiological caps", "all below 0.5%", and a verification described as already done. One reported transition table is incomplete.

## (a) Claims checked

| # | Claim | Verdict | My numbers |
|---|---|---|---|
| 1 | Runs A, B, C: 57 IEM records each, all `complete` | Confirmed | 57/57/57, all `complete`; summaries agree (673 solves each) |
| 2 | One run fingerprint per file, equal to the sha256 of the provenance (v0.3 recipe) | Confirmed | A `51d7bda7…`, B `92385e32…`, C `5b6c8d07…`; one provenance block per file; sha256(json.dumps(provenance, sort_keys=True)) matches for all three (and for v0.3 `c85ebde9…`) |
| 3 | B and C record senses = min; A is max-only with v0.3's provenance fields | Confirmed | B and C `senses: ["min"]`. A has no `senses` key and the same 14 provenance keys as v0.3. A's biomarker records carry 4 extra null `*_min` fields (informational). |
| 4 | B's `lp_bounds_sha256` equals v0.3's `693af5c9…` | Confirmed | Equal, and equal to the sha256 of MATLAB's Harvey global bounds. A and C are both `8b2d70b5…`. |
| 5 | Model and protocol hashes match the files | Confirmed | Harvetta `98810012…`, Harvey `10e6cb6d…`, protocol `6fed6856…`, constraint inputs `c63977ba…`. Both .mat files equal the blobs of COBRA.models commit 75c070d (hash-object vs ls-tree). |
| 6 | Freeze verifies | Confirmed | `freeze_study.py verify`: valid, 22 files, `f3dde4e2…`. My own recomputation of the content fingerprint matches. The manifest plan equals `data/studies/wbm_iem_v0.4_plan.json`. |
| 7 | Freeze (13:39:37Z) precedes any Harvetta v0.4 IEM result or biomarker minimum | Confirmed | Plan commit 350c526 13:39:28Z → manifest 13:39:37Z (snapshot: HEAD 350c526, clean) → freeze commit 9f9ebab 13:39:41Z → Claude-project freeze record 13:39:50.6Z (fingerprint matches) → runs A and B start about 13:39:58Z (file time minus session wall) → M submitted 13:40:27Z. First appearances in git: Harvetta MATLAB IEM results 62b41c0 (5 Oct 18:40:47Z), Harvey minima f595951 (19:51:43Z), A and C c14e702 (6 Oct 02:43:47Z). |
| 8 | Plan and freeze commits are ancestors of the results commits; frozen files unchanged | Confirmed | Ancestry true for A, B, C and M. `git diff 9f9ebab HEAD` over the 22 frozen paths is empty. |
| 9 | Harvetta 218/251 (86.9%) | Confirmed | 218/251 = 86.85%. Only HPC `EX_25aics[u]` is unscored (absent). Stored predicted/expected/correct fields are all consistent. |
| 10 | Harvetta errors: 15 opposite, 18 no change (16 capped, 2 zero) | Confirmed | 15 / 18 (16 capped, 2 zero, 0 other) |
| 11 | IEMs fully correct: Harvetta 35; Harvey v0.2/v0.2b/v0.3 37/38/38 | Confirmed | 35; 37/38/38 |
| 12 | Paper 2 run table (v0.2 220/251, 87.6%, 13/18; v0.2b and v0.3 217/251, 86.5%, 13/21) | Confirmed | Identical |
| 13 | 12 shared opposite errors, 12 shared no-change errors; 10 of Harvetta's 16 capped errors also capped in Harvey | Confirmed | 12; 12 (10 capped plus CYP21D aldosterone and cortisol); 10 |
| 14 | The three Harvetta-only opposite errors and their values | Confirmed (wording: discrepancy 12) | **ASNSD `DM_asn_L[bc]`:** Harvetta healthy 1.34e-6 (MATLAB prints 0), disease 49.71 → Increased; Harvey 130.748 → 49.74, Decreased. **MMA `DM_crn[bc]`:** Harvetta 32.818 → 50.000; Harvey 50 = 50 (capped). **BTD `EX_3hpp[u]`:** Harvetta 1.852e-6 → 0 (Decreased); Harvey 0 = 0. |
| 15 | MMA blood 3-hydroxypropionate is right in Harvetta, on a 0.3% change | Confirmed | 16308.4 → 16355.5 (+0.288%). Harvey: −1.08%, an opposite error. |
| 16 | Paper 2 opposite-error table (directions, expected, "wrong in") | Confirmed | All 13 rows match both models' values |
| 17 | Zero no-change errors: CYP21D aldosterone and cortisol (both models), BTD 3-hydroxypropionate (Harvey) | Confirmed | Exactly these |
| 18 | Effect sizes: Harvey 217/212/210/207/207, Harvetta 218/215/214/213/212 | Confirmed | Identical, with or without zeroing values ≤1e-6. No call is gained at any τ. |
| 19 | Correct calls below 5%: 10 (Harvey) and 5 (Harvetta); the five named Harvetta calls | Confirmed | GACR `DM_gln_L[bc]` 2.48e-5; HPC `DM_C05770[csf]` 1.80e-4; HPC `DM_C05770[bc]` 2.46e-4; MMA `DM_3hpp[bc]` 0.288%; GACR `EX_glu_L[u]` 2.39% |
| 20 | Harvey: five correct calls below 0.1%; smallest is HPC coproporphyrin III at about 10⁻⁷ | Confirmed | 5; HPC `EX_C05770[u]` at 1.32e-7 |
| 21 | Figure 2 bars | Confirmed | Calls lost: Harvey 5/7/10/10, Harvetta 3/4/5/6 |
| 22 | MATLAB on Harvetta: 31 non-finite = 30 without an optimum + 1 absent | Confirmed | 30 with healthy 'NaN' and finite disease; HPC `EX_25aics[u]` is 'NA'/'NA' |
| 23 | Same call on all 221 MATLAB-finite biomarkers | Confirmed | 221/221. All 30 differences have a NaN MATLAB healthy value. |
| 24 | MATLAB 191 correct → 75.8% of 252 and 86.4% of 221 | Confirmed | 191/252 = 0.75794, equal to the stored `Accuracy`; 191/221 = 86.43%. Python is also 191 there. |
| 25 | HiGHS solved all 30 MATLAB failures; Python right on 27; the 3 wrong are LNS folate, PC glucose, BTD 3-hydroxypropionate | Confirmed | All 30 Optimal in both states. LNS and PC are also wrong on Harvey, where MATLAB has values ('0'→69.9995 and '0'→2720.3266) and agrees. |
| 26 | All 291 value pairs above 1e-3 agree within 0.008% | Confirmed | 291 pairs; maximum 0.0079% (BTD `EX_CE2026[u]` healthy, limited by MATLAB printing about 4 significant digits). The 181 smaller pairs differ by ≤9.4e-6. |
| 27 | Per-IEM failure table (PC 9, BTD 6, HPII 5, LNS 3, CIT1 2, XAN1 2, ADSL/ASNSD/DGK 1) | Confirmed | Identical |
| 28 | Harvey MATLAB numbers | Confirmed | **Failures:** 14 without an optimum + 1 absent; 13 call differences, all with NaN healthy. **14th case:** BTD `EX_3hpp[u]`, both implementations Unchanged. **Figures:** 238/251 same call; 204/252 = 80.95%; 204/237 = 86.08%. **Values:** 293 pairs above 1e-3; FED healthy 0.298%, next 0.0298%. **Table:** the 13 are ASNSD 2, BTD 4, EF 1, GMT 1, PHOX1 3, SUCLA 2. |
| 29 | Harvey certification (violations ≤2e-7; independent re-solve within 2.6e-5) | Confirmed from the v0.3 verifier's files | Worst row violation 1.0e-7, worst bound 2.0e-7; maximum absolute difference 2.6e-5 (wording: discrepancy 20) |
| 30 | Harvetta bounds identical to MATLAB after the global constraints | Confirmed | sha256 of MATLAB's global lb+ub = `8b2d70b5…` = runs A and C. My own transcription of runIEM_HH's global step, applied to MATLAB's setup bounds, reproduces MATLAB's global bounds bit for bit (83,521/83,521, lb and ub). Counts 334/34/28 match the provenance. |
| 31 | Harvetta bounds identical after the setup | Confirmed for all entries the global step does not overwrite | No setup change falls on the 368 lb / 62 ub entries the global step overwrites, so the hash match implies equality everywhere else. Those 430 entries rest on the study's own comparison file. |
| 32 | Paper 2 bound-change table, Harvetta (kidney 1,010 ×0.700; CSF 508 ×1.486; 3 BBB; 41 diet; 1,000/562; older GFR 128.64) | Confirmed | **Harvetta:** kidney 1,010 (961 lb + 49 ub) ×0.69963, implying an older GFR of 128.640 ml/min; CSF 508 ub ×1.48571 (= 0.52/0.35); 3 BBB lb; diet 31 lb at −0.1 plus 5 lb (−12) plus 5 ub (−8); totals 1,000 lb / 562 ub. **Harvey:** 1,011 ×0.69365 (GFR 129.749), 517, 3, 41; 1,001/571. |
| 33 | Harvey v0.2 opened "all 261" bile-duct exits; the Toolbox list has 28 | Confirmed | 261 `BileDuct_EX_` reactions in Harvey (262 in Harvetta); 28 set by the global step |
| 34 | Flux ranges: 251 biomarkers per model; Harvey 217/207/217; Harvetta 218/205/218 | Confirmed | Identical; only HPC `EX_25aics[u]` excluded; same with or without zeroing |
| 35 | All R1 losses are Conflicting: disease maximum rises, disease minimum falls to zero | Confirmed | 10/10 and 13/13; disease min = 0 and healthy min >1e-6 in every case |
| 36 | Capped biomarkers with zero minima in both states 17/18 and 13/16; differing minima 1/18 and 3/16 | Confirmed | Identical |
| 37 | R2 turns 1 and 3 errors into "Decreased" and fixes none | Confirmed | Harvey: BTD `EX_2mcit[u]`. Harvetta: the same, plus MMA `DM_c4dc[bc]` and PC `EX_acetone[u]`. All expected Increased; healthy min >0, disease min 0. |
| 38 | Note's R1 transitions "10 / 13 right → conflicting, 0 fixed" | Incomplete (discrepancy 7) | right→wrong 10/13; wrong→right 0/0; wrong→other wrong **1/3** |
| 39 | Readings: R1 harmful, R2 not useful | Confirmed | Under the frozen rules |
| 40 | Post hoc: 217/233 (93.1%) and 218/235 (92.8%) | Confirmed | 93.13% and 92.77% |
| 41 | Wall times 12.9 h, 6.1 h, 6.8 h | Confirmed | Session wall 46,569.8 s, 22,084.3 s and 24,550.3 s; the sums of per-IEM times agree |
| 42 | MATLAB: 3.05 h, 15:28Z–18:31Z; first attempt 13:40Z–14:34Z, ended by a Docker check | Confirmed | **Second attempt:** 15:28:09Z–18:31:09Z (3.050 h). **First attempt:** submitted 13:40:27Z, error at 14:34:26Z ("docker container inspect … timed out"), no automatic retry. **Code:** both attempts used sha256 `3398850d…` = the frozen job script. **Files:** the repository copies are byte-identical to the runner's outbox files. |
| 43 | One process each, no shards | Confirmed (logs only; logs are not in git) | One header per run log. The watcher started C at 19:49:02Z. An unused `logs/shard_c.sh` was written at 00:35Z; nothing shows it ran. |
| 44 | Parameter drift (v0.2b → v0.3) changed no call; urinary maxima scale by 0.69 | Confirmed | 0 call changes; 55 urinary values scale by 0.6936 |
| 45 | v0.2 → v0.2b: 5 changed calls, net −3, values as in the table | Values confirmed; "all below 0.5%" not (discrepancy 5) | **DPYR:** 120.28 → 120.85 (0.47%). **HCYS:** 0.30% each. **HYPRO1:** 28.22 = 28.22, then 0 → 28.22. |
| 46 | Cap examples: 2.592 mmol/day filtration limit; carnitine at 50 | Confirmed | 557 kidney exchange bounds at −2.592 (554 lb, 3 ub) in each model; `Diet_EX_crn[d]` lb is −50 |
| 47 | "Each such value sits at a cap set by the constraints"; uracil 84.11 as a filtration or excretion limit | Not confirmed (discrepancy 6) | Only 6/18 (Harvey) and 7/16 (Harvetta) capped values equal a bound of the same metabolite; 9/18 and 8/16 equal any bound. No bound in Harvey's LP has magnitude 84.1137. |
| 48 | runIEM_HH was run with only two edits | Confirmed | Diff against the Toolbox 67c790d copy: exactly the `edit` line and the Harvetta load line |
| 49 | IEM flux truncated to six decimals; NaN when there is no optimum | Confirmed | `checkIEM_WBM.m` uses `fix(x*1e6)/1e6`; it writes 'NaN' unless solver stat is 1 or 3; it zeroes \|f\| ≤ 1e-6 |
| 50 | (Implicit) minima runs use the same states as the maxima they are compared with | Confirmed | Same IEM reaction sets (57/57); IEM-flux maxima equal within 1.6e-11 |
| 51 | Ref 1 citation; 57 IEMs and 252 biomarkers | Confirmed | Mol Syst Biol 2020;16(5):e8982 |
| 52 | Ref 1: 85.3% Harvey, 84.9% Harvetta, 50.2% Recon3D | **Not confirmed** (discrepancy 1) | Paper: Harvetta 215/252 (85.3%), Harvey 214/252 (84.9%), Recon3D 103/**205** (50.2%) |
| 53 | Refs 2 and 3 | Confirmed | Nat Protoc 2019;14(3):639–702; BMC Syst Biol 2013;7:74 |
| 54 | Ref 4 as the HiGHS citation | Citation correct; insufficient for the method used (discrepancy 15) | Math Program Comput 2018;10(1):119–142 |
| 55 | Ref 5 citation; "recall fell to 0.10 on a clinical database" | Citation confirmed; claim could not be checked | Mol Syst Biol 2009;5:263, PMID 19401675 |
| 56 | Gurobi version 12 | Could not check | Run records say only "LP solver gurobi" |
| 57 | An independent re-implementation matched the port on both models before any MATLAB check | Could not check | Only stated in the v0.3 plan text; no outputs in the repository |

## (b) Discrepancies

**1. Published accuracies swapped (major).**
- **Paper 2, Introduction:** "The direction was right for 85.3% of biomarkers in Harvey and 84.9% in Harvetta, against 50.2% for the generic reconstruction Recon3D \[1\]."
- **Paper 2, Results:** "…right in Harvey (86.5%) and 218 in Harvetta (86.9%). The published figures are 85.3% and 84.9%, from other model and code versions."
- **Paper 2, Discussion:** "The published 85.3% and 84.9% come from other model and code versions \[1\]."
- **Results note:** "The published figure for Harvetta is 84.9% (Thiele et al. 2020, other model and code versions)."
- **Frozen plan:** "the published accuracies: Harvey 85.3% and Harvetta 84.9%" and "The published 84.9% is shown for context only."

What is wrong: Thiele et al. 2020 write "a total of 215 out of 252 (85.3%) … by Harvetta and 214 out of 252 (84.9%) by Harvey". Two copies of the article agree, and so do the repository's own `results/iem_ground_truth/provenance.json` and `docs/evidence/04-clinical-and-rare-disease.md`.

Correct value: Harvey 84.9% (214/252); Harvetta 85.3% (215/252).

**2. "Every run finished on 6 October" is false.**
- **Paper 2, Note for Tim:** "every run finished on 6 October."
- Run B finished on 5 October at 19:48Z and run M on 5 October at 18:31Z. Only A (6 October, 02:36Z) and C (02:38Z) finished on 6 October.

**3. Port reproduction overstated.**
- **Paper 2:** "With those two values the port reproduces every stored bound of both models except three blood–brain-barrier uptakes." The frozen plan has the same wording.
- The "legacy" check (`scripts/check_wbm_constraint_port.py`, `physiology_legacy`) covers only the physiological constraints.
- Re-applying the current diet still changes 41 bounds in each model (31 AGORA-essential uptakes, and 5 bile acids on both bounds). No parameter value explains these, and the very next table lists them.
- Correct: "…reproduces every stored physiological-constraint bound except three blood–brain-barrier uptakes; the diet differs in 41 more."

**4. Verification described as done before it happened.**
- **Paper 2, Methods:** "An independent agent with its own code recomputed every reported number from the result files."
- When the draft was exported (6 October, 06:41Z), the v0.4 check had not been done; the Note for Tim says the numbers "are being checked".
- The v0.3 verification also left one claim not confirmed as stated (certified optima, 12 of 13).
- This check finds the errors in this list.

**5. "Five calls resting on differences below 0.5%."**
- **Paper 2, Abstract:** "one protocol detail flipped five calls resting on differences below 0.5%".
- **Paper 2, Results:** "every one of them rested on a difference below 0.5%".
- Four of the five (DPYR ×2 at 0.47%, HCYS ×2 at 0.30%) did.
- The fifth, HYPRO1 urinary hydroxyproline, went from equal maxima (28.22 = 28.22) to 0 → 28.22, a 100% change caused by the healthy maximum dropping to zero.

**6. Caps presented as established physiological caps.**
- **Paper 2, Abstract:** "values held at a physiological cap".
- **Paper 2, Results:** "Each such value sits at a cap set by the constraints, for example: … a measured filtration or excretion limit (uracil, 84.11 mmol/day)".
- **Paper 2, Discussion heading:** "Caps cause most 'no change' errors".
- The data show equal positive maxima only. 6/18 (Harvey) and 7/16 (Harvetta) capped values equal a bound of the same metabolite (kidney filtration or diet carnitine); 9/18 and 8/16 equal any bound.
- No bound in Harvey's LP has magnitude 84.1137: `EX_ura[u]` ub is 0.0051 and `Kidney_EX_ura` lb is −0.54. The uracil example is therefore unsupported.
- The 2.592 and 50 examples are correct.

**7. R1 transitions incomplete (results note, Q2 table).**
- **Note:** "10 right → conflicting, 0 fixed" and "13 right → conflicting, 0 fixed".
- R1 also turns 1 (Harvey) and 3 (Harvetta) no-change errors into "Decreased" errors. These are the same biomarkers that R2 changes.
- The study's own `v0.4_report.json` records `wrong_to_other_wrong` 1 and 3, and the plan requires the transitions to be reported.

**8. Harvetta results existed before the freeze.**
- **Paper 2, Introduction:** "under a plan fixed before any Harvetta result existed".
- **Paper 2, Note for Tim:** "frozen … before any of its results existed".
- The MATLAB setup-bounds comparison (finding (a), job at 13:32Z) and the port check came before the freeze; the plan itself says Q1(a) "holds already".
- Correct: "before any Harvetta IEM result or biomarker minimum existed".

**9. Figure 1, point 5.**
- **Paper 2:** "Caps: 18 of the 21 no-change errors have the same maximum in both states."
- All 21 have the same maximum in both states; 3 of them are zero in both. It should say "the same positive maximum".

**10. "Always" is contradicted by the data.**
- **Paper 2:** "The model can always send the metabolite elsewhere".
- In 1/18 (Harvey) and 3/16 (Harvetta) capped biomarkers the healthy minimum is positive (0.0061, 17.18 and 0.0042 mmol/day). The results note correctly omits "always".

**11. "Differ only because" overstates.**
- **Paper 2:** "The headline accuracies differ only because the Toolbox scores solver failures as 'no change'."
- runIEM_HH also divides by 252, which includes the absent HPC biomarker.
- With Python's calls on the failures, MATLAB would score 217/252 = 86.1% and 218/252 = 86.5%, not 86.5% and 86.9%.

**12. ASNSD wording.**
- **Note:** "The healthy maximum is zero in Harvetta (MATLAB agrees)".
- **Paper 2:** "a healthy maximum of 131 mmol/day in Harvey but zero in Harvetta, in MATLAB as well as in Python".
- Python's Harvetta value is 1.34e-6, above the 1e-6 zeroing threshold; MATLAB prints 0. The call is unaffected.
- MATLAB has no Harvey value for this biomarker (it is one of Harvey's 14 failures). The 131 is Python's value alone.

**13. A causal explanation stated as fact.**
- **Paper 2, Abstract:** "…cost 10 and 13 correct calls, because the protocol's healthy state forces maximal flux through the affected pathway".
- **Paper 2, Results:** "The second point follows from how the protocol defines health."
- The results note's "Why" section says the same.
- The data show the pattern (healthy minimum >0, disease minimum 0 in every loss), not the mechanism. It should be marked as an interpretation.

**14. Software versions not recorded.**
- **Paper 2:** "with the COBRA Toolbox at commit 67c790d and Gurobi 12".
- **Results note, Runs table:** "(Gurobi 12)".
- The run records show only "0 of 19 checked Toolbox files differ from commit 67c790d" and "LP solver gurobi". No Gurobi version is recorded anywhere in the repository's run records.

**15. Ref 4 does not cover the method used.**
- **Paper 2:** "HiGHS 1.15.1 \[4\], interior point with crossover".
- Huangfu & Hall 2018 is the citation HiGHS asks for, but it describes the dual simplex.
- The runs used `ipm`, HiGHS's interior point solver (IPX: Schork & Gondzio, Math Program Comput 2020;12:603–635). Keep ref 4 and add the IPX reference.

**16. Recon3D denominator.**
- **Paper 2:** "against 50.2% for the generic reconstruction Recon3D".
- 50.2% is 103 of the 205 biomarkers Recon3D could represent, not of 252.

**17. Deviations list incomplete (results note).** It lists the plan's GFR slip but not the plan's swapped published figures (discrepancy 1).

**18. Wrong cross-reference (results note, A2).**
- **Note:** "BTD urinary 3-hydroxypropionate (the threshold case below)".
- In the note that case is described above, under "Opposite-direction errors".

**19. MATLAB run count (Paper 2, Limitations).**
- **Paper 2:** "MATLAB was run once per model".
- The Harvetta job was submitted twice. Correct: one completed run per model.

**20. Certification wording (Paper 2, Harvey).**
- **Paper 2:** "The solutions were certified against the LP the solver held".
- The certified solutions come from a re-run, because the original v0.3 solutions were not saved. That re-run differs from the reported BTD `EX_nh4[u]` value by 8.8e-5 (v0.3 verifier, claim 7).

**21. Omissions (minor).**
- Paper 2's code-availability table leaves out `scripts/v04_report.py` and `scripts/iem_error_anatomy.py`, which produced the v0.4 numbers.
- The results note does not report results by biofluid, although the plan declares them. They exist in `v0.4_range_calls.json` and match mine. R1 drops blood from 102 to 97 and urine from 107 to 103 in Harvey; blood from 105 to 98 and urine from 105 to 100 in Harvetta.

## (c) Could not check, and notes

1. **Shlomi et al. 2009, recall 0.10:** I could not reach the full text. The PMC pages (incl. pmc.ncbi.nlm.nih.gov/articles/PMC2683725) were behind a reCAPTCHA. Europe PMC full text, link.springer.com/article/10.1038/msb.2009.22 and scholar.archive.org returned HTTP 429. NCBI E-utilities was refused (robots.txt; also blocked by the shell's egress policy), and Crossref was blocked by the shell policy. The figure appears to come from the project's own evidence brief, not from the paper. It needs a manual check.
2. **Gurobi version:** not recorded in any run file (discrepancy 14).
3. **The pre-MATLAB independent re-implementation:** no outputs exist in the repository.
4. **Setup bounds on the 430 entries the global step overwrites:** checking them needs the port, which I was not allowed to run.
5. **MATLAB first attempt:** "the queue stopped the container" is not in the job record, which shows only the Docker timeout.
6. **Times:** git and file times are self-reported, so only their order is evidence. Run start times are inferred from file time minus session wall time. `logs/` is excluded from git.
7. **Thiele 2020 publisher page:** embopress.org returned HTTP 429. The figures come from the Europe PMC full text and the Leiden repository PDF, which quote the same sentence.
8. **Not tested:** that COBRApy drops coupling constraints (the model files do contain C, d, dsense and ctrs), and the literature statements about prior robustness work.
9. **Informational:** `v04_report.py` was committed at 22:44Z on 5 October, while A and C were still running and after B and M had finished. It loads the frozen analysis scripts, and its outputs match mine.

Files are in `/home/claude/mma/results/wbm_iem/independent_verification_v0.4/`:
- Scripts: `vcommon.py`, `verify_v04_scores.py`, `verify_v04_matlab.py`, `verify_v04_bounds.py`, `verify_v04_provenance.py`, `verify_v04_caps.py`, `verify_v04_state_pairing.py`, `collect_v04_numbers.py`.
- Main output: `recomputed_numbers_v04.json` (all key numbers).
- Per-check outputs: `scores_v04.json`, `per_biomarker_v04.tsv`, `matlab_v04.json`, `bounds_v04.json`, `provenance_v04.json`, `caps_v04.json`, `state_pairing_v04.json`, `freeze_verify_output.json`.

Sources:
- [Thiele 2020, Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7285886/fullTextXML)
- [Thiele 2020, Leiden repository copy](https://scholarlypublications.universiteitleiden.nl/access/item%3A3134879/download)
- [Shlomi 2009, Europe PMC record](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:10.1038/msb.2009.22&format=json&resultType=core)
- [Heirendt 2019, ORBilu](https://orbilu.uni.lu/handle/10993/39264)
- [Ebrahim 2013, eScholarship](https://escholarship.org/uc/item/2018k3q6)
- [Huangfu & Hall 2018](https://webhomes.maths.ed.ac.uk/hall/HuHa13/)
- [HiGHS website](https://highs.dev)
- [HiGHS on Wikipedia](https://en.wikipedia.org/wiki/HiGHS_optimization_solver)
