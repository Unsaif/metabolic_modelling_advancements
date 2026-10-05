# Independent verification of IEM v0.3, v0.2b and the MATLAB comparison

5 October 2026. A separate agent with its own scripts (`verify_*.py`, `resolve_nan_healthy.py` and `summarise_resolve.py` in this folder) recomputed every number from the raw files. Its scripts import nothing from `gembench` or the project's comparison scripts. The only project code it ran was `scripts/freeze_study.py verify`. The agent's environment blocked it from writing Markdown, so its final report was saved here verbatim by the main session.

| # | Claim | Verdict | Verifier's numbers |
|---|---|---|---|
| 1 | Scoring rule and counts | **Confirmed** | v0.2 220/251, v0.2b 217/251, v0.3 217/251 correct of scored. The one unscored biomarker in every run is HPC `EX_25aics[u]`, which is not in the model. v0.3 errors: 13 opposite direction, 21 "no change" (v0.2: 13 and 18). IEMs with every biomarker scored and correct: 37 / 38 / 38. The stored `predicted`, `expected` and `correct` fields match the recomputation in every record. |
| 2 | Prediction changes | **Confirmed** | v0.2 → v0.2b: exactly 5 changes (DPYR `EX_ura[u]`, `EX_56dura[u]`; HCYS `DM_Lhcystin[bc]`, `DM_hcys_L[bc]`; HYPRO1 `EX_4hpro_LT[u]`). All five are expected Increased, so 4 go from right to wrong and 1 from wrong to right. v0.2b → v0.3: no call changes, although values change. |
| 3 | Effect size (post hoc) | **Confirmed** | v0.3: 217, 212, 210, 207, 207 correct at τ = 0, 0.001, 0.01, 0.05, 0.10. v0.2: 220, 215, 209, 204, 204. v0.2b gives the same as v0.3. |
| 4 | MATLAB comparison | **Confirmed** | 252 MATLAB biomarkers, 15 non-finite: 14 with NaN healthy, plus HPC `EX_25aics[u]` as NA/NA. runIEM_HH accuracy 204/252 = 0.8095, matching the stored value. Same call for 238 of 251. All 13 differences have MATLAB healthy = NaN. On the 237 biomarkers with finite values and a direction, both get 204 correct. Of 488 finite value pairs, 293 exceed 1e-3; only FED `DM_chsterol[bc]` healthy differs by more than 0.1% (0.298%), and the next largest is 0.030%. |
| 5 | Bounds | **Confirmed** (raw hash) | sha256 of the MATLAB global lb + ub bytes is `693af5c9…`, equal to the v0.3 `lp_bounds_sha256`. The version with −0.0 normalised (`d31341ae…`) does not match (note 2). |
| 6 | Provenance and time order | **Confirmed** | One fingerprint per file, each equal to sha256 of that file's provenance. The provenance fields are as stated, and the model, protocol and constraint-input files hash to the recorded values. Freeze verify is valid (17 files, `66bd36b0…`), and the verifier's own recomputation gives the same fingerprint. The plan and freeze commits are ancestors of the first results commit. |
| 7 | Certified optima | **Not confirmed as stated (12/13)** | All six IEMs are present. All 13 healthy solves are Optimal; worst row violation is 1.0e-7 and worst bound violation 2.0e-7 (BTD `DM_acetone[bc]`). The objective check fails for **BTD `EX_nh4[u]`**: the certification re-run gives 43.900102586 against 43.900014763 in v0.3, a difference of 8.8e-5 (2.0e-6 relative). |

**Further checks for claim 5**
- The verifier's own version of runIEM_HH's global-constraint step, applied to the MATLAB setup bounds, reproduces the MATLAB global bounds bit for bit. The counts match the provenance (319 / 24 / 28).
- The model's stored bounds plus that step hash to `bd8de4c4…`, the v0.2b `lp_bounds_sha256`.
- MATLAB setup bounds versus stored bounds: 1,001 lower and 571 upper bounds change. Every change falls in a documented category:
  - kidney: 1,011 bounds, ×0.6936;
  - CSF: 517 bounds, ×1.4857;
  - 3 BBB uptakes;
  - 31 diet uptakes at 0.1;
  - 5 bile acids at −12/−8.
- None of those changes is overwritten by the global-constraint step, so the hash match also implies that the port's setup bounds equal MATLAB's.

**Further checks for claim 6**
- Commit order: plan commit `4012fd3` (4 Oct 10:02:23Z) and freeze `a0b6c51` (10:02:47Z) precede the first results commit `1dea776` (4 Oct 15:27Z) and the final one `c08d0bb` (5 Oct 11:49Z).
- The frozen files are unchanged from the freeze commit to HEAD.
- The 35 v0.3 and 36 v0.2b records in the checkpoint commit are byte-identical in the final files.

**Discrepancies and notes**
1. **Claim 7 miss.** The certification file comes from a fresh re-run, not from the original v0.3 solutions, which were not saved. That re-run reproduced all 23 calls of the six IEMs, but healthy values differ from v0.3 by up to 8.8e-5 (BTD `EX_nh4[u]`). The verifier's own independent solve gives 43.9000104, within 4.3e-6 of v0.3. The call does not change.
2. **Hash with −0.0.** The model file stores −0.0 in 36 lower and 715 upper bounds. Both MATLAB and the port keep them, so only the raw-byte hash matches.
3. **MATLAB parsing detail.** For HPC `EX_25aics[u]`, the Healthy row's column 3 reads "Healthy - Reported:" and the values are "NA". This is how checkIEM_WBM marks a biomarker that is missing from the model. No count changes.
4. **Timing.** Both Python runs were restarted once partway through. The restarted v0.3 session overlapped the MATLAB run (08:17–11:06Z on 5 Oct), and the final results commit came after the MATLAB commit `996deb2`.
   - No setting changed: every record has the same fingerprint, and those fingerprints already appear in the external freeze record (Claude project doc `studies/2026-10-04-wbm-iem-v0.3-plan-freeze.md`, 4 Oct 10:09Z).
   - The frozen files are unchanged and the checkpointed records are identical.
   - The original 02:06Z freeze cannot be checked in git because those commits were lost in the workspace reset.
   - Git timestamps are self-reported, so only the commit order is proven.
5. **Small values (informational).** 46 value pairs below 1e-3 differ by more than 1e-6. The largest is 1.8e-4 (PC `EX_acac[u]` healthy). No call changes.
6. **Tiny effects (informational).** 5 of v0.3's 217 correct calls rest on relative changes below 0.1%:
   - HPC `EX_C05770[u]` (1.3e-7) and two other HPC C05770 biomarkers;
   - GACR `DM_gln_L[bc]`;
   - MSUD `EX_ile_L[u]`.

**Extra: independent re-solve of the 14 healthy values that MATLAB returned as NaN**

The verifier built the LP from the model file and the MATLAB global bounds, following checkIEM_WBM, and solved it with HiGHS.
- The six IEM-flux maxima agree with v0.3 within 1.9e-7, and the IEM reaction sets are identical.
- All 14 solves are Optimal and feasible within 1.3e-7; their values agree with v0.3 within 2.6e-5.
- Dual-based upper bounds are valid but loose: up to about 0.6% above the value, and 78 for `EX_nh4[u]`.
- Together with MATLAB's own finite disease values, these bounds give Python's call for all 13 differing biomarkers.

**Outputs in this folder:** `scores.json`, `per_biomarker_calls.tsv`, `matlab_comparison.json`, `bounds.json`, `provenance.json`, `freeze_verify_output.json`, `certified_check.json`, `misc_checks.json`, `resolve_summary.json`.
