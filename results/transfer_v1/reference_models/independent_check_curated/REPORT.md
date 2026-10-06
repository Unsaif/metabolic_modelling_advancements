# Independent check: curated reference models outside BiGG (iSO783, iGD1575)

6 October 2026. Written by an independent checking agent. It could not save files named REPORT, so the main session
saved its returned report here verbatim, apart from formatting. Scripts and outputs are in this directory; run each
with `python -I`.

## Summary

- **Reported numbers.** Every number in the new `MR1/iSO783` entry of `comparison.json` reproduces exactly from the
  saved matrices, with independent numpy code. Re-running the scoring protocol from scratch also reproduces the
  saved iSO783 matrix: growth values agree to 3e-9, no growth call differs, and the fitness values are identical.
- **Translation is relabelling only, for both models.** Reaction by reaction, stoichiometry, bounds, gene rules and
  objective are identical, except for the documented iGD1575 gene-rule repair. All 290 iGD1575 two-step exchanges
  collapse to equivalent single exchanges (21 are flipped, all correctly). Optimal growth is identical in 109 media,
  and single-gene deletions agree in all 120 tested.
- **iGD1575 gate failure.** Confirmed on the published SBML as read. There are at least 4 independent energy-generating
  cycles, and no single reaction removal closes them. The gate is the same code, applied at the same step, as for
  the draft arms.
- **Model choice and timing.** The choice of models follows `PROVENANCE.md` at 6cabc1e. The plan was committed 13 s
  before the first scoring run.
- **Problems.** None is critical. Four change numbers or wording (problems 1, 2, 3 and 9).

## Problems

### 1. iSO783's gene mapping drops 13 genes, 10 of them with fitness data (major)

- Identity mapping places 769 of 783 genes. 13 of the 14 unmapped genes are in the Fitness Browser under the
  RefSeq form of the same tag: the model's `SO0419` is listed as `SO_0419`. The 13 are SO0419, 0590, 1140, 1223,
  1900, 2304, 2784, 2937, 3068, 3079, 3130, 3442 and 3579.
- The GenPept map and the model's reactions confirm they are the same genes: nadC→NNDPR, serB→PSP_L, prpE→ACCOAL,
  ald→ALAD_L, aroC→CHORS, psd→PSD, dapB→DHDPRy.
- No MR-1 locus is listed in both forms, so a fallback rule is unambiguous.
- The drafts map 6 of these 10 genes through GenPept. This asymmetry takes them out of the common-gene comparison
  and, on the union, counts iSO783 as predicting "no effect" for genes it contains.
- Re-scored with the rule SOnnnn → SO_nnnn (`rescore_iso783.py`):
  - own MCC: 0.5473 (560 genes) → 0.5460 (570 genes);
  - four-way union for iSO783: 0.527 → 0.534;
  - union differences, curated − B0, U′ and M: +0.056 → +0.063 [−0.004, 0.133]; +0.043 → +0.049; +0.051 → +0.058
    [−0.002, 0.125];
  - share of the gap closed (U′ / M): 0.38 / 0.21 → 0.34 / 0.19;
  - common-gene differences change by 0.002–0.003.
- Over 50 bootstrap seeds, the corrected curated − B0 union interval excludes zero in 2 of 50, so it is borderline.
- **Suggested fix:** add the fallback rule and re-score as a disclosed deviation from the plan, which fixed
  "identity".

### 2. D-(-)-tagatose is mapped to `tagat__D`, which is not a BiGG identifier (major)

- `tagat__D` is not in the CarveMe universe, iJN1463, iML1515 or ModelSEED's BiGG aliases. The BiGG identifier is
  `tag__D`, which is in the universe and is cpd00589's BiGG alias.
- iGD1575 contains D-tagatose with a boundary reaction (`cpd00589_e0`). Its BiGG view has `EX_tag__D_e`, so
  "EX_tagat__D_e absent" in the translation report is wrong.
- Only Smeli has this condition, so no Smeli model can ever grow on it. The drafts lack tag__D anyway, so their
  scores do not change.
- The genuinely absent carbon sources are confirmed: D-lactate is cytosol-only in iGD1575, and asn__L, his__L,
  met__L, phe__L and pro__L are cytosol-only in iSO783 under any identifier.
- **Suggested fix:** leave the locked table and record the error. Either map tagatose to tag__D for all Smeli models
  (disclosed), or count Smeli as 32 mappable conditions (the study-wide 150/267 then becomes 150/266).

### 3. Under the protocol, iGD1575 grows with no carbon source (major, latent)

- The protocol closes only `model.exchanges`. The published SBML puts D-proline (`cpd00567_e0`) and its boundary
  species in compartment c0, so its ±1000 boundary is never closed. The same applies to 21 reversible `Demand_`
  sinks (tRNAs, tetradecenoate) and to the glycogen and biomass boundaries.
- With no carbon source the model grows at 57.5 h⁻¹ on 587 units of D-proline. With D-proline closed it still grows
  at 1.28 h⁻¹ through the tRNA-Asp sink. With both closed, growth is 0.
- All 33 Smeli conditions grow at 57.5–59.5 h⁻¹, tagatose included. The published file gives identical values.
- The other boundary reactions of the drafts and of iSO783 are export-only, so they are unaffected.
- **Suggested fix:** before any iGD1575 scoring, close every boundary reaction that is not part of the medium (or
  document a compartment fix), and mention this in the failure note.

### 4. The energy gate is weaker for the curated models than for the drafts (minor)

- Neither curated model has a periplasm, so the `h_p` currency is silently skipped: 4 currencies are tested, against
  5 for the drafts.
- An equivalent check, extracellular-to-cytosolic proton flow (h_e → h_c), gives 0 for iSO783 and all six draft arms
  and 1000 for iGD1575. No outcome changes.
- **Suggested fix:** state which currencies were tested, or add the h_e check for models without a periplasm.

### 5. The iGD1575 gate failure is not recorded anywhere (minor)

- `Smeli/iGD1575/` is empty, and there is no log and no Smeli entry in `comparison.json`. The plan says such a model
  "is reported".
- The published model gives ATP 1000, NADH 125, NADPH 125 and ubiquinol-8 125. The cycles:
  1. rxn00058 (reversed) with rxn10043, feeding rxn10042 (ATP synthase);
  2. rxn00379 with rxn09240 and rxn00237 (both reversed): ADP + Pi → ATP;
  3. rxn00324 with rxn08655 (reversed), rxn10122 and rxn10125 (reversed), driving ATP synthase;
  4. rxn00772 (reversed), rxn00778, rxn01545 (reversed) and rxn01649.

  No single-reaction cut exists.
- **Suggested fix:** write a failure record and cite it.

### 6. The translation script's built-in equivalence check is weak (minor)

- It compares optima only. For iGD1575 the "published bounds" check is the same as the "all exchanges open" check,
  so it could not detect a flipped or mis-bounded collapsed exchange.
- `check_translation.py` found none.
- **Suggested fix:** add the reaction-by-reaction comparison, and optima in the protocol media.

### 7. "Share of gap closed" is not in the plan (minor)

It is a ratio of differences whose intervals include zero. **Suggested fix:** label it exploratory or drop it.

### 8. The curated-only calls use 8 conditions, not iSO783's 9 (minor)

They use the conditions where both U′ and iSO783 grow, giving 43 calls, 16 confirmed. On all 9 conditions there
are 49 calls, 19 confirmed. **Suggested fix:** name the conditions in the text.

### 9. MR-1 does not repeat the *P. putida* pattern (wording)

- **Coverage.** U′, M and iSO783 each grow in 9 of 12 mapped conditions, but not the same ones: iSO783 grows on
  N-acetylglucosamine and not adenosine, and the drafts the reverse. "A second gap, larger for the drafts, is
  coverage" holds for *P. putida* only.
- **Union.** The union favours the curated model (0.527 vs 0.471–0.493), the opposite of "the union view favours the
  drafts".
- **Shared genes.** The curated model leads U′, M and B0 by +0.053, +0.066 and +0.072. Every interval includes zero,
  in all 50 seeds.
- **Error overlap.** 372 of U′'s 569 errors are shared: 65%, or 4.3× chance, not 8×. 348 of the 372 are missed
  important genes. iSO783 misses fewer important genes (362 vs 473) but makes more false "important" calls (169 vs
  96).
- Paper lines 119, 261 and 266 still say these models were not scored.

### 10. Two MR-1 dipeptide conditions are unmapped although BiGG has identifiers for them (minor)

- The carbon table leaves Gly-Glu and Gly-DL-Asp unmapped ("no BiGG identifier"), but the universe has `gly_glu__L`
  and `gly_asp__L`.
- iSO783 contains both, under old identifiers, and grows on them (0.60 and 0.45). The drafts have no exchange for
  either.
- "Level coverage" therefore holds only on the 12 mapped conditions. Including these two would give 11/14 against
  9/14.

### 11. Process (minor notes)

- **Model choice.** The models match PROVENANCE at 6cabc1e, unchanged since.
- **Timeline.**
  - downloads 10:32–10:42; translation 10:51;
  - plan mtime 10:52:48; plan commit 10:53:09;
  - iGD1575 run started 10:53:22; iSO783 card written 10:53:35; `comparison.json` written 10:56:03.
- **Absent exchanges.** The translation printed the absent-exchange lists before the plan, but no fitness outcomes.
  Worth one Methods sentence.
- **Protocol unchanged.** The protocol code and tables are unchanged since 5–6 September.
- **Checksums.** All plan sha256 values match, and the zip's SBML is identical to the extracted file.
- **Uncommitted.** The MR1 results and `comparison.json` were not yet committed.
- **iAH991.** An iAH991 rebuild appeared after 10:58 and was not checked.

## Numbers: reported vs reproduced

All reported values matched the recomputation exactly. The corrected column uses the SOnnnn → SO_nnnn gene map.

| Quantity | Reported = reproduced | Corrected map |
|---|---|---|
| iSO783 own MCC | 0.5473 (560 genes × 9 conditions) | 0.5460 (570 × 9) |
| Own MCC, B0 / U′ / M | 0.4699 / 0.4871 / 0.4785 | same |
| Coverage, iSO783 / B0 / U′ / M | 9 / 8 / 9 / 9 of 12 | same |
| Four-way union: size | 839 genes, 7 conditions, 5,873 cells | 843 genes |
| Four-way union: MCC | B0 0.4713, U′ 0.4925, M 0.4831, iSO783 0.5273 | iSO783 0.5337 |
| Share of gap closed, U′ / M | 0.379 / 0.211 | 0.337 / 0.188 |
| Common genes, curated − B0 | 434 × 7: 0.4754 vs 0.5477, +0.0723 [−0.0137, 0.1528] | +0.0699 |
| Common genes, curated − U′ | 434 × 8: 0.4951 vs 0.5482, +0.0531 [−0.0225, 0.1304] | +0.0513 |
| Common genes, curated − M | 437 × 8: 0.4830 vs 0.5490, +0.0659 [−0.0115, 0.1448] | +0.0641 |
| Union, curated − B0 / U′ / M | +0.0561 [−0.019, 0.132] / +0.0427 [−0.025, 0.115] / +0.0511 [−0.014, 0.122] | +0.0629 / +0.0494 / +0.0578 |
| Error overlap (3,472 cells) | U′ 569, iSO783 531, both 372, expected 87.0 | 580 / 542 / 383 / 89.3 |
| Curated-only genes | 123 genes, 8 conditions, 43 calls, 16 confirmed | 127 genes, 43, 16 |

## Checks that passed

- Every protocol identifier present in either model is the same compound: 83 identifier-compartment pairs in
  iSO783 and 93 in iGD1575.
- All 52 + 1,197 renames, and the q8/q8h2 overrides, are correct.
- iGD1575 maps 1,574 of 1,577 genes, 1,206 of them with fitness data. The three unmapped are the pseudo-genes
  Unknown and Spontaneous, and thiS, which has no Fitness Browser locus.
- In 20 spot-checks, the gene descriptions are consistent with the models' reactions.
- The repository's gate code gives the same values as an independent re-implementation.

## Could not check

- The downloads themselves: BioModels and Springer are blocked here, so the files were checked by checksum only.
- iSO783's fidelity to the original article: the BioModels file has no maintenance requirement (ATPM lower bound 0),
  and this was not compared with the article.
- The plan's literature statements about iMR1_799, iLJ1162 and iGD1348.
- iAH991.
