<!-- Saved verbatim by the main session from the checking agent's final message (6 October 2026). Actions taken are in RESPONSE.md. -->

# Independent check: hand-curated comparison in Paper 1

Every new number reproduces from the raw matrices, cards and Fitness Browser tables. The problems are in the wording:
- **Overstated similarity:** "as good as", "matches" and "brought … to the level" claim more than the data show.
- **Wrong denominators:** "screened" is used where the counts are BiGG-mapped conditions.
- **False exclusivity:** P. putida is not "the one study organism with a usable hand-curated model".
- **Unsupported ceiling:** nothing in the data shows a ceiling near 0.6.

Method: my own numpy/json/pandas code, with no `gembench` or `scripts/` imports, no LP solves and read-only git only. The gene bootstrap resamples gene rows directly, with my own seeds.

## (a) Claims

| # | Claim | Verdict | My numbers |
|---|---|---|---|
| 1 | iML1515 grows in all 32 *E. coli* conditions; MCC 0.59 | Confirmed | 32/32 conditions; MCC 0.5863 on 1,339 genes × 32 conditions (42,848 cells); 95% [0.513, 0.652]. No fitness value is exactly −2, so ≤ and < give the same result. |
| 2 | Shared genes: 764 genes, 28 conditions; B0 0.53, U′ 0.57, M 0.61, iJN1463 0.56 | Confirmed | 0.5294 / 0.5699 / 0.6068 / 0.5642. M's shared set is 765 genes (0.6066 on the 764). Fitness values are identical across the four models in every shared cell. |
| 2b | Curated − B0 +0.035 [−0.063, 0.140]; curated − M −0.042 [−0.139, 0.054] | Agrees in substance | +0.0347 [−0.069, 0.137] and −0.0424 [−0.148, 0.051]. Over 5 seeds the B0 bounds span [−0.070…−0.065, 0.130…0.148] and the M bounds [−0.148…−0.131, 0.048…0.062]. Curated − U′ (not reported in the paper): −0.006 [−0.120, 0.100]. |
| 3 | Union: iJN1463 0.45 against 0.49 to 0.58 | Confirmed | 1,423 genes × 28 conditions: iJN1463 0.447; B0 0.486, U′ 0.541, M 0.577. Not reported: on the union, curated − M = −0.130 [−0.215, −0.034], which excludes 0. |
| 4 | iJN1463 grows in 36 of 43 conditions, drafts in 28; the 28 are a subset of the 36 | Numbers confirmed; "screened" is wrong | 36/43 and 28/43; the same 28 conditions in B0, U′ and M, all inside the 36. P. putida was screened on 57 carbon-source conditions; 43 map to BiGG. All 7 of iJN1463's non-growth conditions, and 10 of the drafts' 15, lack an exchange reaction for the carbon source. |
| 5 | Twelve organisms: 33–75%, median 60%, 150/267; U′ 155/267 | Numbers confirmed; "screened" is wrong | Cards and matrices agree with each other and with Table 1: 150/267 (56.2%), median 59.8%, IQR 48–67%, U′ 155. The cards' `conditions_total` gives 321 screened carbon-source conditions, so 150/321 = 46.7%. |
| 6 | Same protocol; no gap-fill; identity gene mapping; iML1515 hash; Bernstein source | Confirmed, with 2 caveats | See the protocol note below the table. |
| 7 | 94% = (273 + 94)/(274 + 117) | Confirmed | Recomputed from the B0 and UNQ matrices of the 10 new organisms: 274/273, 117/94, 53 genes, 93.86%. |
| 8 | Supplementary Table S1 | Confirmed except one time | Lock times equal the manifests' `created_at`. Commits are 5 s, 0 s, 0 s and 0 s after them; the external records are 11 s, 7 s, 8 s and 24 s after the commits. Downloads ran 10:12:18.565–10:12:40.500 and 13:07:00.845–13:07:19.691; all hashes verify. All 4 manifests verify at their own commits. At HEAD only `filter_report.json` differs: its 9 old entries are unchanged and 5 panel B keys were added in 23d4291, the parent of the panel B lock commit. The 376,620 and 881,911 cell counts also reproduce. |
| 9 | Only P. putida and M. tuberculosis have hand-curated models in BiGG | True for BiGG; wider wording wrong | BiGG API (108 models): iJN1463 and iJN746 (P. putida), iEK1008 and iNJ661 (M. tuberculosis), none for the other 14 study organisms. |
| 10 | References 17 and 18 | Confirmed (Crossref) | Authors, title, journal, year, volume, issue, pages and DOI all match. The article names its model iJN1462. |
| 11a | "the corrected *P. putida* draft is as good as the curated one" | Not supported | 90% interval for curated − U′ is [−0.096, 0.086]. Equivalence would need a margin of ±0.10 (±0.13 for M), 2.6 times the rules' mean gain. B0 is just as indistinguishable (±0.125). |
| 11b | Heading "matches a hand-curated model where both grow" | Overstated | Same as 11a. It also omits "on shared genes": on its own genes iJN1463 scores 0.44 (all 36 conditions), and on the union it scores lowest. |
| 11c | "Much of the remaining disagreement … is therefore shared by curated models" | Inference invalid; conclusion partly true for P. putida | 710 of U′'s 1,190 errors (60%) are also iJN1463 errors, against 85 expected by chance. That includes 294 of U′'s 326 false "important" calls. But this cannot be derived from MCC levels, and it rests on one organism. |
| 11d | "a realistic target for a draft is nearer 0.6 than 1" | Overreach | Replicate screens agree with each other at MCC 0.94 for E. coli (0.88 across sets) and 0.90 for the P. putida shared genes (0.85 across sets). Noise does not cap the scale near 0.6. |
| 11e | "brought the draft's gene-level predictions to the level of a hand-curated model" | Not supported | B0 was already within noise of the curated model (+0.035 [−0.063, 0.140]). U′ − B0 on the same genes is +0.040 [0.007, 0.083]. P. putida is a development organism. U′'s carbamate kinase (CBMKr) patch fixes an error iJN1463 shares (September runs: 0.490 → 0.510). |
| 11f | "The larger remaining gap is coverage" | Interpretation, not shown | The two gaps are on different scales; the ranking rests on 11d and one organism. |
| 11g | "filling these gaps needs no gene-level outcome data" | True as worded; needs qualifiers | It still uses condition-level growth data. Most gaps lack any exchange reaction, and 14 of 57 P. putida carbon sources (54 of 321 across the study) have no BiGG metabolite. |

Protocol note for row 6:
- **Identical:** parameters (carbon uptake, both thresholds, `complete_medium_transport`, `medium_completion_exclude`, GLPK) and the `data/reference` tables, which last changed on 5 September.
- **No changes to the model:** `applied` is empty and `model_genes` equals `browser_genes`.
- **Hashes:** both model files match their cards. The GitHub copy at `dbernste/E_coli_GEM_validation/Models/iML1515.xml` is byte-identical (sha256 9c772d44…de28).
- **Code and data:** the scoring code is unchanged since the lock. Both reference fitness matrices rebuild exactly from the raw Fitness Browser tables (maximum difference 0).
- **Caveat 1:** the reference runs used Python 3.13.16 / numpy 2.5.3; the arms used 3.11.15 / 2.4.4. COBRApy and optlang versions are the same.
- **Caveat 2:** the BW25113 deletions were not applied to iML1515. The project's September runs from the same file did apply them.

## (b) Discrepancies, most serious first

1. **False exclusivity.** Quoted: abstract "In *P. putida*, the one study organism with a usable hand-curated model"; Limits "Only *P. putida* had a usable hand-curated model". Wrong: curated models exist for B. thetaiotaomicron (iAH991, Heinken, Sahoo, Fleming, Thiele 2013), S. oneidensis (iSO783 / iMR1_799) and S. meliloti (iGD1575). The repository's own `data/fitness_browser/PROVENANCE.md` lists all three. They are outside BiGG. Fix: "the only study organism with a hand-curated model in BiGG whose screens the pipeline can represent", and add to Limits that the three exist outside BiGG and were not scored.
2. **Similarity stated as equivalence.** Quoted: abstract "scored as well as the curated model"; heading "matches"; bullet "is as good as the curated one"; Discussion "did as well as a hand-curated model". Wrong: absence of a difference is presented as equivalence; the intervals are about ±0.1 wide, and the uncorrected draft is also inside them. Fix: "scored about the same as the curated model on shared genes and conditions (…; all differences within about ±0.1)"; report curated − U′ = −0.006 [−0.108, 0.103]; heading "On shared genes, the corrected *P. putida* draft scores about as well as a hand-curated model, which grows in more conditions".
3. **Abstract last paragraph.** Quoted: "in the one organism where it could be checked they brought the draft's gene-level predictions to the level of a hand-curated model". Wrong: other organisms could have been checked; it is a development organism, placed next to "carry over to new bacteria"; the uncorrected draft was already statistically indistinguishable. Fix: "In *P. putida*, a development organism, the corrected draft's gene-level agreement was similar to that of a hand-curated model (so, within the uncertainty, was the uncorrected draft's)."
4. **Denominators and scope.** Quoted: "28 of the 43 conditions in which the bacterium was screened" and similar. Wrong: 43 and 267 count conditions whose carbon source maps to BiGG; P. putida was screened on 57 conditions, and the twelve organisms on 321. Fix: name the mappable conditions and give the totals.
5. **Unsupported ceiling.** Quoted: "so on this benchmark a realistic target for a draft is nearer 0.6 than 1"; "First, in practice the scale does not run to 1." Wrong: one model's score is not a ceiling; replicate screens agree at MCC 0.85–0.94. Fix: "Even iML1515 reaches only 0.59, although replicate screens agree with each other at 0.94; the limit lies in the models and the binary growth call, not in noise in the screens."
6. **Shared disagreement.** Quoted: "Much of the remaining disagreement between models and screens is therefore shared by curated models." Wrong: a shared error set does not follow from similar MCC values. Fix: state the overlap (710 of 1,190, about eight times chance) or delete.
7. **Union dismissed as unfair.** Quoted: "so shared genes are the fair basis for models with different gene sets". Wrong: on the 384 genes only iJN1463 contains, 1,235 of its 1,288 "important" calls are unconfirmed; the union is the study's primary metric, and on it M beats the curated model (−0.130 [−0.215, −0.034]). Fix: state that cause; say the shared-gene view isolates the genes both represent while the union favours the drafts.
8. **"Larger" gap.** Quoted: abstract "The larger remaining gap is coverage"; Discussion heading "Coverage is the larger gap." Fix: "A second gap is coverage".
9. **Coverage gaps and Next tests.** Add that the fix uses condition-level growth data, that most gaps lack any exchange reaction, and that some carbon sources cannot be mapped; restrict the Next tests bullet to carbon sources that map to BiGG.
10. **Protocol description.** Add the runtime-environment difference and that BW25113's deletions were not applied.
11. **Missing disclosure of circularity.** iJN1463 itself was consulted during development (first curated arm in Sprint 3; quinone audits used its biomass coefficient and pathway), and U′'s CBMKr patch fixes an error iJN1463 shares.
12. **Gene count for M.** "(765 for M)".
13. **Table S1 and the Figure 1 caption.** Row 6 should end 13:07:19; the caption promises exact times the table lacks (panel declared 07:44:00, drawn 07:47:45, runs, and the H3 arms M_noIonR6 for MR1, Cola and R. palustris at 12:28:54–12:29:25, after the panel A analysis and before the 12:31:09 plan lock). Fix: correct the time and add rows, or narrow the caption to locks and downloads.
14. **Housekeeping.** Note for Tim out of date; Introduction "hand-curated models" should be singular for the draft comparison; Data availability should name the source of the gitignored iML1515 file.
15. **Reference 18.** Note that BiGG's iJN1463 is its version of the article's iJN1462.

## (c) What I could not check

- **PubMed:** captcha; Europe PMC rate-limited; references 17–18 checked against Crossref only. The Wiley page returned 403, so the iJN1462 naming rests on the project's NOTICE and the literature.
- **Current BiGG iJN1463 file:** not hashed (downloads from bigg.ucsd.edu blocked); the local file matches its card and the NOTICE.
- **External-record times:** taken as given; only their gaps after the commits were checked.
- **Effect sizes needing LP solves:** the BW25113 deletions and the CBMKr patch under this protocol. Direct effect of the deletions: 3 of the 7 genes are scored, all as "no effect" with no fitness defect.
- **GLPK build:** not recorded in either environment.
- **Writing time of the rules:** git shows only the first commit (09:40:45, after the panel was drawn at 07:47:45) and the first run (07:55:37).
- **Run-once and the 5 October agent check:** documentary only.
