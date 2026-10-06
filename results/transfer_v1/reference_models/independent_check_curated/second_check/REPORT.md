# Second independent check: iAH991, corrected iSO783, Paper 1 text, references 19–27

6 October 2026. Written by the same independent checking agent as `../REPORT.md`. It could not save files named
REPORT, so the main session saved its returned report here, verbatim apart from formatting. Scripts and outputs
are in this directory; regenerate the PDF re-read with `reparse_pdf_tables.py`.

## Summary

There are no critical problems. Every reported number for iAH991 (exploratory run) and the corrected iSO783
reproduces exactly. Four problems change wording or one rounded value, and the paper does not say that the
rebuild failed the plan's own acceptance rule.

- **Rebuild fidelity.** The agent re-read Tables S10a, S10b and S12 itself, using the drawn cell borders rather than
  the rebuild's extractor.
  - All 1,487 visible S10a rows match the model in stoichiometry and gene rules, with 0 differences.
  - The cut-off biomass row matches the S12 printing term for term (97 metabolites), and the 36 terms visible in
    S10a agree with it.
  - A second PDF text engine (poppler) agrees on 30 of 30 spot-checked reactions.
  - The manual interventions were checked against page images. ARAI is printed −1.000 in S10a and −1000 in S12.
    ATPM is printed 84.300 in S10a and 8,43 in S12, and only 8.43 reproduces Table S3: ATPM 0 gives 0.326, 8.43
    gives 0.237 and 84.3 is infeasible.
  - Recomputed validations match: glucose 0.2369, arabinan 0.2261, dextran 0.2552 and levan 0.2369, and 61
    rich-medium essential genes (published 61).
  - **Verdict:** a faithful rendering of the printed network, fit to score (but see problem 1).
- **Translation of iAH991.** It is relabelling only: 0 mismatches, and 142 growth tests are equal, including 60
  knockouts.
  - Every Btheta medium component, carbon source and energy currency is the right compound by name and formula,
    except arabinan (problem 7).
  - The energy gate passes, with and without the supplement, and there is no proton import from nothing.
- **iSO783.** Every number now in `comparison.json` matches the corrected-map recomputation exactly, and the
  re-scored matrix matches an independent rescoring with 0 call differences.
- **References 19–27.** All correct.

## Problems

1. **(Major) The iAH991 rebuild did not meet the plan's acceptance rule, and the paper does not say so.** The
   results against the plan's criteria:

   | Criterion | Result |
   |---|---|
   | Published counts | pass |
   | Table S3 growth rates | 42/45 as printed |
   | Table S1 sole carbon sources | 62/63 (PA6 fails) |
   | Table S14 knockouts | 15/17 |
   | Essential genes on glucose | 205 vs 204 |
   | Essential genes on tryptone–yeast extract–glucose | 127 vs 116 (unexplained) |

   **Fix:** state this in Methods and Limits, list the misses, and call scoring a judged deviation.

2. **(Major) The coverage conclusions rest on a medium artefact.**
   - iAH991's 0 of 25 is caused only by B12: its biomass needs adenosylcobalamin, which it cannot make.
   - The authors' own simulation media (Table S8a, and S8c "corresponding to Varel & Bryant") include
     cob(I)alamin and hydrogen sulfide.
   - With B12 at the trace rate alone, iAH991 grows in 24 of 25 conditions; the drafts grow in 14.
   - iSO783 also grows on the two unmapped dipeptides.

   **Fix:** report 0 as published and 24 with B12, and reword "Coverage did not consistently favour curation" and
   "Curated models are not a reliable fix".

3. **(Major) The union results are stated without intervals, and one value is mis-rounded.** Curated − arm on the
   union excludes zero only against M:

   | Organism | vs B0 | vs U′ | vs M |
   |---|---|---|---|
   | *P. putida* | −0.040 [−0.144, 0.059] | −0.095 [−0.203, 0.004] | −0.130 [−0.227, −0.039] |
   | *S. oneidensis* | +0.063 [−0.004, 0.133] | +0.049 [−0.011, 0.112] | +0.058 [−0.002, 0.125] |
   | *B. thetaiotaomicron* | −0.048 [−0.128, 0.035] | −0.059 [−0.138, 0.027] | −0.112 [−0.188, −0.039] |

   In the run checked, iAH991's four-way union MCC was 0.4448, which rounds to 0.44, not 0.45.

4. **(Major, wording) "What hand curation buys" and the abstract's last inference overstate.**
   - "Not better calls": no difference was detected, but the intervals allow the curated model to be better by up
     to about 0.13–0.15.
   - "Shared errors point to…" is untested. Most shared errors are genes both models call dispensable but the screen
     calls important: 416 of 710, 359 of 383, 303 of 478.
   - iAH991's extra wrong calls come mostly from one curation choice. The capsular polysaccharide in its biomass
     needs all 21 capsule genes at once, which gives 266 of its 568 curated-only "important" calls, none confirmed.
     CBLAT (BT_2760) is called essential only because the supplement is cob(I)alamin.
   - The abstract's inference that further gains will come from changing how the models predict was not tested.

5. **(Minor) Sulfide was not needed and changes nothing.**
   - With B12 only, iAH991 grows in the same 24 conditions and makes the same calls (0 of 2,718 "important" calls
     change).
   - Confirmed:
     - there is no route from sulfate or from methionine to cysteine;
     - B12 is absent from the FEBA recipe, and cysteine is 8.25 mM in it;
     - both sinks are needed for growth but supply no carbon (no growth without a carbon source);
     - the demands are export-only.
   - **Fix:** name methionine too, and change "removes the medium's sulfur source" to "limits it to a trace".

6. **(Minor) iAH991's gene map misses two loci with fitness data, the mirror image of iSO783.** BT_0823 and BT_2070
   are listed with the underscore. Effect: own MCC 0.4485 → 0.4487, and curated − U′ +0.014 → +0.013.

7. **(Minor) The arabinan override maps a different molecule representation.** iAH991's arabinan101 is a ~101-unit
   polymer (C509H812O410), while BiGG's araban__L has 3 subunits. At the same uptake that gives about 34× the carbon
   per molecule. The screen disagrees with iAH991's dispensable arabinose genes (BT_0350 −2.7, BT_0353 −2.1), and the
   override was decided after the plan. The effect on the numbers is negligible.

8. **(Minor) "Largely the same calls" needs numbers.** Agreement is 90–94% of cells (κ 0.57–0.66), but only about
   half of the "important" calls coincide (Jaccard 0.43–0.56).

9. **(Minor) The S3 misses are not "three evident misprints".** Two rates break the table's own scaling rule
   (rhamnose 12, GlcNAc-Core 2 at 2.7273). The third is an ambiguous substrate (amylopectin).

10. **(Minor) Precision.**
    - Translation decisions after the plan need acknowledging, and so does "the identifier translation had already
      listed…", which is true only for iSO783 and iGD1575.
    - "No model file is public" should be "we found none".
    - Reference 27: be specific about where iGD1575 comes from (ModelSEED).
    - The Table S12 agreement excludes 11 gene rules that S12 prints malformed.
    - L-fucose and L-rhamnose grow at only 0.0012/h; excluding them gives own MCC 0.453 and curated − U′ +0.023.

11. **Fairness check (no problem).** Given the same supplement, the drafts' calls are unchanged: 0 of 7,098 (B0),
    7,098 (U′) and 7,154 (M) calls change, and coverage is unchanged. The exploratory comparison is reportable as
    exploratory.

## References 19–27

All are correct: authors, title, journal, year, volume, pages or article number, and DOI. The statements about
iMR1_799, iLJ1162 and iGD1348 (references 22–24) are correct. Some index services could not be read: PubMed and PMC
asked for a captcha, Crossref, Semantic Scholar and OpenAlex rate-limited, and Taylor & Francis returned 403.
Publisher and repository pages were used instead.
