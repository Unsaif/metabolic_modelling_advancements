# Main session's response to the second check

6 October 2026. Paper 1 doc revision 71; repository copy `docs/paper/paper1-draft.md`.

| # | Problem | Response |
|---|---|---|
| 1 | Rebuild did not meet the plan's acceptance rule | **Stated.** Results says the rebuild "did not pass every check we had set in advance". Methods lists every miss (42/45 growth rates, 62/63 carbon sources, 15/17 knockouts, essential genes 61/61, 205 vs 204 and 127 vs 116) and calls scoring a judged deviation, with the reason. Limits and the deviations list repeat it. |
| 2 | Coverage conclusions rest on the B12 artefact | **Reworded** in the abstract, Results, Discussion and Note: "the curated models grew in at least as many conditions (36/28, 9/9, 24/14), the *B. thetaiotaomicron* model only once given vitamin B12 … as published … none". The Discussion now says that curation helped on coverage. iSO783's growth on the two unmapped dipeptides is stated. |
| 3 | Union intervals; 0.44 | Union intervals added: only the differences from M exclude zero (*P. putida* −0.130 [−0.227, −0.039], *B. thetaiotaomicron* −0.112 [−0.187, −0.040]). In the re-scored iAH991 run (B12 only, corrected gene map), the four-way union MCC is 0.4451, so "0.45" is correct for the current run. |
| 4 | Overstatement | "No clear difference" now appears with the upper limits (up to 0.10–0.13). The reasons are labelled "we can only suggest reasons, untested here". The shared-error composition (416, 359 and 303 missed important genes) and the capsule choice (266 of 568 calls) are added. The abstract's untested inference is replaced by the descriptive "much of what limits these predictions is common to automatic and hand-curated models rather than specific to drafts". |
| 5 | Sulfide unnecessary | **Re-run with B12 only** (`../../Btheta/iAH991_exploratory_B12/`). The B12+sulfide run is kept in `../../superseded/`. Methods: "it cannot make cysteine from sulfate or methionine … adding sulfide raises its growth but changes none of its calls". Limits: "is then left with a trace of sulfur". |
| 6 | iAH991 gene map missed two loci | **Fixed** with `--gene-normalize identity_then_no_underscore` (identity, else without the underscore; 991 of 991 mapped) and re-scored. Own MCC is 0.4487; curated − U′ is +0.014 [−0.073, 0.102]. It is disclosed in the deviations list. |
| 7 | Arabinan representation | **Disclosed** as a deviation: a 101-sugar polymer against BiGG's three sugars, about 34× the carbon. Without the condition, curated − U′ is +0.013 [−0.074, 0.100]. |
| 8 | "Largely the same calls" | Numbers added: agreement 90–94% of calls, κ 0.57–0.66, and 43–56% of the "important" calls made by both (recomputed on the current runs). |
| 9 | "Three evident misprints" | Now: "two of the misses come from rates that break the paper's own dosing rule, one from an ambiguous substrate". |
| 10 | Precision points | All addressed: later decisions are listed as deviations; the translation-timing sentence names iSO783 and iGD1575; "we found no public model file"; reference 27 is tied to the ModelSEED reaction database, which iGD1575 uses; the 11 malformed S12 rules are noted; the fucose and rhamnose sensitivity is in Methods (+0.023 [−0.063, 0.111]). |
| 11 | Fairness check | Added to Results: "giving the drafts the same addition changes none of their calls". |
