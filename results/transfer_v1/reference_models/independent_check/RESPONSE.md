# Main session's response to the check of the hand-curated comparison

6 October 2026. All numbers in the report reproduced; every wording problem was fixed in the Paper 1 doc (repository copy `docs/paper/paper1-draft.md`).

| # | Problem | Fix |
|---|---|---|
| 1 | "The one study organism with a usable hand-curated model" | Now "only *P. putida* and *M. tuberculosis* have hand-curated models in BiGG, and the pipeline cannot represent the *M. tuberculosis* screens"; published curated models of *B. thetaiotaomicron*, *S. oneidensis* and *S. meliloti* exist outside BiGG and were not scored (Results, Limits, Next tests). |
| 2 | Similarity stated as equivalence | "About the same … (differences within about ±0.1)", before and after correction; curated − U′ −0.006 [−0.108, 0.103] added; heading now "A hand-curated *P. putida* model scores about the same on shared genes but grows in more conditions". |
| 3 | Abstract last paragraph | Rewritten: the gain is small on a scale where neither drafts nor the curated models come close to the screens' reproducibility; coverage is a second gap. |
| 4 | Denominators | "Screened conditions the pipeline can map" throughout; 57 screened for *P. putida*, 321 for the twelve organisms (150 of 267 mapped = 47% of 321). |
| 5 | Unsupported ceiling | Removed. Replicate screens agree at MCC 0.90–0.94 (0.85–0.88 between experiment sets); the limit lies in the models and the binary growth call. Method described. |
| 6 | Shared disagreement | Replaced by the measured overlap (710 of 1,190 wrong calls, about eight times chance), computed again by `scripts/compare_reference_models.py`. |
| 7 | Union dismissed | Union result reported with its cause (on the 383 genes only iJN1463 contains, 1,235 of 1,287 "important" calls are unconfirmed); both views described. The main session counts 383 genes and 1,287 calls (genes absent from all three draft arms); the checker's 384 and 1,288 use one arm. |
| 8 | "Larger" gap | "A second gap is coverage". |
| 9 | Coverage fix qualifiers | Condition-level growth, missing exchange reactions (10 of 15 in *P. putida*), 14 of 57 unmappable carbon sources; Next tests restricted to carbon sources that map to BiGG. |
| 10 | Protocol description | Python 3.13 vs 3.11 and the unapplied BW25113 deletions stated in Methods. |
| 11 | Circularity | Disclosed: iJN1463 was consulted in development, and the carbamate-kinase patch fixes an error it shares. |
| 12 | M's gene count | "(765 for M)". |
| 13 | Table S1 and caption | 13:07:19; H3 exploration row added; caption narrowed to locks and downloads. |
| 14 | Housekeeping | Note for Tim updated; Introduction singular; iML1515 source named in Data availability. |
| 15 | Reference 18 | Methods: iJN1463 is BiGG's version of the article's iJN1462. |

Also corrected in conversation: an earlier remark to Tim that the benchmark's realistic ceiling is about 0.6–0.7 was wrong for the same reason as item 5.
