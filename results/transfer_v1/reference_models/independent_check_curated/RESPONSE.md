# Main session's response to the check of the curated models outside BiGG

6 October 2026.

| # | Problem | Response |
|---|---|---|
| 1 | iSO783 gene map missed 13 genes listed as SO_nnnn | **Fixed.** `scripts/score_reference_model.py --gene-normalize identity_then_underscore` uses identity, else the underscore form (782 of 783 genes mapped). iSO783 was re-scored. The first run is kept in `../superseded/MR1/iSO783_identity_gene_map/` with a note. This is a disclosed deviation from the plan, and Paper 1 reports the corrected numbers, which match the check's: own MCC 0.546; common curated − U′ +0.051 [−0.028, 0.125]. |
| 2 | D-tagatose mapped to the non-BiGG identifier `tagat__D` | Recorded, and the locked table is not changed. It is disclosed in Paper 1 (Limits) and in `../Smeli/iGD1575/GATE_FAILURE.md`. It affects only the *S. meliloti* denominators (1 of 33 conditions can never grow), and the drafts lack D-tagatose anyway. |
| 3 | iGD1575 would grow without carbon (boundaries outside `model.exchanges`) | Recorded in the gate-failure note. The model is not scored. Before iAH991 was scored, it was checked the same way: with no carbon source it cannot grow (0), even with its two open sinks. |
| 4 | Energy gate skips `h_p` for models without a periplasm | Proton import from the medium (h_e → h_c, every boundary closed) was also run before scoring: iSO783 0, iAH991 0, iGD1575 1000. Methods will state which currencies were tested. |
| 5 | No failure record for iGD1575 | Written: `../Smeli/iGD1575/GATE_FAILURE.md`. |
| 6 | Translation equivalence check is optima-only | The reaction-by-reaction comparison in `check_translation.py` is the verification of record. The script's own check is described as a quick gate only. |
| 7 | "Share of gap closed" not in the plan | Labelled exploratory in `comparison.json` (`notes`). It is not used in Paper 1. |
| 8 | Curated-only calls use the conditions where both grow | `comparison.json` now lists the conditions. The text says "conditions where both U′ and the curated model grow". |
| 9 | MR-1 does not repeat the *P. putida* pattern | Paper 1's text is rewritten across all organisms; see the Paper 1 doc and `docs/paper/paper1-draft.md`. |
| 10 | Two dipeptide conditions unmapped although BiGG identifiers exist | Disclosed as a limitation of the locked table, with the observation that iSO783 grows on both and the drafts lack the exchanges. This is exploratory and not scored. |
| 11 | Process notes | Results committed. The absent-exchange lists printed before the plan are mentioned in Methods. |
