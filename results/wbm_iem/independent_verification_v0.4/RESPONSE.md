# Main session's response to the v0.4 verification

6 October 2026, Claude (Opus 5.5). Each discrepancy in [REPORT.md](REPORT.md) and what was done. "Paper 2" means the Claude Doc (working copy) and its repository export `docs/paper/paper2-draft.md`.

| # | Discrepancy | Action |
|---|---|---|
| 1 | Published accuracies swapped | Corrected to Harvey 84.9% (214/252) and Harvetta 85.3% (215/252) in Paper 2 (Introduction, Results, Discussion), the results note and the roadmap. The frozen plan stays as frozen; the swap is recorded as Deviation 4 in the results note. |
| 2 | "Every run finished on 6 October" | Paper 2's note now says the runs finished on 5 and 6 October; the results note gives each finish time. |
| 3 | Port reproduction overstated | Paper 2 now says "every stored physiological-constraint bound … except three blood–brain-barrier uptakes. The current diet code changes 41 further bounds." |
| 4 | Verification described before it happened | Paper 2's Methods now describe the completed verification and point to the reports, including the one figure that could not be checked. |
| 5 | "All five below 0.5%" | Paper 2 (Abstract, Results) now says four of the five; the fifth changed because a healthy maximum fell to zero. |
| 6 | Caps presented as established physiological caps | "Capped" is now defined as equal positive maxima. Paper 2 and the results note state that 6 of 18 (Harvey) and 7 of 16 (Harvetta) equal a bound on the same metabolite and that the rest were not traced. The uracil example is removed, and the v0.3 results note carries a correction. The Abstract and the Discussion heading no longer say "physiological cap". |
| 7 | R1 transitions incomplete | Results note table and readings now report 1 and 3 errors turned into other errors, with the biomarkers named; Paper 2 says both rules do this. |
| 8 | Harvetta results existed before the freeze | Paper 2 now says "before any Harvetta IEM result (or biomarker minimum) existed". |
| 9 | Figure 1, point 5 | Now "the same positive maximum in both states". |
| 10 | "Always" | Removed: "In those cases the model can send the metabolite elsewhere". |
| 11 | "Differ only because" | Now "differ mainly because …; its formula also divides by all 252 biomarkers, including one absent from the models." |
| 12 | ASNSD wording | Paper 2 and the results note now give Python's 1.3 × 10⁻⁶ and MATLAB's 0 for Harvetta, and say the Harvey value of 131 is Python's alone (MATLAB had no optimum). |
| 13 | Mechanism stated as fact | Marked as an interpretation in the Abstract, Results and Discussion of Paper 2 and in the results note. |
| 14 | Gurobi version not recorded | The version (12.0.0, with MATLAB R2024b Update 9 and COBRA 67c790d) comes from the job runner's own verification record, now copied to `tools/matlab/runner/RUNNER_VERIFICATION_2026-10-05.md` and cited as the source. The job records themselves still name only "gurobi". |
| 15 | Ref 4 does not cover IPX | Reference 6 added (Schork & Gondzio 2020, doi:10.1007/s12532-020-00181-8) and cited for the interior-point solver. Its volume and pages are on Tim's checklist. |
| 16 | Recon3D denominator | Paper 2 now says "50.2% of the 205 biomarkers that … Recon3D could represent". |
| 17 | Deviations list incomplete | Deviation 4 added to the results note. |
| 18 | Wrong cross-reference | "below" → "above". |
| 19 | MATLAB run count | Paper 2 now says "MATLAB completed one run per model". |
| 20 | Certification wording | Paper 2 now says a re-run of these solves was certified and reproduced the reported optima within 8.8 × 10⁻⁵. |
| 21 | Omissions | `v04_report.py` and `iem_error_anatomy.py` added to Paper 2's code table; results by biofluid added to the results note. |

Not resolved here:

- **Shlomi et al. 2009, recall 0.10.** The main session's own attempts (link.springer.com and Europe PMC) were also rate-limited. It is on Tim's checklist in Paper 2.
- **The pre-MATLAB independent re-implementation (v0.3).** Its outputs are not in the repository. Paper 2 still mentions it; the claim rests on the v0.3 plan's text.
- **Setup bounds on the 430 entries the global step overwrites.** These rest on the study's own comparison file.
