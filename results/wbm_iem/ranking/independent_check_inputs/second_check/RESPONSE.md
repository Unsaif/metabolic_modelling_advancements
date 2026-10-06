# Main session's response to the second inputs check

6 October 2026, before the main run. No cross-disease values had been computed beyond the HIS feasibility run, which
was not examined.

| # | Problem | Response |
|---|---|---|
| 1 | SUCLA linked to the SUCLA2 form | **Fixed.** SUCLA is linked to Orphanet 17 alone. HP:0002912 Methylmalonic acidemia → DM_HC00900[bc] Increased is added, and the five other terms are added as exclusions. |
| 2 | Stale counts in the plan | **Fixed.** |
| 3 | Hyperglycemia note | **Fixed.** The note now names PC's annotation and the conflict it causes. |
| 4 | Bilirubin | **Excluded.** By its label, neonatal hyperbilirubinaemia is total bilirubin, a class, as for conjugated hyperbilirubinaemia. DM_bilirub[bc] leaves the panel, which is now 187 readouts (162 + 25). |
| 5 | Exclusion rules not listed; erythrocyte glutathione | **Rules listed** in the plan. The erythrocyte-content exclusion stands: the protocol's readouts are blood plasma, urine and CSF. |
| 6 | Link reasons | **MMA** relinked to the complete-deficiency entry 289916 (mut0), as for FED; the profile is unchanged. Reasons rewritten for **TETB** (only GCH1, PAH and TH reactions exist), **2OAA** (linked by name; the knockout is the SLC25A21 carrier), **STAR** (CYP11A1 is the modelling proxy; 168558 not added, recorded as a judgment) and **DGK**. |
| 7 | Obsolete HPO ids; unused rows; label wording | The builder now follows `replaced_by` and raises on an obsolete id without a replacement. HP:0040087 (reached from the obsolete HP:0012335) is added as an exclusion. Unused rows are harmless and kept. The labels are left in Orphanet's wording; the map is keyed by HPO id. |
| 8 | Supplement file | **Removed.** The plan states that the main run's panel is the protocol's readouts plus `iem_ranking_extra_readouts_v0.2.txt`. |

New counts after the fixes:

| | Profiles | Tuples |
|---|---|---|
| hpo | 39 | 142 |
| hpo_frequent | 38 | 119 |
| lab_hpo_corroborated | 29 | 63 |

There are 25 extra readouts. 288 annotations were considered, of which 149 are used.
