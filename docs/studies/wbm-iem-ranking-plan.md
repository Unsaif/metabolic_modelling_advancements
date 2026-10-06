# IEM disease ranking: plan (fixed before the prediction matrix is computed)

6 October 2026. Paper 2 extension. Written and committed before any cross-disease prediction was computed for the
analysis, and before any ranking. Later changes go in `wbm-iem-ranking-deviations.md`, with when and why.

## Question

Given the known biomarker profile of an inborn error of metabolism (IEM), do the whole-body model's predictions rank
that disease above the other simulated IEMs?

v0.3 showed that 86.5% of each IEM's own biomarker calls are correct (Harvey 1.03d). That does not show that the
calls are specific: if other diseases' knockouts produce the same changes, the predictions cannot tell diseases apart.
This study measures that specificity.

## What was seen before this plan

- **Feasibility runs.** They were used only for timing and for reproducing own biomarkers. Cross-disease values were
  compared between runs only by program (maximum difference, call agreement), never printed or examined. The runs are
  kept as feasibility records (`results/wbm_iem/feasibility/`) and are not used in the analysis.
  - **HIS with HiGHS, one IEM at a time.** Single warm-started solves took up to 29 minutes, and the run was stopped.
  - **HIS with Gurobi, one IEM at a time** (Tim's Mac).
    - All five own biomarkers reproduce v0.3: same calls, values within 6e-8.
    - It took 45 minutes for one IEM, because primal warm starts needed about 12,800 iterations on average.
  - **Gurobi, one readout at a time, timing test.** 12 IEMs × 6 readouts, described under "Solver" below.
- **HPO maps.**
  - v0.1 was written from the list of HPO terms and Harvey's metabolite list, before any cross-disease value existed.
  - v0.2 (links, term map, profiles) was written from Orphanet, HPO and Harvey's metabolite list. Only the HIS
    feasibility run existed then, and its cross-disease values were not examined.
- **Already known.** The own-biomarker calls for all 57 IEMs (v0.3, v0.4) are published in the Paper 2 draft.

## Model and protocol

These are the same as v0.3:

- Harvey 1.03d, file sha256 10e6cb6d….
- The Toolbox setup: physiological constraints and the EU average diet re-applied as runIEM_HH does, using
  `wbm_constraint_inputs_v0.3.json`.
- runIEM_HH's global constraints, with the 28 listed bile-duct exits.
- LP bounds sha256 693af5c9…, which the runner checks and records.
- Protocol `data/iem/iem_protocol_v0.2.json`: 57 IEMs, with 252 biomarker tuples (215 increased, 37 decreased).

## Readout panel

187 reactions:

- **The protocol's biomarker reactions (162).** One of them, `EX_25aics[u]`, is not in Harvey, so it is NA for every
  disease.
- **25 readouts that only the HPO profiles need.** They are listed in `data/iem/iem_ranking_extra_readouts_v0.2.txt`
  and include the 15 that the v0.1 HPO profiles needed. The main run's panel is the protocol's readouts plus this list.

If a correction after the independent check of the HPO profiles needs a readout outside the panel, it is computed in a
supplementary run with identical settings, before any ranking. Each readout's LPs are independent of which other
readouts are computed, so a supplementary run gives the same values the main run would have.

Readouts follow the protocol's convention: blood `DM_<met>[bc]`, urine `EX_<met>[u]`, CSF `DM_<met>[csf]`.

## Prediction matrix

For each of the 57 IEMs, `scripts/run_wbm_iem_cross.py` sets up the IEM exactly as the protocol does:

- its bound tweaks;
- the healthy reference flux, the maximum summed IEM flux with the auxiliary bound ±1e5;
- the healthy pin (summed IEM flux ≥ that maximum, truncated to 6 decimals);
- the disease state, with the IEM reactions fixed at zero;
- the check that the whole-body objective is still feasible in the disease state.

It then maximises every panel readout, first in the healthy state and then in the disease state. Each readout's
upper bound is raised to 1e5 for its own solve, as the protocol does for a biomarker. Calls use the protocol's rule:

- a value with |f| ≤ 1e-6 counts as 0;
- disease − healthy > 1e-6 is Increased, < −1e-6 is Decreased, and anything else is Unchanged;
- a failed solve is NA.

### Context (primary: "protocol")

The IEM's own demand sinks (its `demand_metabolites` and DM biomarkers) are open at ub 1000 throughout, as in the
protocol. Any other demand sink is closed except while it is being maximised. The IEM's own-biomarker readouts are
then the same LPs as v0.3.

### Order and solver

- **Order: one readout at a time** (`--order readout`).
  - First, the three protocol solves for every IEM (healthy reference flux, disease check, whole-body check).
  - Then, for each readout, every IEM's disease state, then every IEM's healthy state.
  - Between consecutive solves the objective stays the same and only bounds change: the IEM's reactions, its pin, its
    own sinks and its bound tweaks. Each LP is identical to the one-IEM-at-a-time order, which was checked on a small
    model with both solvers.
- **Primary solver: Gurobi 13.0.1.** It runs on Tim's MacBook with his academic licence, and Tim starts the command.
  - The protocol solves and the first solve of each readout use barrier with crossover.
  - The other solves start from the previous basis with dual simplex.
  - FeasibilityTol = OptimalityTol = 1e-7, with a time limit of 1,800 s per solve.
  - A warm solve that does not end optimal is solved again by barrier, and this is recorded.
  - About 1 in 100 warm solves is re-solved from scratch by barrier as a check (`--recheck-every 100`). They are chosen
    by a hash of IEM, state and readout.
- **Fallback: HiGHS 1.15.1, same settings and order.** It runs in the cloud and is slower.
- **On the real model.** The HIS readouts of the timing test are compared by program with the HIS feasibility run
  (one IEM at a time). We report the maximum difference and call agreement.

Command, from the repository root:

```
caffeinate -i python3 scripts/run_wbm_iem_cross.py Harvey_1_03d --backend gurobi --order readout --warm dual \
  --context protocol --extra-readouts data/iem/iem_ranking_extra_readouts_v0.2.txt --recheck-every 100 \
  --out-suffix _ranking_v1 --check-against results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json --quiet
```

The output is `results/wbm_iem/Harvey_1_03d_iem_cross_ranking_v1.json`. It is written after every readout and can be
resumed per readout under the same fingerprint.

### Checks before any ranking

Each check is reported whatever it finds.

1. **Completeness.** All 57 IEMs reach a final status. Any partial IEM and any NA readout are listed.
2. **Own biomarkers against v0.3.** The 252 own-biomarker calls are compared with v0.3 (HiGHS, interior point), and
   value differences are reported.
   - If more than 5 calls differ, we stop and trace them before ranking.
   - Any differing call is listed with its values.
3. **Rechecks.**
   - Every recheck should agree with its warm solve within 1e-6 absolute, or 1e-6 relative for values above 1.
   - If any recheck would change a call, that IEM is recomputed with barrier for every solve (`--warm ipm`), which is
     recorded as a deviation.

## Profiles

The candidates are always the 57 simulated IEMs. CRFD, GSD6, HFI, MNGIE and OROA are in the lab table but are not
simulated by the protocol, so they are neither candidates nor profiles.

| Set | Source | Profiles | Tuples |
|---|---|---|---|
| **lab (primary)** | The protocol's biomarker tuples and expected directions (215 increased, 37 decreased) | 57 | 252 |
| hpo | HPO metabolite-concentration annotations of each disorder in Orphanet, mapped to readouts (v0.2, rules below) | 39 | 142 |
| hpo_frequent | As hpo, limited to Obligate, Very frequent and Frequent annotations | 38 | 119 |
| lab_hpo_corroborated | Lab tuples whose readout and direction also appear in the disorder's hpo profile | 29 | 63 |

Files:

- **`scripts/build_iem_ranking_profiles.py`** builds the profiles. `--version v0.2` adds the Orphanet build.
- **`data/iem/iem_ranking_profiles_v0.2.json`** holds them. The v0.1 file is kept for the record.
- **One conflict.** In PC, hypoglycaemia and hyperglycaemia are both annotated, so DM_glc_D[bc] is dropped from PC's
  hpo profile.

### HPO profiles v0.2

v0.1 took its HPO annotations from the v0.1 linking tables. An independent check
(`results/wbm_iem/ranking/independent_check_inputs/`) found that this loses annotations:

- the linker kept only the first matching term per lab row;
- it dropped terms matched by another IEM with the same Orphanet code;
- its term parser missed some directed terms;
- some IEMs are linked to the wrong Orphanet entry (HLYS1, ADSL), to none (STAR, 2OAA), or to an entry without HPO
  annotations.

v0.2 is therefore built directly from Orphanet, by these rules:

1. **Sources.** Orphanet's `en_product4.xml` (JDBOR 2026-06-23, sha256 4f44e8a6…, from the Orphadata_aggregated
   repository) and HPO `hp.obo` (release 2026-09-01, sha256 93dace95…). They are kept outside the repository, with URLs
   and checksums recorded.
2. **Links.** Each IEM is linked to the Orphanet disorder (or disorders) for the same disease.
   - The v0.1 link is kept if it is right and annotated.
   - Otherwise it is replaced by the annotated Orphanet disorder for the same disease, chosen by name and gene.
   - Every link and every change is listed with its reason in `data/iem/iem_orphanet_links_v0.2.tsv`.
3. **Terms.** All of a linked disorder's annotations are candidates if they lie under "Abnormality of
   metabolism/homeostasis" (HP:0001939), or under the branches holding the two v0.1 metabolite terms outside it
   (HP:0003117 Abnormal circulating hormone concentration; HP:0040085 Abnormal circulating aldosterone concentration). Each distinct
   term is mapped to one model metabolite, biofluid and direction, or excluded with a reason, by the v0.1 rules:
   - the term's label decides when label and definition disagree, and the disagreement is noted;
   - excluded are classes of compounds, proteins, enzymes and enzyme activities, electrolytes (Na⁺, K⁺, Ca²⁺, Mg²⁺,
     Cl⁻, phosphate) and the anion gap, and groups of vitamers;
   - also excluded are acid-base states, ratios of two metabolites, terms without a direction ("abnormal …"), signs
     that are not concentrations, and erythrocyte contents, which lie outside the protocol's blood, urine and CSF
     readouts;
   - and compounds without a blood, urine or CSF form in Harvey;
   - an obsolete HPO id is replaced by its `replaced_by` term.
4. **Frequency, conflicts and corroboration.** Annotations with frequency "Excluded (0%)" are dropped. A readout with
   both directions in a profile is dropped. Corroboration of a lab tuple means the same readout and direction.

v0.2 was built and committed with this plan, before the main run.

- **Links.** 35 v0.1 links were kept, 14 changed and 3 added. One was removed: HMET's v0.1 entry is a different
  enzyme, and GNMT deficiency has no annotated entry. Four IEMs have no annotated Orphanet entry.
- **Terms.** 193 terms are in the map, 97 of them mapped. Of the 288 annotations considered, 149 are used.

A second independent check of the links, the map and the build was made before the main run, without access to any
matrix (`results/wbm_iem/ranking/independent_check_inputs/second_check/`). Its fixes are included above:

- SUCLA linked to Orphanet 17 rather than 1933;
- MMA linked to the complete-deficiency entry;
- neonatal hyperbilirubinaemia excluded as a class;
- obsolete HPO ids followed to their replacements;
- the exclusion rules listed in full.

## Scoring and ranking

- **Call.** c(d, r) ∈ {+1, 0, −1} for candidate d and readout r; NA counts as 0.
- **Score.** S(d | P) = Σ over the profile's tuples (r, s) of s · c(d, r), where s = +1 for Increased and −1 for
  Decreased. It is matches minus contradictions.
- **Rank of the true disease.** Let a be the number of candidates scoring higher and t the number tied, including
  itself. With ties broken at random:
  - the expected rank is a + (t + 1)/2;
  - the expected reciprocal rank is (1/t) Σ from i = a+1 to a+t of 1/i;
  - P(rank ≤ k) = clip((k − a)/t, 0, 1).

## Metrics and test

- **Primary metric.** The mean reciprocal rank (MRR) over the 57 lab profiles.
- **Also reported:**
  - the expected number of profiles with the true disease first (top 1) and in the top 5;
  - the median expected rank;
  - the number of profiles where no candidate scores higher than the true disease;
  - each profile's rank and the candidates scoring at least as high.
- **Null.** The true disease is replaced by a random candidate. Its expected values are exact: an MRR of
  H₅₇/57 = 0.0812, top 1 of 1/57 per profile and top 5 of 5/57 per profile.
- **Test.** One-sided Monte Carlo on the MRR, with 100,000 draws and seed 20261006.
  - p = (1 + number of draws ≥ the observed value) / (1 + 100,000).
  - The same test is applied to the top-1 and top-5 counts.
- **Strata.** Results are reported by profile size (1, 2–3, 4–6, ≥7), descriptively.

`scripts/iem_disease_ranking.py` implements all of this. It was written and tested on synthetic matrices before the
matrix existed.

## Secondary and sensitivity analyses

All of these are fixed now. Each reports the same metrics.

1. **HPO profiles.** hpo (39 profiles) and hpo_frequent (38).
2. **HPO-corroborated lab tuples.** 29 profiles.
3. **Promiscuity-adjusted score.**
   - S_adj(d | P) = S(d | P) − (n↑(P) − n↓(P)) · (p↑(d) − p↓(d)).
   - p↑ and p↓ are the fractions of d's available panel readouts called Increased and Decreased.
   - This discounts candidates whose knockout raises (or lowers) many readouts at once.
4. **Material changes only.** A call needs |disease − healthy| > 1e-3 and a relative change above 5%. The relative
   threshold is v0.3's effect-size analysis. The absolute floor removes the protocol's pin-truncation effects, which
   are of order 1e-6 to 1e-5.
5. **Confusions.** For each lab profile, the candidates scoring higher than the true disease, and those tied with it,
   are listed in full. For each candidate we count:
   - the profiles where it scores strictly higher than the true disease;
   - the profiles where it ties with a true disease whose score is positive.

   This is descriptive.
6. **Replications, if they are computed.** Neither is needed for the primary result.
   - The "minimal" context: no demand sink open except the one being maximised, for every IEM alike.
   - Harvetta 1.03d.

   For each, we report call agreement with the primary matrix and the same metrics.

## Interpretation

- **The claim.** We say the predictions rank the right disease above chance only if the primary test gives p < 0.05.
  The size of the effect is described by the MRR and the top-1 and top-5 counts against their null values.
- **Not diagnostic accuracy.** The profiles are literature biomarker lists, not patient measurements.
- **Not held out.**
  - The protocol and the lab biomarker table were used throughout model development (v0.2–v0.4), so the lab-profile
    analysis is not a held-out test.
  - The HPO profiles come from an independent annotation source, but the model and the readout panel are the same.
- **What the ranking adds.** Each true disease's own calls carry the protocol's known accuracy. The ranking adds
  whether other knockouts produce the same calls.

## Independent check

An agent that did not compute the matrix will:

- recompute the ranking metrics from the matrix with its own code;
- check the profile builder and the HPO map against the source tables;
- review the write-up.

## Files fixed with this plan

- `scripts/run_wbm_iem_cross.py`
- `gembench/wbm_iem_gurobi.py`
- `scripts/build_iem_ranking_profiles.py`
- `scripts/iem_disease_ranking.py`
- `data/iem/iem_orphanet_links_v0.2.tsv`
- `data/iem/hpo_term_readout_map_v0.2.tsv` (v0.1 kept for the record)
- `data/iem/iem_ranking_profiles_v0.2.json` (v0.1 kept for the record)
- `data/iem/iem_ranking_extra_readouts_v0.2.txt`
- tests `tests/test_wbm_iem_cross.py` and `tests/test_iem_disease_ranking.py`
