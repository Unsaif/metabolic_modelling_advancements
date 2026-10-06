# Independent check of the ranking inputs (HPO map, profile builder, ranking script)

6 October 2026, before the prediction matrix was computed. Written by an independent checking agent. It returned the
report as a message, and the main session saved it here verbatim apart from formatting. Its scripts are in this
directory. They read the Orphanet and HPO files from the session's scratch space; see `RESPONSE.md` for sources.

## Summary

**0 critical, 4 major, 4 minor.**

- **Sound:** the primary analysis. The lab profiles (57 profiles, 252 tuples) and the ranking code (scores, ties,
  metrics, Monte Carlo) reproduce exactly.
- **Problems:** in the HPO-based secondary profile sets, which inherit errors from the v0.1 Orphanet/HPO linking, and
  in one descriptive analysis the script does not implement.
- **Decide before the run:** items 2 and 6a would change the readout panel.

Sources I used: the Orphanet phenotype file en_product4 (dated 2026-06-23, GitHub mirror; it agrees with every other
annotation in the repository tables) and hp.obo release 2026-09-01.

## Problems

### 1. Major: four HPO annotations are silently lost

The builder reassembles each disorder's HPO annotations from `hpo_evidence` plus
`hpo_metabolite_tuples_not_in_lab_table.tsv`. Two things in the v0.1 linker (`scripts/link_iem_ground_truth.py`) make
that lossy:

- it records only the first term whose name substring-matches a lab row;
- it leaves out of the not-in-lab file any term that matches a lab row of the same Orphanet code, even when that row
  belongs to another IEM.

What is lost:

- **DPYR:** HP:6000118 (urinary dihydrouracil up) and HP:6000331 (urinary thymine up), both Frequent. The EX_56dura[u]
  and EX_thym[u] rows recorded Uraciluria and dihydrothymine as their evidence instead.
- **HLYS1:** Citrullinuria and Elevated plasma citrulline (Occasional). These went only to HLYS2, which shares Orphanet
  code 3124.

The correct counts are hpo 29/83 and hpo_frequent 29/68, with no new readouts.

**Fix:** build the HPO sets directly from the Orphanet annotations of each linked code, and add map rows
HP:6000118 → EX_56dura[u] Increased and HP:6000331 → EX_thym[u] Increased.

### 2. Major, time-critical: some v0.1 Orphanet links are wrong or missing

- **HLYS1** ("Hyperlysinemia I, familial") is linked to 3124 Saccharopinuria rather than 2203 Hyperlysinemia. 2203 has
  9 directed metabolite terms, including argininuria, low ornithine, low urine alpha-ketoglutarate, high CSF lysine and
  low CSF arginine.
- **STAR** is unlinked because its header Entrez id (1538) does not map. Orphanet 90790 has decreased cortisol
  (Frequent) and hypoglycaemia (Occasional).
- **2OAA** is unlinked because its header gene is SLC25A21. Orphanet 79154 has three: alpha-aminoadipic aciduria (Very
  frequent), 2-hydroxyadipic aciduria (Frequent) and raised circulating 2-aminoadipate (Frequent).
- **ADSL** is linked to 31 Oxoglutaric aciduria; this has no effect at present.

Relinking needs readouts the panel lacks (DM_arg_L[csf], DM_lys_L[csf], EX_L2aadp[u]). Either decide before the
matrix run, or state "as linked in v0.1" as a limitation.

### 3. Major: the lab_hpo_corroborated set rests on substring matching

I counted a lab tuple as corroborated only when an HPO annotation maps, through the map, to the same readout and
direction.

- **Four committed tuples are not corroborated:**
  - GA1 EX_c5dc[u]: glutarylcarnitine is not glutaric acid.
  - HYPRO1 EX_4hpro_LT[u]: hydroxyproline is not proline.
  - IVA EX_3ivcrn[u]: the carnitine ester is not 3-hydroxyisovaleric acid.
  - HPC EX_C05770[u]: it rests on Porphyrinuria, a class term the map itself excludes.
- **Six true corroborations are missed:** HLYS1 and HLYS2 EX_lys_L[u] (Hyperlysinuria), LNS DM_urate[bc], NAGS and OTC
  DM_nh4[bc], and ASA EX_argsuc[u]. The causes are name mismatches: "Urate" vs uric acid, "Ammonium" vs ammonia, the
  VMH spelling "Arginosuccinic", and "Hyperlysinuria" being parsed as "hyperlysin".

The set should be 18 profiles/39 tuples, not 19/37.

**Fix:** define corroboration by readout identity, which is also how the hpo sets are built.

### 4. Major (fidelity to the plan): the confusion analysis is not implemented

Secondary analysis 5 is not implemented, and the list of candidates scoring at least as high is cut to 10. The script
keeps only `better[:10]`, with ties ordered by call_index. Profiles where every candidate scores 0 have 56 tied
candidates, so any confusion count made from the output would be biased.

**Fix:** output the full list and implement the count.

### 5. Minor: the v0.1 term parser misses directed metabolite terms

These are HMG nonketotic hypoglycaemia (Very frequent) and recurrent hypoglycaemia (Frequent), GA1 fasting
hypoglycaemia (Occasional), and IVA lactic acidosis (Frequent). With item 1, adding them gives hpo 86 and hpo_frequent
70 tuples, with no new readouts. Otherwise, document this as a limitation.

### 6. Minor: map judgment calls where the HPO definition and label disagree

- **a. HP:0025436 (11-deoxycortisol).** The definition describes deoxycorticosterone, "a precursor to aldosterone",
  which fits CMO1 better. Switching would replace the extra readout DM_11docrtsl[bc] with DM_11docrtstrn[bc], which is
  already in the panel, so decide before the run.
- **b. HP:0032164 (blood folate).** The definition says "folic acid", and DM_fol[bc] is already in the panel. Mapping
  it would add a FIGLU profile. The vitamer exclusion is defensible but should cite this.
- **c. HP:0002160 (hyperhomocystinemia).** Mapping to Lhcystin follows the definition. HMET's own protocol biomarkers
  use hcys_L, which is closer to the clinical total-homocysteine measurement. Low impact.
- **d. HP:0003344 (3-methylglutaric aciduria).** The definition says 3-hydroxy-3-methylglutaric acid, which has no free
  form in Harvey, so mapping by the label is right.

### 7. Minor: wording of the exclusion rules

- "Inorganic ions excluded" conflicts with including ammonium via Hyperammonemia, which appears in 12 of the 29 hpo
  profiles. Say "electrolytes (K⁺, Na⁺, Ca²⁺)" instead.
- L-cystine also occurs in Harvey's gut lumen and faeces, not only in organ cytosols. The conclusion (no blood or urine
  form) still holds.

### 8. Minor: robustness and cosmetics

- `build_iem_ranking_profiles.py` skips any evidence string it cannot parse without saying so; it should raise. All 51
  strings parse today.
- The builder keys corroboration by (IEM, reaction) and ignores call_index. MMA has two calls in the lab table; no
  effect now.
- The ranking script's default output name becomes `Harvey_1_03d_ranking_ranking_v1.json`, not the name in its
  docstring.

## Verified correct

- **HPO map.**
  - All 55 included rows name the right compound by Harvey's metNames and HMDB ids, with the biofluid and direction the
    term implies.
  - Every readout follows the convention: EX_ reactions are in rxns, and every DM_ metabolite is in mets.
  - The only duplicate Harvey id is 3hivac/CE2028. Both are connected, and 3hivac is the protocol's own id.
  - The exclusions are applied consistently, and the 77 map rows match the 77 annotated terms exactly.
- **Profile builder.** My own builder reproduces all four sets and the 15 extra readouts exactly: lab 57/252
  (215 increased, 37 decreased), corroborated 19/37, hpo 29/79, hpo_frequent 29/66.
  - There are no direction conflicts, and only the 57 protocol IEMs are used.
  - The "Excluded (0%)" rule drops only EF and TYR3, and frequencies agree between the two source files.
  - Each evidence string holds exactly one term.
  - The duplicated HLYS rows are deduplicated, and EX_25aics[u] is the only readout absent from Harvey.
- **Ranking script.**
  - **Synthetic check.** I compared an exact-fraction reimplementation with the script on 60 synthetic matrices
    (57 candidates × 177 readouts with NA and absent readouts, plus small cases with many ties). They agree on every
    rank statistic (a, t, expected rank, expected reciprocal rank, P(top-1) and P(top-5)), on plain and adjusted scores
    under both call rules, and on the summaries and strata. The maximum difference is 1e-16.
  - **Monte Carlo.** The test is unbiased. On the script's own draws, its float tie tolerance classifies every draw
    exactly as integer arithmetic does. The p-values match full enumeration within Monte Carlo error, the null mean is
    H_n/n, and the fixed seed reproduces the result.
  - **Edge cases.**
    - A profile that is NA for every candidate gives expected rank 29 and reciprocal rank H₅₇/57.
    - Floating-point near-ties of 1e-16 count as ties. Distinct adjusted scores always differ by at least 1/177², so the
      1e-9 tolerance never merges them.
    - The material rule is strict at both 1e-3 and 5%, and uses (d−h)/max(|h|,|d|), as v0.3's effect-size analysis did.
    - The field names match the runner's output.

I opened no "cross" files and did not run the runner. The repository is unchanged: one git-ignored .pyc file that my
import created has been deleted.
