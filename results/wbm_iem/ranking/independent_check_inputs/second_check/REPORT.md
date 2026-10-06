# Second independent check: HPO profiles v0.2 (links, term map, build)

6 October 2026, before the main matrix run. Written by an independent checking agent with no access to any
cross-disease prediction. The main session saved its returned report here verbatim apart from formatting. Its scripts
are in this directory.

## Summary

**0 critical, 1 major and 7 minor problems.**

- **Sound:** the build reproduces exactly, and every term-map row names the right compound and readout.
- **Main error:** the SUCLA link points to the wrong subtype.
- **Panel:** none of the fixes needs a readout outside the 188-readout panel, except the optional glutathione mapping in
  problem 5.

## Problems

### 1. Major: SUCLA keeps a link to the wrong subtype (1933)

**Evidence.**

- The protocol names the IEM "…Lactic Acidosis, Fatal Infantile".
- Its knockout removes both the ADP-forming ligases (grRules 8802 SUCLG1 + 8803 SUCLA2) and the GDP-forming ones
  (8802 + 8801 SUCLG2). That is a loss of the shared SUCLG1 subunit.
- ORPHA:17 "Fatal infantile lactic acidosis with methylmalonic aciduria" (64 annotations) matches by name and enzyme.
  Its terms include lactate in blood and CSF, methylmalonic acid in blood and urine, 3-methylglutaconic aciduria and
  hypoglycaemia.
- ORPHA:1933 (encephalomyopathic mtDNA depletion syndrome with MMA) is the SUCLA2, ADP-only form. Its one used term is
  methylmalonic aciduria, so neither of SUCLA's lab tuples (lactate, pyruvate) is corroborated.

**Fix.** Link 17, alone or together with 1933 (both give the same result).

- Add HP:0002912 Methylmalonic acidemia → DM_HC00900[bc] Increased.
- Add exclusion rows for HP:0012087, HP:0011923, HP:0011924, HP:0008347 and HP:0001397.
- New counts: hpo 39/143, hpo_frequent 38/119, lab_hpo_corroborated 29/63 (DM_lac_L[bc] becomes corroborated).

**Outside the panel:** none.

### 2. Minor: stale counts in the plan

Under "Secondary and sensitivity analyses", items 1–2 still give hpo 29, hpo_frequent 29 and corroborated 19 (the v0.1
numbers). The Profiles table gives 39, 38 and 28.

### 3. Minor: wrong note in the map

The note on HP:0003074 Hyperglycemia says its only annotation is EF's. PC (3008) also has it (Occasional), and it causes
the PC conflict that the plan describes.

### 4. Minor: bilirubin is handled inconsistently

- HP:0002908 Conjugated hyperbilirubinemia (OTC) is excluded as a class.
- HP:0003265 Neonatal hyperbilirubinemia (PC, Occasional) means total bilirubin by its label, but is mapped to
  unconjugated bilirubin on outside knowledge. It alone brings DM_bilirub[bc] into the extra readouts.

**Fix:** exclude it (hpo 39/136, 25 extras) or record it as an explicit exception. No new readouts.

### 5. Minor: the map uses exclusion rules the plan does not list

- These are: acid-base states, terms without a direction, signs that are not concentrations, ratios, Cl⁻, and
  "erythrocyte content".
- The last rule drops HP:0034738 Reduced erythrocyte glutathione (OXOP, Very frequent), the defining finding of
  glutathione synthetase deficiency, although Harvey has gthrd[bc] and RBC_gthrd[c].

**Fix:** list these rules in the plan. Mapping the term as DM_gthrd[bc] Decreased instead **would add a readout outside
the panel**.

### 6. Minor: some link reasons are incomplete or inaccurate

None of these changes a profile unless a number is given.

- **MMA:** the reason omits the annotated mut0 entry 289916. A complete MUT knockout is mut0, the same reasoning used for
  FED. Its only used term is the same as 27's (Hyperammonemia, Occasional).
- **TETB:** "knocks out the BH4 pathway as a whole" is not what happens in Harvey. Only GTPCI (GCH1), PHETHPTOX2 (PAH)
  and TYR3MO2 (TH) exist; 7 of the 10 patterns match nothing. 238583 is still the best annotated choice.
- **2OAA:** the knockout is the mitochondrial carrier SLC25A21 (_2OXOADPTm/_2AMADPTm, Entrez 89874). 79154 is the
  DHTKD1 disease, so the link holds by name only, and the reason should say so. Dropping the link gives hpo 38/134.
- **STAR:** the knocked-out reaction is CYP11A1 (P45011A1m, gene 1583); the header id 1538 may be the same digits
  transposed. The annotated entry 168558 (complete CYP11A1 deficiency) matches the gene. Adding it gives hpo 39/138 and
  hpo_frequent 38/117, with no new readouts. This is a judgment call; record it.
- **DGK:** 329314 is the adult-onset multiple-deletion form. The infantile hepatocerebral form has no annotated entry.
  329314 has no metabolic annotations, so DGK gets no profile either way.

### 7. Minor: builder robustness and cosmetics

- Annotations that use obsolete HPO ids (5, e.g. FIGLU HP:0012335 → HP:0040087) are silently treated as outside the
  roots. This changes nothing now; the builder should follow replaced_by or raise.
- Five map rows are used by no linked disorder: HP:0000848 and HP:0000870 (outside the roots), and HP:0003073,
  HP:0003235 and HP:0011900.
- Seven map labels are Orphanet's wording rather than hp.obo's.

### 8. Minor: the supplement file is unexplained

The plan's command uses the 26-readout v0.2 extras list. `iem_ranking_supplement_readouts_v0.2.txt` (11 readouts: the
v0.2 extras not in the v0.1 extras) is listed among the plan's files but never described. State which list the main run
uses.

## Verified correct

### Links

- 57 rows: 36 kept, 13 changed, 3 added, 1 removed, 4 none.
- Every code, name and annotation count matches en_product4 (sha256 4f44e8a6…), and every "no HPO annotations" claim
  holds.
- No annotated entry exists for HYCARO, LTC4S, HYPVLI, ASNSD or GNMT.
- These are right: HLYS1→2203, ADSL→46, CYP21D, GA2, OXOP, FED, MSUD, PKU, CIT1, XAN1, PC, EP, and the removal of HMET.
- The other kept links match the reaction names of the knockouts in Harvey.

### Map

- All 97 included rows name the right Harvey compound (metNames and HMDB ids), with the biofluid and direction the label
  gives. Every readout exists and follows the readout convention.
- The exclusion facts hold in Harvey:
  - L-cystine has no blood, urine or CSF form.
  - 17-hydroxyprogesterone exists only in the endoplasmic reticulum.
  - There is no 7-biopterin and no free 3-hydroxy-3-methylglutaric acid.
  - Conjugated bilirubin is two compounds.
  - Folate and vitamin B6 are groups of vitamers.
- The label rule is applied in all three label/definition disagreements, and each is noted.

### Build

- My own parser gives 274 considered and 143 used annotations. My sets are identical to the JSON: hpo 39/137,
  hpo_frequent 38/115, corroborated 28/62, lab 57/252, and one conflict (PC DM_glc_D[bc]).
- The 26 extra readouts are complete, include the 15 v0.1 extras, and all exist in Harvey.
- Re-run on a scratch copy, the builder reproduces all three output files byte for byte.

### Roots

Outside the three roots, the linked disorders' only metabolite-like annotations are:

- brain MRS terms (PC, AGAT);
- hypoglycaemic coma (HMG, Very rare), whose readout HMG's profile already has;
- stones and anaemia.

None of them would change a profile.

### Plan

The link counts, 186/97, 274/143, the Profiles table, 188 = 162 + 26, and both checksums all match the files.

I opened no "cross" or feasibility files and used no web. The repository is unchanged; I re-ran the builder only on a
copy.
