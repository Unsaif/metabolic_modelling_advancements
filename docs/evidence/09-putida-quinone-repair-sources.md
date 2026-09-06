# Evidence brief 09 — KT2440 quinone pathway candidates

6 September 2026. This source audit supports a **provisional biochemical pathway hypothesis**, not a validated repair. It uses primary literature, the KT2440 RefSeq genome and published sequence-classification metadata. No benchmark outcome was used to choose genes, reactions or coefficients. Historical patches and frozen studies are unchanged.

## What the evidence supports

The strongest current assignment is a candidate **UbiI–UbiH–Coq7** aerobic hydroxylation system. Several curated-model gene rules retain older functional assignments. Copying the five missing reactions can test whether that representation restores net quinone synthesis; it cannot by itself establish that the copied gene rules or oxygenase chemistry are correct.

The [gene evidence JSON](../../results/quinone_repair_2026_09_06/evidence/putida_gene_evidence.json) contains 13 candidates, accession versions, source hashes, evidence grades and unresolved questions. The [RefSeq excerpts](../../results/quinone_repair_2026_09_06/evidence/putida_refseq_features.txt) retain their annotation provenance. Every selected CDS explicitly reports computational protein-homology annotation. No direct KT2440 catalytic experiment was located for these biosynthetic genes in this bounded search; this does not establish that such experiments do not exist. [NC_002947.4](https://www.ncbi.nlm.nih.gov/nuccore/NC_002947.4).

Evidence grades distinguish **A**, direct KT2440 catalytic/native-product measurement; **B**, experiments on another organism's protein; **C**, KT2440 sequence annotation or accession-specific primary computational classification; and **D**, an existing model's representation. B and C can jointly support a candidate without becoming A.

| Model locus | Current RefSeq locus / protein | Candidate role | Evidence and limit |
|---|---|---|---|
| `PP_5317` | `PP_RS27695` / `WP_010955810.1` | UbiC-family chorismate lyase; `CHRPL` | C; precursor formation |
| `PP_0528` | `PP_RS02795` / `WP_010951810.1` | IspA; `DMATT`, `GRTT` | C; precursor supply, not final tail length |
| `PP_0687` | `PP_RS03670` / `WP_010951940.1` | Polyprenyl synthetase; source-model `OCTDPS` | C; RefSeq does not specify eight versus nine isoprene units |
| `PP_5318` | `PP_RS27700` / `WP_004575206.1` | UbiA; `HBZOPT` | C; prenyl attachment, native chain length unmeasured here |
| `PP_5213` | `PP_RS27165` / `WP_003253663.1` | UbiD; `OPHBDC` | C plus B chemistry; requires prenyl-FMN |
| `PP_0548` | `PP_RS02895` / `WP_003255348.1` | UbiX; prenyl-FMN production | C plus B chemistry; not an alternative UbiD catalyst |
| `PP_5197` | `PP_RS27075` / `WP_164721831.1` | UbiI candidate; early `OPHHX` | C, primary accession-specific classification; conflicts with RefSeq's UbiF-like product name |
| `PP_1765` | `PP_RS09070` / `WP_003252662.1` | UbiG; `OHPHM` and `DMQMT` | C plus B; two SAM-dependent O-methylations |
| `PP_5199` | `PP_RS27085` / `WP_010955723.1` | UbiH; `OMPHHX` | C, including primary accession-specific classification |
| `PP_5011` | `PP_RS26140` / `WP_004576729.1` | UbiE; `OMBZLM` | C plus B; SAM-dependent C-methylation |
| `PP_0427` | `PP_RS02245` / `WP_010951761.1` | Coq7 candidate; terminal `OMMBLHX` | C plus B; does not establish redundancy with `PP_5197` |
| `PP_5012` | `PP_RS26145` / `WP_010955574.1` | UbiJ accessory factor | C; no catalytic reaction assigned here |
| `PP_5013` | `PP_RS26150` / `WP_049587492.1` | UbiB regulatory protein | C; current RefSeq calls it a regulatory kinase, not the `OPHHX` catalyst |

## Hydroxylase identity and alternatives

Kazemzadeh et al.'s primary comparative-genomics data classify the exact three accession versions above as UbiI, UbiH and Coq7. The [extracted rows](../../results/quinone_repair_2026_09_06/evidence/hydroxylase_family_rows.tsv) were decoded from the legend in `iTol_annotations/1_annot_genes_itol.txt`. The author taxon labels name *Pseudomonas* sp. KBS0802: nonredundant WP proteins can occur in several genomes. Their presence in KT2440 is established separately by its RefSeq CDS records. None of these accessions appears in the supplement's experimental-highlight list. This is specific computational support, not a KT2440 enzyme assay. [Kazemzadeh et al., 2023](https://doi.org/10.1093/molbev/msad219); [Data Set S2](https://doi.org/10.6084/m9.figshare.23230562).

The relevant biochemical precedent is the experimentally established C5 role of *E. coli* VisC, renamed UbiI. It supports reassessing `PP_5197` as the early hydroxylase and rejecting the old UbiB-only catalytic assignment. It does not prove KT2440-specific activity or rule out weak promiscuity. [Hajj Chehade et al., 2013](https://doi.org/10.1074/jbc.M113.480368).

Coq7 proteins from *P. aeruginosa* and *Thiobacillus ferrooxidans* complemented an *E. coli* hydroxylation-defective mutant. This supplies genus-adjacent experimental precedent for the terminal role proposed for `PP_0427`, rather than evidence that KT2440 has a redundant UbiF enzyme. [Stenmark et al., 2001](https://doi.org/10.1074/jbc.C100346200).

UbiL and UbiM can provide broader regioselectivity: primary experiments demonstrated two-position hydroxylation by *Rhodospirillum rubrum* UbiL and three-position hydroxylation by *Neisseria meningitidis* UbiM. Those results do not license assigning all missing KT2440 hydroxylations to an arbitrary monooxygenase. No supported KT2440 UbiL/UbiM candidate was identified here; an exhaustive genome-wide exclusion was not performed. [Pelosi et al., 2016](https://doi.org/10.1128/mSystems.00091-16).

## Decarboxylation and methylation

The `PP_5213 or PP_0548` source-model rule treats UbiD and UbiX as interchangeable decarboxylases. Primary biochemical work instead established that UbiX produces prenylated FMN, a cofactor required by UbiD-family decarboxylation. The experiment used *P. aeruginosa* UbiX, while the KT2440 accession is a homology assignment. A model rule `PP_5213 and PP_0548` would encode an inferred cofactor dependency, **not a claim that both proteins directly catalyze decarboxylation in a heteromeric enzyme**. Explicit cofactor production and use is a different representation needing its own assessment. [White et al., 2015](https://doi.org/10.1038/nature14559).

A KT2440 engineering paper discusses its native `PP_0548` UbiX homologue as a possible source of cofactor for introduced AroY. That discussion is useful strain-specific context but does not directly establish native quinone-pathway flux or UbiD dependency. [Johnson et al., 2016, Discussion](https://pmc.ncbi.nlm.nih.gov/articles/PMC5779730/).

The UbiG assignment to both O-methylations is consistent with *E. coli* experiments establishing this bifunctionality. UbiE's C-methylation role is supported by *E. coli* gene disruption, complementation and accumulated quinone intermediates. These are biochemical precedents for the KT2440 homology assignments, not measurements of its enzymes. UbiE's shared role in quinone classes also does not establish a complete native menaquinone pathway. [Hsu et al., 1996](https://doi.org/10.1021/bi9602932); [Lee et al., 1997](https://doi.org/10.1128/jb.179.5.1748-1754.1997).

## Chemistry and the two separate development arms

For the source-model hydroxylations, `substrate + 0.5 O2 → hydroxylated substrate` can conserve atoms. That does not make it a complete monooxygenase turnover. Class A flavin monooxygenases use FAD, molecular oxygen and a reducing donor. The usual net formulation for an aromatic hydroxylation includes `NAD(P)H + H+ + O2`, producing `NAD(P)+ + H2O` alongside the hydroxylated substrate. Thus a half-oxygen abstraction can omit donor demand while passing an elemental-balance test. The actual donor specificity of KT2440 UbiI/UbiH and the electron-transfer requirements of its distinct Coq7 enzyme remain unresolved here. Do not assign NADH versus NADPH by benchmark fit. [Pelosi et al., 2016, Introduction](https://doi.org/10.1128/mSystems.00091-16).

Keep two interpretations separate:

1. **Exact source transfer, D:** copy `OHPHM → OMPHHX → OMBZLM → OMMBLHX → DMQMT` from the specific local curated-model file, retaining and labeling its original equations and rules. Restoration of synthesis would show that this representation connects the pool. It would not validate its chemistry or gene assignments.
2. **Sequence-informed hypothesis, B/C:** independently declare `PP_5197` at early `OPHHX`, `PP_5199` at `OMPHHX`, `PP_0427` at terminal `OMMBLHX`, and the UbiD/UbiX dependency. UbiB's accessory role and donor stoichiometry require explicit treatment. This changes more than the five missing reactions and must not be presented as the same arm.

Tail chemistry is a separate uncertainty. Q9 was extracted from *P. putida* IAM1219, whose prenyltransferase also accepted several donor-chain lengths. This does not measure KT2440 composition. Neither `OCTDPS` nor `q8h2` in an existing model establishes its native eight-unit product. Retain the model's representation as a stated assumption pending direct strain data. [Kawahara et al., 1991](https://www.jstage.jst.go.jp/article/bbb1961/55/9/55_9_2307/_pdf).

## Small next checks

For computational development, compare declared arms using net quinone production from ordinary nutrients, atom/charge and donor accounting, and energy consistency; an artificial quinone source is only a diagnostic control. These checks establish properties of the representation, not biological truth.

To discriminate the main biological hypothesis, test the KT2440 `PP_5197` protein against defined C5- and C6-hydroxylation-defective backgrounds and measure quinone intermediates, with `PP_0427` as the C6 candidate. Native KT2440 quinone extraction with chain-resolving standards would separately determine which tail species the model should represent. These are proposed experiments, not performed results.

## Provenance

The RefSeq response is `NC_002947.4`, record date **04-FEB-2026**, downloaded 6 September 2026. The primary supplementary archive is Figshare file **42305007**, `Dataset_S2_trees+aln.zip`; only its sequence-family, taxon-label and experimental-highlight metadata support the accession-specific claims here. The JSON records full-response/archive SHA-256 values and individual metadata-member hashes. Relevant rows and RefSeq qualifiers are retained locally; genome sequence and archive remain retrievable from their recorded public sources.
