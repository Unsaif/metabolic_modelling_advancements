# A specific ribosyl-disposal candidate, with an unresolved evidence gap

**Assessment, 2026-09-06:** PP_4248 is a credible candidate for a missing nucleoside phosphorolysis route. Its current annotation is supported by enzyme-family evidence, but this search did not identify a direct functional experiment on the KT2440 protein. It merits sequence and biochemical review, not immediate acceptance as a repair. No model was changed, optimized or scored in this investigation.

This is a **post-outcome** follow-up to the [dependency trace](../../results/quinone_repair_2026_09_06/evidence/postrun_dependency_trace.md). That trace proves a necessary RBK dependency in the stored model; disagreement with pooled mutant fitness does not independently establish the phenotype of a complete, isolated ribokinase deletion. The [source ledger](../../results/ribosyl_disposal_2026_09_06/evidence/source_evidence.json) records evidence categories, access limits and retrieved annotation versions. A separate [model context](../../results/ribosyl_disposal_2026_09_06/evidence/model_context.json) records the unchanged reaction network.

## The strongest candidate: PP_4248 / PpnP

The current reviewed UniProt entry **Q88F51**, explicitly assigned to **KT2440 PP_4248**, annotates a 94-residue pyrimidine/purine nucleoside phosphorylase. Its functional evidence is **HAMAP rule MF_01537, ECO:0000255**, and protein existence is inferred from homology. The listed KT2440 publication establishes genome sequence, not this catalytic activity. The local annotation still describes an unknown UPF0345 protein. [KT2440 annotation](https://www.uniprot.org/uniprotkb/Q88F51/entry).

The corresponding **E. coli PpnP/YaiE, P0C037**, has experimentally attributed activity for nucleoside phosphorolysis from Sévin and colleagues. The original study used a broad in-vitro metabolomics screen and targeted validation; its abstract was available, but its full assay supplement was not inspected here. The specific substrate assignments are traceable through the curated entry to that paper. A subsequent primary structural study describes the PpnP enzyme class and reports structures from E. coli and other organisms; its Pseudomonas structure is **P. aeruginosa**, not KT2440. [Sévin et al., 2017](https://doi.org/10.1038/nmeth.4103), [E. coli evidence record](https://www.uniprot.org/uniprotkb/P0C037/entry), [Wen et al., 2022 / deposited P. aeruginosa structure](https://www.rcsb.org/structure/7EYP).

The relevant candidate reactions are:

```text
adenosine + phosphate ⇌ adenine + ribose 1-phosphate
inosine + phosphate ⇌ hypoxanthine + ribose 1-phosphate
```

The enzyme intercepts a nucleoside **before free ribose is formed**. It is not a substitute ribose kinase: the E. coli evidence record specifically reports no formation of ribose 1-phosphate from free ribose plus phosphate. The proposed connection is therefore SAH → adenosine (existing AHCi), then candidate phosphorolysis directly or after existing ADA, then existing PPM to ribose 5-phosphate. This is a structural hypothesis, not a newly demonstrated feasible flux or a native physiological pathway.

The unchanged saved model already contains reversible **PPM**, assigned to **PP_1777**, connecting ribose 1-phosphate and ribose 5-phosphate. It has no PP_4248 model identity and no adenosine/inosine phosphorolysis reaction. Existing PNP, PYNP2r and NP1 have different substrates; reaction-name similarity is insufficient. PPM's particular gene/substrate assignment also needs evidence review before claiming that the entire alternative route is biologically established.

## SAH cleavage: distinguish the two chemistries

| Candidate or represented step | Evidence found | What it supports |
|---|---|---|
| **PP_4976 / AhcY**, represented by AHCi | Current KT2440 annotation assigns SAH hydrolase by homology. A primary transcriptomic study detects an ahcY response in **EM42**, a reduced-genome KT2440 derivative. | A supported annotation and expression context; no purified KT2440 enzyme assay was identified here. |
| **PP_3254**, represented by AHCYSNS | Local annotation says putative nucleosidase. UniProt Q88HU9 is predicted, with genome/reannotation references and no specific catalytic evidence in the retrieved record. | Insufficient evidence to establish SAH specificity or a complete native nucleosidase branch. |
| **PP_2460 / Nuh**, represented in nucleoside hydrolysis | Retrieved UniProt Q88K32 is predicted and names a ribonucleoside hydrolase. | A plausible annotation, but this search did not trace a KT2440-specific functional paper behind stronger labels in aggregated databases. |
| **Gene-less RHCYS** | The old EC 3.3.1.3 designation was transferred to EC 4.4.1.21, LuxS. Established LuxS chemistry produces homocysteine and DPD. | The stored ribose-producing reaction must not be justified by the modern LuxS name or EC number alone. No KT2440 LuxS assignment was established. |

Sources: [AhcY annotation](https://www.uniprot.org/uniprotkb/A0A140FWS3/entry), [Turlin et al., 2023](https://doi.org/10.1128/msystems.00004-23), [PP_3254 annotation](https://www.uniprot.org/uniprotkb/Q88HU9/entry), [PP_2460 annotation](https://www.uniprot.org/uniprotkb/Q88K32/entry), [transferred EC entry](https://enzyme.expasy.org/EC/3.3.1.3), [current EC reaction](https://enzyme.expasy.org/EC/4.4.1.21).

AhcY hydrolysis gives **adenosine + homocysteine** from SAH and water. An older P. putida whole-cell study reports the reverse condensation, but the accessible abstract does not establish KT2440 identity. That result cannot be attached to PP_4976 as a strain-specific assay. [1986 primary study](https://doi.org/10.1016/0168-1656(86)90020-9).

By contrast, a nucleosidase cleaves adenine off SAH, leaving S-ribosylhomocysteine. In a 1966 E. coli extract study the remaining sugar product was only tentatively identified as ribose. Later purified-protein work established the Pfs/LuxS sequence in other bacteria, with a chemically different sugar-derived product, **4,5-dihydroxy-2,3-pentanedione (DPD)**. These findings expose an annotation/chemistry concern; they do not establish that KT2440 uses LuxS or justify introducing a DPD sink. [Duerre and Miller, 1966](https://pubmed.ncbi.nlm.nih.gov/5326098/), [Schauder et al., 2001](https://doi.org/10.1046/j.1365-2958.2001.02532.x).

## Ribose disposal and transport remain incomplete explanations

A direct KT2440 study observed oxidation of supplied ribose to ribonate during anoxic electrode-driven cultivation, alongside no growth in those aldose cultures. This supports extracellular/periplasmic ribose processing under those conditions. It does **not** show export of cytoplasmic ribose, rescue of an RBK deletion, or disposal of SAH-derived carbon. The ribose experiment had one biological cultivation; the main text's explicit Gcd-knockout comparisons name other aldoses, so a ribose-specific Gcd genetic conclusion should not be inferred from the general pathway discussion. [Nguyen et al., 2021, Results and Figure 1](https://doi.org/10.1111/1751-7915.13862).

No experimentally demonstrated KT2440 cytoplasmic ribose or adenosine/inosine export route was identified in this bounded search. Ribose ABC uptake annotation and outer-membrane permeability do not establish that missing direction across the inner membrane. A KT2440 rehydration study discusses possible nutrient release after osmotic shock while measuring transcript responses; it does not establish a constitutive ribose exporter. [Rehydration study](https://doi.org/10.1186/s13213-020-01596-3).

Strain differences also matter: the 2025 study of **P. putida ATCC 17536** detected nucleoside hydrolase activity but not uridine phosphorylase activity in its extracts. This is neither KT2440 evidence nor a test of PP_4248 on purine nucleosides. [Fatima and West, 2025](https://doi.org/10.1139/cjm-2025-0161).

## The next discriminating check

1. **Resolve the candidate sequence and annotation provenance.** Match the repository's PP_4248 coding sequence to Q88F51; examine the HAMAP family assignment against experimentally characterized proteins and retrieve the original substrate-assay details. Independently check PP_1777's ribose-phosphate mutase activity. No fitness values are needed for these checks.
2. **Test the chemistry before using benchmark agreement as evidence.** Evidence distinguishing phosphate-dependent formation of ribose 1-phosphate from hydrolytic formation of free ribose would directly separate the PpnP and Nuh explanations. A specific functional attribution to PP_4248 would be stronger than a generic cell-extract activity.
3. **Resolve the SAH branch separately.** Establish whether KT2440 routes SAH through AhcY, whether PP_3254 acts on SAH, and whether any evidence supports the gene-less RHCYS reaction. Do not replace that reaction merely because a different product would remove the score disagreement.

Only after that review would a separately declared, source-supported model intervention be useful. Improved agreement on these already examined development cases would remain a diagnostic result, not independent validation of the pathway or a field-level advance.
