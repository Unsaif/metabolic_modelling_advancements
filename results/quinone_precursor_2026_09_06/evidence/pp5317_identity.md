# PP_5317 identity and evidence for quinone-precursor synthesis

Source audit, 6 September 2026. This work reads gene metadata, reference sequences and primary literature. It makes no model changes, performs no optimization and reads no numerical benchmark phenotypes. It follows the previously exposed development result; it is not independent validation of a selected repair.

**The existing PP_5317 identity is correct. Its UbiC assignment has strong homology support, but this audit does not establish direct catalytic characterization of the exact KT2440 protein. No evidence presently justifies replacing or broadening its gene–reaction rule.** A distinct, uncharacterized KT2440 protein is a research lead, with unresolved product specificity.

## Exact identity

The local gene metadata, historical GenPept map and current feature table agree on the same forward-strand interval, **NC_002947.4:6063112–6063669**, with old locus **PP_5317**, current locus **PP_RS27695**, historical protein **NP_747418.1** and current protein **WP_010955810.1**. The model already uses its unique historical alias **NP_747418_1**. These are identifiers for one gene, not interchangeable evidence for several enzymes.

Translation of the newly retrieved 558-base reference CDS, its deposited translation, both RefSeq protein sequences and **UniProt Q88C66 sequence version 1** are identical across all 185 amino acids. The initiator is genomic **GTG**, translated as methionine under bacterial genetic code 11; it is not a missing start codon. The extraction explicitly excludes a partial downstream **ubiA/PP_5318** feature also returned by NCBI. Local metadata match the retrieved reference, but no original local genomic FASTA was available for an independent original-genome sequence comparison. [NCBI protein](https://www.ncbi.nlm.nih.gov/protein/WP_010955810.1), [UniProt Q88C66](https://www.uniprot.org/uniprotkb/Q88C66/entry).

The offline [extraction script](pp5317_identity.py) reproduces the exact matches, primer checks and input hashes in [the JSON record](pp5317_identity.json). It accepts `--out` only for a fresh file. Source retrieval URLs, dates, versions and byte hashes are preserved in [the source directory](pp5317_identity.sources/). Full copyrighted papers are not redistributed.

## What supports the catalytic assignment?

Q88C66 is reviewed, yet its protein existence is **inferred from homology**. Its probable UbiC name and chorismate → 4-hydroxybenzoate + pyruvate reaction are assigned by **HAMAP MF_01632**, not a cited KT2440 enzyme assay. The saved rule, version 27, uses **E. coli UbiC P26602** as its template and explicitly retains “Probable” outside Enterobacterales. P26602 has purification and enzyme-characterization references; that establishes the template's activity, not an experiment on every homolog. [HAMAP rule](https://hamap.expasy.org/rule/MF_01632), [E. coli template](https://www.uniprot.org/uniprotkb/P26602/entry).

**Kitade et al. (2018) provide stronger, species-level evidence with an identity limit.** Their screen expressed genomic PCR products from multiple donors in wild-type *C. glutamicum*. A *P. putida* construct showed crude-extract activity (Table 1: 142 ± 12 nmol mg⁻¹ min⁻¹, n = 5). Table 5 primers 67/68 match both ends of the KT2440 CDS exactly. The NdeI cloning site would replace genomic GTG with ATG without changing initiator methionine. However, the screen does not state its *P. putida* donor strain/accession or full construct sequence. S12 is explicitly identified for a separate tolerance experiment and cannot safely be assigned as the enzyme donor. End matches therefore support PP_5317 compatibility, not exact donor or internal sequence identity. [Kitade et al., methods and Tables 1/5](https://pmc.ncbi.nlm.nih.gov/articles/PMC5835730/).

Two studies using KT2440 as a host do **not** close this gap. Yu et al. use heterologous **E. coli K-12 W3110 UbiC, CAA40681.1**. Jha et al. engineer **E. coli** UbiC; their discussion describes native PP_5317 as an annotated homolog and offers an expression hypothesis. Neither demonstrates the catalytic activity of purified KT2440 PP_5317. [Yu et al. 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5124731/), [Jha et al. 2019](https://doi.org/10.1021/acssynbio.8b00465).

## Alternatives: real chemistry, unresolved KT2440 assignment

**XanB2 is a genuine alternative enzyme in another organism.** Experiments on *Xanthomonas campestris* pv. *campestris* **XanB2/Xcc4014** support production of both 3- and 4-hydroxybenzoate and supply of the ubiquinone precursor. This does not identify a native KT2440 XanB2 or justify borrowing the Xanthomonas gene rule. [Zhou et al. 2013](https://doi.org/10.1111/mmi.12084).

An official KT2440 annotation search also retrieves **PP_3784/Q88GE0**, a distinct 383-residue protein, current locus **PP_RS19690**, at 4311281–4312432. Local metadata call it a conserved protein of unknown function. UniProt assigns an **FkbO/Hyg5-like N-terminal domain** through Pfam evidence and provides no catalytic reaction. It is not an alias or established paralog of UbiC. [Q88GE0](https://www.uniprot.org/uniprotkb/Q88GE0/entry).

That family annotation is insufficient: related FkbO, Hyg5 and XanB2 enzymes have different product outcomes despite their related folds. Experimental mechanistic work identifies specificity determinants; a broad domain hit does not tell us which product PP_3784 makes. The annotation query is also not an exhaustive homology search and cannot exclude further unannotated routes. [InterPro family description](https://www.ebi.ac.uk/interpro/entry/InterPro/IPR049368/), [Hubrich et al. 2015](https://doi.org/10.1021/jacs.5b05559).

## Discriminating next checks

1. **Source work:** recover the exact Kitade donor/construct record if available; compare PP_3784 with experimentally characterized family members using full-length alignment and validated specificity positions. This could strengthen or reject a candidate without phenotype-based reaction selection.
2. **Enzyme evidence:** test exact KT2440 PP_5317 and, separately, PP_3784 with chorismate, resolving 4-hydroxybenzoate from 3-hydroxybenzoate and other products by standards-based analysis. Include no-enzyme controls for spontaneous chorismate conversion. Activity of one protein does not establish that it substitutes for the other in cells.
3. **Native requirement:** independently constructed deletions and complementation, combined with quinone and precursor measurements, could test compensation in defined medium. A controlled 4-hydroxybenzoate rescue would test precursor limitation; it would not by itself identify a second enzyme. These are proposals, not experiments performed here.

For the next modeling decision, retain PP_5317 as the supported UbiC candidate and treat PP_3784 as unresolved source work. This audit supports neither a new alternative gene rule nor a claim that PP_5317 must be the sole biological route.
