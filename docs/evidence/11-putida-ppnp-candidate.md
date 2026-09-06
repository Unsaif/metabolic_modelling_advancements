# PP_4248/PpnP candidate: sequence identity and limits of the downstream PPM assignment

Date: 2026-09-06. Status: **provisional, source-driven model hypothesis**. This follow-up was prompted by the previously observed RBK disagreement and the exact structural dependency described in the [post-run trace](../../results/quinone_repair_2026_09_06/evidence/postrun_dependency_trace.md). It is development work using an already exposed organism. No fitness values were inspected, model optimization performed, or model changes made for this evidence audit.

The identity of PP_4248 is well supported. Its assignment to nucleoside phosphorolysis comes from an enzyme-family annotation backed by experiments on an E. coli protein. A proposed downstream connection through the model's existing phosphopentomutase reaction has weaker gene-specific support. These are sufficient grounds for a transparent computational hypothesis test, not for declaring a biologically validated repair.

## Exact sequence and identifier checks

The [reproducible identity report](../../results/ppnp_repair_2026_09_06/evidence/sequence_identity.json) checks local gene metadata and the archived NCBI feature table against newly retrieved, versioned reference CDS intervals. It independently translates the complete CDS using genetic code 11 and compares the result with both the deposited GenBank translation and the separately retrieved RefSeq protein and UniProt sequences. All comparisons pass, including coordinates, strand, old locus tag, protein accession, complete terminal stop and amino-acid sequence.

| Local locus | NC_002947.4 interval, forward strand | Current locus | RefSeq protein | UniProt, sequence version | Exact protein length |
|---|---|---|---|---|---|
| PP_4248 | 4,834,852–4,835,136 | PP_RS22065 | WP_003254278.1 | Q88F51, version 1 | 94 aa |
| PP_1777 | 1,985,382–1,986,743 | PP_RS09130 | WP_010952814.1 | Q88LZ9, version 1 | 453 aa |

Local `genes.tsv` names the scaffold as unversioned AE015451. The archived feature table identifies assembly GCF_000007565.2 and reference NC_002947.4; its coordinates and old locus tags match the retrieved reference CDS. No pre-existing local genomic/CDS FASTA was available for an independent comparison with an original downloaded genome. The audit therefore establishes this precise metadata-to-reference chain, without claiming an additional original-genome comparison. Raw source bytes, normalized DNA/protein sequences, record versions, retrieval URLs and SHA256 hashes are preserved in the report and its [source folder](../../results/ppnp_repair_2026_09_06/evidence/sequence_identity.sources/).

The reference records identify the same proteins despite different annotation wording. The local PP_4248 description is an unknown UPF0345 protein; the archived feature table calls it a nucleoside phosphorylase. For PP_1777, the archived feature table says phosphomannomutase/phosphoglucomutase, the retrieved reference says phosphohexomutase domain-containing protein, and UniProt says phosphomannomutase. Sequence equality does not resolve substrate specificity.

## What supports PP_4248 as a provisional PpnP gene

[UniProt Q88F51](https://www.uniprot.org/uniprotkb/Q88F51/entry) is reviewed but labels protein existence as inferred from homology. Its gene name and adenosine/inosine phosphorolysis activities carry **ECO:0000255, HAMAP-Rule MF_01537**. Its only cited paper is the genome publication; the record supplies no direct KT2440 biochemical experiment.

The archived [HAMAP rule MF_01537](https://hamap.expasy.org/rule/MF_01537) is version 25, updated 2024-09-03. It explicitly uses **P0C037, E. coli PpnP/YaiE (b0391)** as its template. The associated family profile is version 6, with a 2017-05-10 data update and 2025-02-05 information update. Q88F51 records a match to this profile. We did not independently run the profile or infer a numerical match score.

The template's adenosine and inosine activities carry **ECO:0000269, PMID 27941785**, referring to [Sévin et al.](https://doi.org/10.1038/nmeth.4103). The independently completed [substrate evidence audit](../../results/ppnp_repair_2026_09_06/evidence/substrate_evidence.md) inspected the original supplementary experiments and identifies forward phosphorolysis support, relevant negative controls, and a conflict concerning the reverse direction. This supports provisionally testing the following forward chemistry for PP_4248, with compound protonation handled separately by the source-reaction audit:

- Adenosine + phosphate → adenine + ribose-1-phosphate.
- Inosine + phosphate → hypoxanthine + ribose-1-phosphate.

The HAMAP rule propagates reversibility, but the primary supplement's Table 9 and Figure S8 do not agree on reverse catalysis. A forward-only implementation is a conservative design choice for this test, not proof that the enzyme is irreversible. PP_4248 activity on these substrates, cellular capacity and contribution to methyl-cycle disposal in KT2440 remain homology-based hypotheses. These reactions intercept nucleosides before hydrolysis; they do **not** convert free ribose plus phosphate into ribose-1-phosphate.

## The PPM step is a separate, weaker assumption

The saved-model [metadata context](../../results/ribosyl_disposal_2026_09_06/evidence/model_context.json) contains reversible `PPM: r1p_c ⇌ r5p_c`, assigned to `NP_743933_1`, mapped to PP_1777. It also shows that PP_4248 has no existing mapped model gene. The identity audit verifies PP_1777's reference sequence, but its official [UniProt Q88LZ9](https://www.uniprot.org/uniprotkb/Q88LZ9/entry) entry is unreviewed and assigns **cpsG, phosphomannomutase EC 5.4.2.8**. Its ARBA annotation specifies mannose-1-phosphate/mannose-6-phosphate interconversion, not ribose-phosphate interconversion. Its cited papers concern the genome, rather than a PP_1777 substrate assay.

A related-enzyme precedent makes ribose-phosphate mutase activity plausible without validating this gene assignment. [Ye, Zielinski and Chakrabarty (1994)](https://doi.org/10.1128/jb.176.16.4851-4857.1994) purified **P. aeruginosa AlgC** and reported ribose-1-phosphate conversion at substantially lower activity than the preferred phosphohexose substrates. We inspected the original paper's abstract, preserved through its [PubMed record](https://pubmed.ncbi.nlm.nih.gov/8050998/); no quantitative ribose kinetic estimate is asserted here. This is another protein and species, so it cannot establish adequate PP_1777 activity or physiological flux in KT2440.

The inherited PPM edge and its GPR must remain explicitly provisional. Blocking PPM in a diagnostic model can reveal whether a proposed rescue relies on that edge. It cannot establish that PP_1777 catalyzes the reaction, and an unconstrained reaction bound cannot establish that a slow side activity has enough capacity.

## Discriminating follow-up

The strongest molecular check would assay the sequence-verified KT2440 PP_4248 protein with adenosine/inosine plus phosphate, including omission-of-phosphate and free-ribose controls, and identify ribose-1-phosphate directly. A separate PP_1777 assay should compare ribose-1-phosphate conversion with its annotated mannose-phosphate activity and measure capacity. A computational intervention can test the consequences of these explicit assumptions, with PpnP and PPM blocking controls. Recovery of an already inspected development score would show consistency with those data, not prove either enzyme assignment or independent predictive progress.

To reproduce the identity report offline from the repository root:

```sh
.venv/bin/python results/ppnp_repair_2026_09_06/evidence/sequence_identity.py
```

The script prints deterministic JSON, verifies source hashes before comparison, reads no phenotype measurements and performs no optimization. `--out` accepts a fresh output path and refuses to overwrite an existing file. Failed initial retrievals remain recorded alongside the successful later retrieval.
