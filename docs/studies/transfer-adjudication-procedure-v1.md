# Blind gene-rule adjudication — procedure v1 (transfer study v1)

You are curating the gene–reaction rules of an automatically reconstructed (CarveMe) genome-scale metabolic model of one bacterium. Your decisions will be applied to the model and later scored against experimental data you will not see. **Work from biochemistry and genome annotation only.** Do not look for, open or use any phenotype data for this organism (fitness tables, essentiality screens, growth phenotypes, papers reporting mutant phenotypes), and do not open any file outside the packet directory you are given.

## What you receive

- `candidates.json` — reactions that the model needs for growth in at least one tested condition, each with its equation (identifiers and names), any EC/database annotation, its current gene rule, and for every gene: the locus tag, the Fitness Browser description, the RefSeq protein definition, and the genomic neighbourhood (three genes either side, with strand and description).
  - Type `A_geneless`: the reaction has no gene.
  - Type `B_or_rule`: the reaction has alternative genes joined by OR (possibly with AND-complexes).
- `genes.tsv` — the organism's complete gene table (locusId, sysName, scaffold, coordinates, strand, description). Search it freely.

## Decisions allowed

**A. Gene-less reactions — rule R6 (assign).** Assign a gene (or an AND-complex of subunits) only when the genome annotation names an enzyme or transporter for *exactly this reaction*: the same enzyme name or EC number, or the specific transporter for this substrate. Use operon context as supporting evidence, never as the only evidence. Abstain when the reaction is spontaneous, a diffusion step, a lumped/pseudo reaction, a biomass or bookkeeping reaction, or when only a family-level or generic annotation matches (e.g. "ABC transporter permease", "MFS transporter", "aminotransferase", "dehydrogenase", "hydrolase"). If several genes are annotated as the same enzyme and could each do it, assign them as OR alternatives.

**B. OR-rule reactions — rules R2 (join) and R1 (remove).**

- **R2 — join subunits with AND** when the alternatives are subunits of one obligate complex: explicitly annotated subunits of the same enzyme ("large/small subunit", "subunit A/B", "alpha/beta", "E1/E2/E3 component" of the same complex, carA/carB, glcD/glcE/glcF, sucC/sucD, sdhA/B/C/D …). A subunit listed alone as an alternative to its complex is removed.
- **R1 — remove an alternative** only when it is clearly not a catalyst of this reaction: a transcriptional regulator, a sensor, a ribosomal or translation protein, a "hypothetical protein" without any enzymatic annotation, a transporter listed for an enzymatic reaction (or an enzyme listed for a transport reaction), or an enzyme with a clearly different specific activity (different EC number, different substrate class). **Never remove** a same-name paralog, a gene with a family-level or broad-specificity annotation (aminotransferases, phosphatases, aldolases, dehydrogenases, esterases, acyl-CoA enzymes, transporters of the right substrate class), or a gene whose annotation is compatible with the reaction. Never apply R1 if the genes that would remain all carry generic annotations.
- Do not add new alternative genes to OR-rules in this version (no isozyme additions).

**When in doubt, abstain.** An abstention leaves the rule unchanged. Development experience behind this rule: removing apparently false alternatives from essential reactions was usually wrong, because genuinely redundant enzymes, broad-specificity enzymes and unannotated isozymes are common; joining true complex subunits and assigning a specifically annotated gene to a gene-less step were usually right.

## Output

Write `decisions.json` in the packet directory: a JSON object `{"org": ..., "procedure": "transfer-adjudication-procedure-v1", "decisions": [...]}` with **one entry per candidate**, in candidate order:

```json
{"reaction": "DHORD6", "decision": "R6",  "new_rule": "BT0892 and BT0891",
 "evidence": "BT0892 'dihydroorotate dehydrogenase' + adjacent BT0891 'dihydroorotate dehydrogenase electron transfer subunit' (PyrD/PyrK, NAD-dependent family 1B)"}
{"reaction": "ACLS", "decision": "R2", "new_rule": "BT2075 and BT2076", "evidence": "..."}
{"reaction": "PHETA1", "decision": "abstain", "new_rule": null, "evidence": "generic aminotransferases; redundancy plausible"}
```

- `decision` is one of `R6`, `R2`, `R1`, `R2+R1`, `abstain`.
- `new_rule` uses the Fitness Browser locus tags (`sysName` column of genes.tsv) joined with `and`/`or` and parentheses, and must name every gene that should remain. For a gene with no locus tag in the packet, keep its model gene id exactly as given.
- `evidence` is one or two sentences citing the annotations you relied on.

Check before finishing: every candidate has exactly one entry; every locus tag in `new_rule` exists in genes.tsv (or is a model gene id from the packet); no rule names a gene you have not seen annotated.
