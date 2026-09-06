# What the existing benchmark can measure

Only **1 of 12** curated genes covering the immediate quinone synthesis pathway and prenyl precursor support has a row in the downloaded Putida fitness export: **PP_5317**, assigned to chorismate pyruvate-lyase (`CHRPL`). All 12 occur in the genome annotation table. None of the five distinct genes assigned to the five missing terminal reactions has an exported fitness row.

| Curated gene | Reactions in scope | Exported fitness row |
|---|---|---|
| PP_5317 | CHRPL | Present |
| PP_0528 | DMATT, GRTT | Absent |
| PP_0687 | OCTDPS | Absent |
| PP_5318 | HBZOPT | Absent |
| PP_5213, PP_0548 | OPHBDC | Both absent |
| PP_5013 | OPHHX | Absent |
| PP_1765 | OHPHM, DMQMT | Absent |
| PP_5199 | OMPHHX | Absent |
| PP_5011 | OMBZLM | Absent |
| PP_0427, PP_5197 | OMMBLHX | Both absent |

The exported data contain 4,778 gene rows versus 5,661 annotated genes. The prepared baseline has 1,301 model genes, 1,300 mapped genes, and 1,047 genes shared with the fitness export. This last set exactly matches the saved baseline simulation gene axis. The audit loaded only five metadata columns from the fitness export and only gene labels from the simulation archive; it did not load numeric fitness or knockout arrays.

**Absent rows are unscored, not zero fitness and not experimental nonessentiality.** The local files have no insertion-site, barcode, strain-usage, or gene-essentiality inventory establishing why each row is absent. Actual absence from the original mutant library remains unknown. Library insertion coverage and genes with fitness estimates are separate quantities in the original method; data processing and initial abundance also affect eligibility. [Wetmore et al. 2015](https://journals.asm.org/doi/10.1128/mbio.00306-15) and the [Price et al. data documentation](https://genomics.lbl.gov/supplemental/bigfit/) support this distinction. These sources describe the method; they are not Putida-specific insertion evidence.

Adding the five missing terminal reactions with their curated GPRs cannot create measurements for their unscored genes. Later changes in aggregate development-set agreement may reveal indirect effects on covered genes, but cannot directly validate these terminal gene assignments. A row that is present also does not establish reliable measurements for every condition; per-condition finite coverage was intentionally not examined here. The deeper MEP and shikimate precursor pathways are outside this bounded coverage table.

A synthesis repair can be tested for mass balance, producibility, energy artifacts, and consistency with source reactions without claiming experimental validation. Direct biological validation would require additional pathway-specific evidence. A copied iJN1463 pathway remains evidence from that reconstruction, and agreement with the same reconstruction is not an independent experiment. Putida remains an exposed development organism.

Reproduce the metadata audit from the repository root, using a new output path:

```sh
.venv/bin/python results/quinone_repair_2026_09_06/evidence/recompute_benchmark_coverage.py --out /tmp/putida-quinone-coverage-new.json
```

`benchmark_coverage.json` records exact input hashes, reaction GPRs, annotation records, model mappings, explicit unknown/unscored states, and counts. The script refuses to overwrite an existing result. No old pipeline, manifest, or matrix was modified.
