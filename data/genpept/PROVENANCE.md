# GenPept-derived gene identifier maps

Purpose: map EMBL GEM (CarveMe) gene identifiers, which are RefSeq protein accessions, to the locus tags used by the
Fitness Browser.

Source: NCBI E-utilities efetch (db=protein, rettype=gp, retmode=text), queried 5 September 2026 through the Claude
desktop browser pane, for every protein accession in the EMBL GEM of each organism. NCBI records are public domain.

- `Btheta_batch{0,1,2}.gp` — raw GenPept records for Bacteroides thetaiotaomicron VPI-5482 (674 records, three GET batches).
- `<org>_genpept_map.tsv` — one line per record, parsed in the browser (tools/genpept_parse.js):
  `version  locus_tags  old_locus_tags  gene_names  coded_by  definition` (multiple values ';'-separated).
  Btheta (674), Putida = Pseudomonas putida KT2440 (1291), MR1 = Shewanella oneidensis MR-1 (894),
  Smeli = Sinorhizobium meliloti 1021 (1195). Raw records for Putida/MR1/Smeli were parsed on the fly and not kept.

Locus tags are normalised to Fitness Browser sysNames by removing underscores where the Browser omits them
(BT_0009 -> BT0009, SO_0020 -> SO0020; PP_0001 and SMc02791 are unchanged).
