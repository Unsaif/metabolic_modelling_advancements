# Superseded reference-model runs (kept for the record)

- `MR1/iSO783_identity_gene_map/`: the first iSO783 run (6 October 2026, 10:53 UTC). It mapped model genes to
  Fitness Browser sysNames by identity, as the plan said, which missed 13 genes that the Fitness Browser lists in
  the RefSeq form (SO0419 is listed as SO_0419; no MR-1 locus is listed in both forms). The independent check
  (`../independent_check_curated/`) found this. The current run (`../MR1/iSO783/`) maps by identity, else with the
  underscore (782 of 783 genes), a disclosed deviation from the plan. Effect: own MCC 0.5473 -> about 0.546;
  differences from the draft arms change by less than 0.01.
