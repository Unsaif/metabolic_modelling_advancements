# Superseded reference-model runs (kept for the record)

- `MR1/iSO783_identity_gene_map/`: the first iSO783 run (6 October 2026, 10:53 UTC). It mapped model genes to
  Fitness Browser sysNames by identity, as the plan said, which missed 13 genes that the Fitness Browser lists in
  the RefSeq form (SO0419 is listed as SO_0419; no MR-1 locus is listed in both forms). The independent check
  (`../independent_check_curated/`) found this. The current run (`../MR1/iSO783/`) maps by identity, else with the
  underscore (782 of 783 genes), a disclosed deviation from the plan. Effect: own MCC 0.5473 -> about 0.546;
  differences from the draft arms change by less than 0.01.
- `Btheta/iAH991_no_underscore_gene_map/` and `Btheta/iAH991_exploratory_B12_H2S/`: the first iAH991 runs
  (6 October 2026). Their gene map removed the underscore from every locus tag, which missed five loci that the
  Fitness Browser lists with it (BT_0823, BT_2070, BT_4101, BT_4253, BT_4366; two have fitness data), found by the
  second independent check (`../independent_check_curated/second_check/`). The exploratory run also added sulfide,
  which changes no call (the same check). Current runs: `../Btheta/iAH991/` (protocol) and
  `../Btheta/iAH991_exploratory_B12/` (vitamin B12 added at the trace rate only).
