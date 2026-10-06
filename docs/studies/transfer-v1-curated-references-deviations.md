# Curated models outside BiGG: deviations from the plan

Plan: `transfer-v1-curated-references-plan.md` (commit 260b92f, 6 October 2026, before any of these models was scored).
The plan file is left as committed. Each deviation below is disclosed in Paper 1 (Methods, "Deviations from the
plan").

| # | Deviation | When decided | Why | Effect |
|---|---|---|---|---|
| 1 | iSO783 genes mapped by identity, else with an underscore (`identity_then_underscore`) | After the first scoring; found by the first independent check | The Fitness Browser lists 13 MR-1 loci in the RefSeq form (SO_0419). The drafts' GenPept mapping already handled this. | Own MCC 0.547 → 0.546; differences from the drafts change by less than 0.01 |
| 2 | iAH991 genes mapped by identity, else without the underscore (`identity_then_no_underscore`) | After the first scoring; found by the second independent check | Five *B. thetaiotaomicron* loci are listed with the underscore (BT_0823) | Own MCC (exploratory) 0.4485 → 0.4487; curated − U′ +0.014 [−0.071, 0.097] → +0.014 [−0.073, 0.102] |
| 3 | iAH991 scored although the rebuild did not pass every check the plan set | Before gene-level scoring | The printed network is reproduced exactly; most mismatches trace to the paper's simulation conditions | The rebuild reproduces 42/45 growth rates, 62/63 carbon sources, 15/17 knockouts and essential genes 61/61, 205 vs 204, 127 vs 116 |
| 4 | Exploratory iAH991 run with vitamin B12 at the trace rate | After seeing 0/25 growth under the protocol; before any gene-level result | Its biomass requires adenosylcobalamin, and the screen medium has none (the bacterium grew without it) | Coverage 0 → 24 of 25. The same addition changes none of the drafts' calls. A first exploratory run also added sulfide; that changes no call and is kept in `results/transfer_v1/reference_models/superseded/` |
| 5 | Arabinan: the screen's sugar-beet arabinan mapped to iAH991's `arabinan101` (about 101 sugars), whereas the drafts get BiGG's three-sugar `araban__L` | Before gene-level scoring | Same product (Megazyme sugar-beet arabinan); the representations differ | About 34× the carbon per molecule; curated − U′ +0.013 without the condition |
| 6 | Energy gate supplemented, for models without a periplasm, by an extracellular-proton check | Before scoring iAH991 | The gate's `h_p` currency is skipped when there is no periplasm | None (0 for iSO783 and iAH991; iGD1575 already failed) |

Not deviations, but recorded for the record:

- the iGD1575 gene-rule repair and the ubiquinone overrides (in the plan's mechanism);
- the locked table's D-tagatose identifier error;
- the two dipeptides the table leaves unmapped.
