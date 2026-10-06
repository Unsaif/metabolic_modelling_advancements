# iGD1575 (S. meliloti) was not scored: it fails the protocol's energy gate

6 October 2026. Plan: `docs/studies/transfer-v1-curated-references-plan.md` ("A model that fails the gate is
reported and not scored"). No card is written, so `scripts/compare_reference_models.py` skips this model.

## What failed

The protocol (`gembench/protocols/carbon_fitness_generic.py`, step 2c) closes every boundary reaction and asks whether
the model can still regenerate energy currencies (`gembench/checks.py`, `energy_from_nothing`). Every draft arm passed
this gate. For iGD1575 the maxima are:

| Currency | Maximal flux with every boundary closed |
|---|---|
| ATP hydrolysis | 1000 |
| NADH oxidation | 125 |
| NADPH oxidation | 125 |
| Ubiquinol-8 oxidation | 125 |

The same holds for the published SBML as read, before any identifier translation (ATP 1000 with the ModelSEED
identifiers), so the cycles are in the distributed file. The reactions involved are reversible in the file
(`reversible="true"`, bounds −1000/1000). At least four independent cycles exist, and no single reaction removal
closes them all (independent check, `../../independent_check_curated/check_energy_gate.json`):

1. Cytochrome c oxidase: the non-pumping form (rxn00058_c0) runs backwards against the pumping form (rxn10043_c0),
   and the protons drive ATP synthase (rxn10042_c0).
2. Sulfate adenylyltransferase: the ATP form (rxn00379_c0) and the GTP-coupled form (rxn09240_c0, reversed), with
   nucleoside-diphosphate kinase (rxn00237_c0, reversed), turn ADP and phosphate into ATP.
3. Glycolate oxidation: glycolate dehydrogenase (NADP, rxn00324_c0) and glycolate oxidase (quinone, rxn08655_c0,
   reversed), NADH dehydrogenase (rxn10122_c0) and transhydrogenase (rxn10125_c0, reversed) pump protons for ATP
   synthase.
4. Ribose salvage: ribokinase (rxn00772_c0) runs backwards (ADP + ribose 5-phosphate → ATP + ribose), with
   phosphopentomutase (rxn00778_c0) and nucleoside phosphorylase or hydrolase steps closing the loop.

The authors' own simulations may have constrained some of these directions in their scripts ("available from the
authors on request"); the statement here is about the SBML file as distributed.

## Other problems found on the way (would matter only if the model were scored)

- **Boundaries outside the exchanges.** The protocol closes `model.exchanges` and applies the medium, but some of the
  file's boundaries are not exchanges:
  - D-proline's boundary species sits in the cytosol compartment;
  - 21 reversible `Demand_` reactions (tRNAs, tetradecenoate) are open;
  - the glycogen and biomass boundaries are open.

  With these left open the model grows without any carbon source (independent check, problem 3).
- **A table error affecting all *S. meliloti* models.** The study's carbon-source table maps D-tagatose to `tagat__D`,
  which is not a BiGG identifier (BiGG: `tag__D`). No *S. meliloti* model can therefore grow on that condition. The
  table is locked and is not changed. The drafts lack D-tagatose anyway, so their scores are unaffected.
