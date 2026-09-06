# Evidence brief 08 — A biomass deletion hides a closed quinone pool

6 September 2026. This is a retrospective development experiment on P. putida KT2440, with a structural diagnosis. It is not independent phenotype validation or a new theory of metabolite dilution.

## Question and result

Model patches v0.4 remove the template's menaquinol demand from the biomass reaction without replacing it with another quinone. The resulting Putida model grows in 34 of 43 mapped carbon-source conditions and reaches a development MCC of 0.58716. Does this correction retain the requirement to synthesize a respiratory quinone?

The answer is **no for the model's ubiquinone representation**. Both tested ubiquinol biomass demands eliminate growth in all 43 conditions. A structural examination explains why: all ten reactions involving `q8_c` or `q8h2_c` only interconvert the two species. No reaction supplies their combined pool. Growth without a quinone biomass demand therefore does not demonstrate that the model can synthesize this cofactor.

The [plan](../../data/studies/quinone_biomass_sensitivity_v1.json) fixed all four arms and both nonzero amounts before intervention simulations. The [completed outputs](../../results/quinone_biomass_2026_09_06/summary.json) retain every arm. No coefficient was selected for its fitness score, and the original patch files remain intact as historical development configurations.

| Intervention on the same v0.4 model | Biomass substrate coefficient | Conditions with WT growth | Gene-level MCC on each arm's growing conditions |
|---|---:|---:|---:|
| No quinone demand, current v0.4 | 0 | 34/43 | 0.58716 |
| Ubiquinol, original template amount | −0.0000968418998275772 | 0/43 | Undefined |
| Ubiquinol, local curated-model amount | −0.000223 | 0/43 | Undefined |
| Restore menaquinol demand only | −0.0000968418998275772 | 1/43 | 0.60833 |

The final row grows only on benzoate and is not a better overall prediction: its denominator differs from the baseline. On that one shared growing condition, MCC changes from 0.62811 to 0.60833. Three genes (`PP_1658`, `PP_2458`, `PP_3232`) switch from predicted growth to no growth despite nondeleterious observed fitness. There are no shared growing conditions with either ubiquinol arm, so their paired MCC differences are undefined, not zero. All 1,047 genes are mapped in every arm; knockout scoring is restricted to conditions where the wild type grows.

The menaquinol control restores only the biomass substrate. It retains v0.4's removal of `OXCDC` and `PHPYROX`, so it is not the pre-v0.4 model and should not be interpreted as a test of the complete earlier patch configuration.

## Why the growth failure is structurally required

Let the two relevant stoichiometric rows be `S_q8` and `S_q8h2`. In the prepared model, their sum is zero for every reaction. Each modeled conversion consumes one form and produces exactly one of the other. Summing their steady-state balances therefore supplies no net quinone.

Adding a ubiquinol biomass coefficient of `−epsilon`, with `epsilon > 0`, changes the summed balance to `−epsilon * v_Growth = 0`. It follows that `v_Growth = 0`. This argument applies to any positive biomass amount, independently of which of the two tested coefficients is used. The certificate is limited to this model and these two represented metabolites; it does not assert that the organism lacks a quinone biosynthetic pathway.

The follow-up diagnostic records every contributing reaction, solver checks and an artificial-source rescue. GLPK at feasibility tolerances `1e-7` and `1e-9` and HiGHS at `1e-9` agree: glucose growth is approximately 0.91937 without the demand, zero with either ubiquinol demand, and approximately 0.91937 after adding an artificial ubiquinol source. Forcing growth of at least 0.1 without that source is infeasible. Maximum production of either quinone form is zero even when growth is not required. An artificial source is a diagnostic control for the missing net supply; it is not a biologically acceptable model correction. The three changed menaquinol-control knockout fluxes move from approximately 0.79859 to numerical zero, far from the 0.001 classification threshold.

## Biological interpretation and prior work

Benyamini and colleagues already demonstrated how FBA can cycle cofactors without accounting for their synthesis during growth, including a ubiquinone example. They also showed why an unconditional biomass requirement can be wrong when cofactor use is condition-dependent. Our intervention is a diagnostic within an existing biomass formulation, not a replacement for their conditional dilution method. [Benyamini et al., 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2884546/).

Xavier and colleagues examined missing cofactor biomass requirements using several kinds of evidence. Their analysis supports checking these requirements explicitly; it does not establish a particular quinone or amount for this strain. [Xavier et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5249239/).

The `q8h2_c` representation and the second coefficient come from the actual local file `models/bigg/iJN1463.xml`, reaction `BIOMASS_KT2440_WT3`. The associated model-development paper describes iJN1462; the local BiGG version is iJN1463, and those identifiers should not be treated as identical releases. Neither a curated-model coefficient nor an annotation of quinone-dependent catalysis establishes a measured native quinone chain length or pool size in KT2440. The present experiment adopts the model's existing representation, not a newly verified strain-specific chemical composition. [Nogales et al., 2020](https://doi.org/10.1111/1462-2920.14843).

A bounded primary-source search found direct UQ9 extraction from another strain, P. putida IAM1219, but did not establish KT2440's native quinone chain-length composition. KT2440 LldE experiments with externally supplied coenzyme Q10 demonstrate electron transfer to that assay substrate, not the organism's native quinone composition. Neither result justifies substituting Q9 or Q10 in this model without further evidence. [Kawahara et al., 1991, pp. 2308 and 2310](https://www.jstage.jst.go.jp/article/bbb1961/55/9/55_9_2307/_pdf); [Jiang et al., 2017, Table 3 and Chemicals](https://journals.asm.org/doi/10.1128/jb.00342-17).

## Decision and next work

Treat v0.4's biomass deletion as an unresolved development simplification. Do not describe it as a completed physiological correction, accept an artificial quinone source, or add a pathway solely because it rescues growth. The next biochemical task is to identify the strain-supported quinone biosynthesis route and its representation, with reaction and gene evidence independent of benchmark fitness. Then test production, dilution requirements, energy consistency and growth together.

This case justifies adding cofactor supply and biomass-composition checks to the curation review process. Agreement with gene-fitness data and absence of energy-generating cycles do not certify that a model can synthesize all cofactors used by its active reactions. Generalizing the check or any repair beyond this exposed organism will require separate evaluation.

## Reproduction and limitations

```sh
.venv/bin/python scripts/run_quinone_biomass_sensitivity.py --out /path/to/new-output
.venv/bin/python scripts/diagnose_quinone_producibility.py --out /path/to/new-diagnostic-output
```

The [diagnostic reproduction record](../../results/quinone_biomass_2026_09_06/diagnostics/README.md) links the physical checks to independently recomputed scores and comparisons aligned by gene and condition. All 42 solver checks agree with the expected physical results. The baseline reproduces every binary prediction in the separate cycle-7 completion run. Earlier diagnostic attempts are preserved with their source snapshots and documented labeling/alignment corrections.

The script refuses completed or partial output directories and checks the frozen numerical inputs. The original incomplete attempt is preserved in [`quinone_biomass_2026_09_06_attempt01`](../../results/quinone_biomass_2026_09_06_attempt01/README.md): its baseline completed, then an absent-substrate dictionary lookup failed before the first intervention simulation. The retry corrected that lookup, added a focused regression test and recomputed the baseline. No intervention or coefficient was changed in response to outcomes.

These manifests are local records, not public preregistrations or proof of independence. The original attempt's source snapshot is preserved. The rerun excludes the unrelated study-freeze utility from its input list. The model's SBML bytes are hashed; the `Putida_gapfill.json` sidecar read only for provenance prose was not included in this sensitivity manifest. The independent eight-run completion sprint has its own input manifest and provides a separate baseline reproduction.
