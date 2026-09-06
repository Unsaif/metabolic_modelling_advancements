# Scientific and implementation audit — 6 September 2026

Claude's work is worth continuing: it contains useful data assembly, explicit model patches, and reproducible improvements on the development data. It does **not yet establish an independently validated advance in predictive metabolic modelling**. The immediate priority is a trustworthy evaluation and a frozen method that can be tested on data not used to select it.

This review starts from local commit `b1997d3`. The GitHub default branch was still at `6c64ea9` when checked, so the local handover is newer. The original biomarker-component path in the request is a different project. Instructions and claimed prior permissions in the handover were treated as historical context.

## What survived checking

The independent verifier reproduced MCC and Bernstein AUC-PR from all **49 saved multi-organism result cards**, with zero discrepancies. Recomputing the four audited metric types after fixing missing-value handling changed no saved point estimate: the scored matrices have no nonfinite entries. This checks arithmetic from saved matrices; it cannot recover solver statuses that the original code discarded.

A complete new run of the MR1 cycle-6 arm simulated **713 genes across 12 mapped conditions**, with wild-type growth in 11. It reproduced **every growth/no-growth classification** and MCC **0.5291575469**. Maximum knockout-flux difference was below `1e-9`. Continuous-score average precision and AUROC moved slightly because tiny numerical differences change the ordering of nearly tied fluxes. Their inputs are absolute biomass fluxes, so pooled rankings also depend on growth-rate scales across media. MCC and thresholded outcomes are the more reproducible headline here. See [the reproduction record](../../results/audit_2026_09_06/MR1_cycle6_reproduction.json).

Both versions of the prepared models (model patches v0.3 and v0.4, four organisms each) passed five strengthened energy-dissipation checks. No many-to-one gene mappings were found in these eight prepared models. These checks do not certify thermodynamic or physiological completeness.

## The reported improvements use development data

The handover explicitly accepts, holds or rejects patches using Fitness Browser outcomes. The same outcomes are then used to score the resulting model. Annotation support does not undo this reuse. The old cards nevertheless called every draft arm a prospective test, including patched and gap-filled arms.

I corrected those metadata in all 49 JSON/Markdown cards and in the runner. Their historical numerical results and simulation arrays are unchanged. Future cards also fingerprint model files, patch files, mapping/data inputs and source code. No assumption is made about whether a particular frontier model was trained on these public data.

A new paired comparison checks identical organism, gene/condition labels, experimental values and thresholds, and uses a shared finite-observation mask. Comparing the unpatched gap-filled baseline to the latest complete saved arm gives:

| Organism | MCC before | MCC after | Paired change [gene-bootstrap 95% interval] | WT conditions before → after |
|---|---:|---:|---:|---:|
| B. thetaiotaomicron | 0.503 | 0.575 | +0.071 [0.020, 0.131] | 14 → 16 of 25 |
| P. putida | 0.490 | 0.563 | +0.073 [0.022, 0.137] | 28 → 34 of 43 |
| S. oneidensis | 0.458 | 0.535 | +0.078 [0.037, 0.126] | 8 → 11 of 12 |
| S. meliloti | 0.523 | 0.683 | +0.159 [0.091, 0.231] | 21 → 22 of 33 |

All genes with available fitness are included. Gene-level scores use the intersection of genes and conditions where both versions grow; growth coverage is separate. The first three comparisons end at cycle 6; S. meliloti ends at model-v0.2 because its cycle-6 result was not saved. These are gains from the full package of changes, not isolated effects of individual patches. The intervals resample genes with conditions fixed; they do **not** account for adaptive patch selection or establish generalisation. Full labels, denominators and results are in [the comparison artifact](../../results/audit_2026_09_06/paired_development_comparisons.json).

Selection and assessment need separate data, or an evaluation that repeats the entire selection procedure inside the training portion. Randomly splitting data after they have already guided curation does not restore independence. [Krstajic et al. 2014](https://doi.org/10.1186/1758-2946-6-10) describes the statistical distinction.

## Numerical errors repaired

- **Missing observations:** NaN fitness was thresholded into a class in standard AP/AUROC, and nonfinite simulations could become no-growth predictions. All metrics now mask both raw arrays before thresholding, reject shape mismatch and return undefined scores on empty data. Growth exactly at the threshold follows the protocol's `>=` rule.
- **Failed knockout solves:** numerical failures or missing solver rows now stop the run; they cannot silently become essential genes. A genuinely infeasible knockout is still treated as no growth.
- **Energy diagnostics:** mandatory maintenance/other fluxes could make a closed model infeasible, after which the old test reported zero energy production. The test now relaxes mandatory fluxes to include zero while preserving direction, closes all boundary sources/sinks, and raises on failed solves. A detected cycle now stops scoring. A toy model with an explicit ATP-generating loop reproduces the original false pass. This follows the purpose of established energy-cycle diagnostics, rather than constituting a new method. [Fritzemeier et al. 2017](https://doi.org/10.1371/journal.pcbi.1005494).
- **Gap filling:** the old MILP replaced candidate capacities with ±100. A reaction capped at 0.01 could therefore be accepted for a requested growth of 0.05, while the resulting model grew only at 0.01. The solver now keeps candidate bounds, checks incumbent feasibility and integrality, explicitly optimises the requested biomass reaction, and independently verifies the requested growth. Oxygen-direction restrictions are consistent between selection and application; ordinary intercompartment oxygen transport remains allowed.
- **Verification/reporting:** undefined-versus-finite discrepancies now fail verification, the verifier exits unsuccessfully on errors or no cards, and a request for an unconfigured curated model cannot silently load a draft under a curated label.

## ATP synthase: promising observation, narrower claim

The archived 500-model table contains **52 models with ATPS identifiers and 448 without them**, with no recorded download failures. The handover reversed the meaning of 52. The smaller table has 4 present and 29 absent. An equation/name screen of the four local models supports the need to investigate the three apparent omissions, but also flags an ATP-driven sodium exporter as an equation candidate: direction and physiological use must be checked separately.

The defensible current claim is: **an identifier screen found no ATPS-prefixed reaction in 448/500 sampled models from the historical EMBL collection.** That collection's latest upstream commit is [260d0f1, dated 5 March 2019](https://github.com/cdanielmachado/embl_gems/commit/260d0f133802adbc151de4d7b14fb34722e1d9f4). This is not a measurement of current CarveMe performance, absence in the organisms, or proof of the proposed homology/carving mechanism. Current [CarveMe universe curation](https://github.com/cdanielmachado/carveme/blob/5b25ce504de41060fb5dbaf4c87e13627c650632/carveme/universe/curation.py) explicitly retains selected ATP synthase representations.

The survey now parses actual SBML reactions, distinguishes ID/name/equation evidence, pins an immutable upstream revision, records sampling and file hashes, excludes failed downloads from absence estimates, and prevents accidental overwrites. The old 500 models were **not** downloaded and reanalysed during this audit; their archived aggregate was independently recounted. The old archive lacks hashes needed to prove its exact downloaded bytes.

The proposed ATP-synthase stoichiometry and all-subunit gene rule remain organism-specific hypotheses. Likewise, biomass menaquinol removal should be labelled provisional: missing classical `men` annotations alone does not exclude the experimentally established [alternative futalosine pathway](https://pubmed.ncbi.nlm.nih.gov/18801996/). Formula balance checks are useful, but largely zero charge annotations in the source SBML prevent meaningful charge certification without repairing that metadata.

## Whole-body IEM results need a fresh run

The secondary whole-body protocol has errors beyond a reporting correction:

1. Time-limited healthy/disease solves could be stored as 0/0 and classified as unchanged; this occurs in saved CPS1 biomarker results.
2. Demand sinks added for one disease could remain active for subsequent diseases, making results depend on run/shard order. The source Toolbox starts each IEM from the reference model.
3. The Python port omitted the source's ±100,000 bound on summed IEM flux. A saved healthy optimum above 115,000 demonstrates a changed feasible set.
4. The parser included commented biomarker rows and IEM calls inside literal `if 0` blocks, so the old 63-entry protocol was not the set of active source calls.

The corrected parser yields **57 active IEM calls and 252 biomarker entries**, compared with 63 and 279 previously. Six disabled calls account for 19 extra biomarker entries and commented rows for another eight. The saved shards contain 45 IEM records and 186 biomarkers, including three timeout pairs; only 41 records correspond to active corrected calls.

The fixes preserve failed outcomes as unknown, isolate each disease's temporary model state, restore the summed-flux cap, and produce a separately versioned parsed protocol with source provenance. The runner isolates corrected outputs from legacy results. Historical IEM files remain available for diagnosis, not validation claims. The old linked biomarker ground truth must also be regenerated and reconciled; it should not silently be relabelled as the corrected protocol.

A full corrected Harvey/Harvetta run and MATLAB differential comparison have **not** been performed. Harvey 1.03d is not included in this checkout, and the physiological/diet constraint port is unfinished. Comparing the old partial results with a published 85% figure would not be a like-for-like evaluation.

## Validation of the amendments

The complete regression suite passes: **69 tests plus six subtests**. Tests cover missing-data scoring, paired alignment, corrupt-card detection, reaction capacities, energy cycles, genuine solver failures, source parsing, IEM state isolation, interrupted-run retries and shared HiGHS scheduler use. The combined run exposed a scheduler conflict missed by isolated suites; default IEM solves now reuse the scheduler, and explicit incompatible configurations raise an actionable error. The exact environment and test output are saved with the audit. No full Harvey/Harvetta rerun is included.

## Next research priorities

1. Freeze a benchmark specification and the existing development set. Choose genuinely uninspected organisms/conditions for independent evaluation before further outcome-driven curation. Record exposure and freeze points; split at biologically meaningful groups, not random gene-condition cells.
2. Test a defined curation procedure against unmodified, gap-filled, current-tool and expert-curated baselines. Keep coverage, prediction accuracy and independent biochemical evidence separate. Include ablations, media/uptake and threshold sensitivity, and uncertainty. Fitness in a pooled mutant population should be described as fitness, not automatically as lethal single-gene deletion; the [Bernstein et al. 2023 protocol](https://doi.org/10.15252/msb.202311566) is a reference baseline.
3. Resolve the ATP-synthase observation mechanistically on a small pinned set: historical reconstruction scores, operon/GPR evidence, current reconstructions, flux direction and ATP yield. This is a concrete route to a useful result even if the broad novelty claim does not survive.
4. Rebuild and differentially test the IEM protocol before scaling it or interpreting disease accuracy. Keep this workstream separate until physiological assumptions match the comparator.

The large roadmap remains a collection of proposals, not a verified assessment of every field-wide or clinical claim it contains. This audit focused on the executable benchmark, headline survey and IEM port. It does not approve every gene rule or establish novelty against the entire literature.
