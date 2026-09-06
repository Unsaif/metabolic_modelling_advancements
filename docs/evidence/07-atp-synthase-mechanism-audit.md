# Evidence brief 07 — Historical ATP synthase omissions: mechanism audit

Date: 6 September 2026. Scope: source history, reaction metadata and genome annotations for the four existing development organisms (Btheta, Putida, MR1, Smeli). No phenotype outcomes for additional organisms were opened. No model or patch assumptions were changed.

The most specific supported explanation is a **historical missing-subunit scoring failure that upstream explicitly changed in January 2018**. The affected EMBL models were already uploaded in November 2017. This establishes a plausible mechanism and its chronology; it does **not** establish which alignment was missing in each organism or prove that this mechanism explains all 448 identifier-negative models. It also changes the novelty claim: the evidence concerns an old model collection, not a newly demonstrated defect in current CarveMe.

## 1. The models predate the upstream change

The original EMBL model upload is commit [`e14497c154dfb6e9e003729e58f9dffdf8f65579`](https://github.com/cdanielmachado/embl_gems/commit/e14497c154dfb6e9e003729e58f9dffdf8f65579), dated 11 November 2017. The later snapshot [`260d0f133802adbc151de4d7b14fb34722e1d9f4`](https://github.com/cdanielmachado/embl_gems/commit/260d0f133802adbc151de4d7b14fb34722e1d9f4), dated 5 March 2019, adds formula/charge metadata. The [collection description](https://github.com/cdanielmachado/embl_gems/blob/260d0f133802adbc151de4d7b14fb34722e1d9f4/README.md) identifies RefSeq release 84.

We downloaded the four development models at **both immutable commits** and compared all reaction identifiers and the XML reactant/product lists. Both matched exactly for every model. The ATP identifier result therefore predates the 2019 metadata update:

| Organism | Reactions at each snapshot | ATPS identifier, 2017 and 2019 |
|---|---:|---|
| B. thetaiotaomicron VPI-5482 | 1,487 | None |
| P. putida KT2440 | 1,895 | `R_ATPS4rpp` |
| S. oneidensis MR-1 | 1,972 | None |
| S. meliloti 1021 | 1,992 | None |

Full download URLs and SHA-256 hashes of all eight compressed files are in [`results/atp_mechanism/model_snapshot_comparison.json`](../../results/atp_mechanism/model_snapshot_comparison.json). This comparison does not assert equality of every annotation or flux bound. The previously archived 500-model identifier screen remains 52 present and 448 absent; it was not rerun for this mechanism audit.

## 2. A specific historical algorithmic failure is demonstrated

The latest source commit dated before the original model upload is [`86e1ef8a39f3f1af32309162894e4f8d999bf955`](https://github.com/cdanielmachado/carveme/commit/86e1ef8a39f3f1af32309162894e4f8d999bf955), 9 November 2017. It is a temporally appropriate comparison source, **not a recovered build identifier for the archived models**.

In its [scoring implementation](https://github.com/cdanielmachado/carveme/blob/86e1ef8a39f3f1af32309162894e4f8d999bf955/carveme/reconstruction/scoring.py):

- A reference complex loses its GPR if any subunit has no alignment hit.
- Its score is the minimum subunit score with missing values propagated.
- Alternative reference complexes are combined by maximum score; an intact alternative can therefore rescue the reaction.
- If every reference complex loses its evidence, the reaction is absent from the reaction-score table.

The corresponding [carving implementation](https://github.com/cdanielmachado/carveme/blob/86e1ef8a39f3f1af32309162894e4f8d999bf955/carveme/reconstruction/carving.py) gives ordinary unscored internal reactions a default penalty of −1. It optimizes total annotation agreement under network/growth constraints. It does not require ATP synthase specifically. **Loss of score creates pressure to remove a reaction; removal is not inevitable**, because a penalized reaction may still be needed for feasible growth or other rewarded reactions.

On 31 January 2018, upstream committed [`5038ece6e445862eb5e73ddf1bb9258e67a780ff`](https://github.com/cdanielmachado/carveme/commit/5038ece6e445862eb5e73ddf1bb9258e67a780ff), titled “don't ignore complexes with missing subunits.” It replaces the missing-sensitive minimum with a mean in which missing scores contribute zero, and constructs GPRs from the remaining matched genes. The currently inspected [v1.6.6 scoring source](https://github.com/cdanielmachado/carveme/blob/5b25ce504de41060fb5dbaf4c87e13627c650632/carveme/reconstruction/scoring.py) retains those mean/maximum aggregation choices.

There is a methods/version discrepancy worth preserving: the [2018 primary paper, Reaction scoring section](https://doi.org/10.1093/nar/gky537), describes minimum subunit scores and a sum over isozymes. Inspected November 2017 code uses minimum/maximum; January 2018 code uses mean/maximum. The paper alone cannot identify which implementation generated a particular downloadable model.

## 3. Controlled experiment on the unmodified source

[`scripts/diagnose_atp_subunit_scoring.py`](../../scripts/diagnose_atp_subunit_scoring.py) executes the two pinned scoring modules on synthetic alignment/GPR inputs: two alternative eight-subunit reference complexes, plus an unrelated control reaction. It invokes the actual `reaction_scoring` functions, not a local approximation of them. The same inputs are passed to both versions. It does not run an optimizer or use experimental fitness.

| Synthetic evidence in each reference complex | November 2017 raw ATP score | January 2018 raw ATP score |
|---|---:|---:|
| All eight subunits score 80 | 80 | 80 |
| Seven score 80; the eighth has no hit | No scored reaction | 70 |
| Seven score 80; the eighth scores 1 | 1 | 70.125 |
| No subunit hits | No scored reaction | No scored reaction |
| First template missing one; alternative template complete | 80 | 80 |

All ten pre-specified score checks passed. The source hashes, resulting GPRs and pandas version are archived in [`scoring_counterexample.json`](../../results/atp_mechanism/scoring_counterexample.json). Reproduce with a local CarveMe Git clone containing both source snapshots:

```sh
python scripts/diagnose_atp_subunit_scoring.py --carveme-repo /path/to/carveme
```

The script uses no network and refuses to overwrite a specified output. A partial clone must already contain the required source/data blobs. Its conclusions concern **scoring behavior**, not model retention or physiological ATP production.

The fixed code also illustrates a separate limitation: when the eighth subunit is missing, it writes a seven-gene GPR. Better reaction retention does not by itself establish a biologically complete complex or correct essential-subunit predictions. This audit does not adopt that partial GPR as a biological correction.

## 4. Reference coverage supports a testable hypothesis

The original [GPR database](https://github.com/cdanielmachado/carveme/blob/86e1ef8a39f3f1af32309162894e4f8d999bf955/carveme/data/generated/bigg_gprs.csv.gz), joined to its [model-organism metadata](https://github.com/cdanielmachado/carveme/blob/86e1ef8a39f3f1af32309162894e4f8d999bf955/carveme/data/input/bigg_models.csv), provides 1,011 ATP-synthase gene rows in 119 reference complexes from 60 models. These are highly concentrated: Escherichia 49, Shigella 7, and one each of Klebsiella, Salmonella, Yersinia and Pseudomonas. The Pseudomonas reference is **P. putida KT2440 itself** (`iJN746`). These counts are independently regenerated by the diagnostic.

This makes the proposed mechanism specific: relatively divergent or difficult-to-align subunits could eliminate every complete reference-complex match under the old rule, even when several catalytic subunits align well. The positive P. putida development model has same-strain template coverage; it is not an independent test of transfer to a distant organism. These observations are consistent with the hypothesis, not proof that distance caused the three omissions.

Current cached RefSeq feature tables identify the expected eight core F-type subunits in Btheta and MR1 and nine entries in Smeli because of its two B-subunit annotations. Putida has the eight core entries as well. Accessory subunit I and Btheta's V-type system should not be confused with the core F-type set. These are present-day annotation observations, not reconstructed 2017 alignment evidence. Their accessions and download provenance are under `data/ncbi_feature_tables/`.

## 5. What remains unproven, and the next small experiment

Still missing are the archived build command/version, exact 2017 protein inputs, DIAMOND version/options and hit tables. No claim is yet justified that a particular subunit, taxonomic distance, solver choice or energy pathway caused the observed omission. The original 500-model identifier count cannot be converted into a biological error rate. An existing upstream fix also prevents presenting the broad scoring mechanism as newly invented here.

**Next experiment, fixed in advance:** use only Btheta, MR1 and Smeli, with Putida as the positive development control. Recover and hash the best available release-84 protein FASTAs and the November 2017 BiGG protein database. First run annotation only, retaining all per-subunit hit scores, reference IDs, coverage and identities. No phenotype-guided choice of thresholds or GPRs is allowed.

For each organism, evaluate the identical hit table with both pinned scoring versions. Record whether every historical ATP reference complex has at least one missing hit, which subunits fail, whether ATP score disappears under the old rule, and whether it returns under the fix. Preserve complete-template controls. If the exact historical inputs cannot be recovered, label the run a contemporary-input mechanism test rather than a historical reconstruction.

Only if that first stage supports the hypothesis, perform a small paired reconstruction holding the universe, alignments, solver, bounds and other scores fixed while changing the ATP-synthase scoring treatment. Compare ATP retention and the MILP objective tradeoff before any fitness outcomes. If old scoring already provides a complete positive ATP score, investigate network consistency/optimization or input-version differences instead. A forced-retention comparison can measure an objective cost; it must remain a diagnostic and must not become an accepted biological patch merely because growth improves.
