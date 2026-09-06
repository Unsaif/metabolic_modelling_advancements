# GEM-Bench v0 (Sprints 1–3, September 2026)

> Start with [the scientific audit](docs/reviews/2026-09-06-scientific-audit.md). Use `requirements-audit.txt` for the tested environment, a new `--output-dir` for carbon-fitness runs, and the versioned v0.2 IEM protocol. Commands and claims below describe the historical setup where they differ.

A small, explicit benchmark harness for genome-scale metabolic models, built for the
"Metabolic modelling improvements" project. Every dataset carries provenance and a
licence note; every protocol exposes its thresholds; every result is written with a
benchmark card (dataset, model, protocol, leakage, software versions).

## Layout

```
gembench/            package: media.py, datasets.py, metrics.py, cards.py, frog.py, protocols/,
                     fitness_browser.py (any-organism Fitness Browser loader + BiGG mapping), gene_mapping.py,
                     gapfill.py (minimal gap-fill on HiGHS), wbm.py / wbm_iem.py (whole-body models, IEM protocol)
scripts/             run_ecoli_fitness.py, run_yeast_essential.py, run_frog_cross_solver.py,
                     parse_iem_ground_truth.py, link_iem_ground_truth.py, run_wbm_solvers.py, run_wbm_iem.py,
                     gapfill_embl_models.py, run_carbon_fitness_generic.py, verify_carbon_fitness_multi.py
results/             benchmark cards (.card.json/.md), arrays (.npz), FROG fixtures, IEM tables
data/                downloaded reference data (HPO, Orphanet, Human-GEM gene table) — re-fetch, see below
external/            git clones used as data/model sources — re-fetch, see below
models/              symlinks/copies of the models used
```

## Explicit medium preparation for new studies

The opt-in `gembench.medium_preparation_v3` API returns a copy and a boundary
report. It closes validated environmental exchanges independently of inferred
compartment labels. Declare the model's compartment roles from its source:

```python
from gembench.media import Medium
from gembench.medium_preparation_v3 import CompartmentPolicy, prepare_medium

policy = CompartmentPolicy(external="C_e", cytoplasm="C_c",
                           external_aliases=("e",), cytoplasm_aliases=("c",))
medium = Medium("example", "Illustrative uptake capacities", {"EX_glc__D_e": -10.0})
prepared = prepare_medium(model, medium, policy=policy)
model_for_deletions = prepared.model
boundary_report = prepared.report
```

This example is not a biological medium recipe. Uptake values are flux bounds,
not concentrations. Missing components raise by default; `missing_policy="report"`
must be selected explicitly to allow them. Optional `completion_media` requires
an explicit protocol because adding transport changes the model. Prepare media
before gene deletion. Unsupported exchange orientations and undeclared naming
exceptions raise without changing the input. Internal demands/sinks are retained
and disclosed. The historical benchmark commands below retain their frozen
helpers; they have not been migrated. See the
[conformance sprint](docs/sprints/2026-09-06-medium-preparation.md), including
the curated-model numerical limitation.

## Reproducing the environment

```
pip install cobra highspy swiglpk python-libsbml osqp scikit-learn pandas
git clone --depth 1 https://github.com/dbernste/E_coli_GEM_validation external/E_coli_GEM_validation
git clone --depth 1 --filter=blob:none --sparse https://github.com/SysBioChalmers/yeast-GEM external/yeast-GEM
  (cd external/yeast-GEM && git sparse-checkout set model data code)
git clone --depth 1 --filter=blob:none --no-checkout https://github.com/opencobra/COBRA.models external/COBRA.models
  (cd external/COBRA.models && git checkout HEAD -- xml/ mat/Recon3DModel_301.mat mat/iMM904.mat ...)
git clone --depth 1 --filter=blob:none --no-checkout https://github.com/opencobra/cobratoolbox external/cobratoolbox
  (cd external/cobratoolbox && git checkout HEAD -- src/analysis/wholeBody/PSCMToolbox/runIEM_HH.m)
# reference data
curl -L -o data/hp.obo https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp.obo
curl -L -o data/phenotype.hpoa https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype.hpoa
curl -L -o data/humangem_genes.tsv https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/genes.tsv
# Orphanet product1/4/6 (CC-BY-4.0) via the Orphadata GitHub mirror (LFS files, served from media.githubusercontent.com)
```

## Running the multi-organism benchmark (Sprint 3)

```
# EMBL draft models (CC-BY 4.0) into models/embl/ from github.com/cdanielmachado/embl_gems (models/<x>/<genus>/<species>_<strain>.xml.gz)
# CarveMe universe into external/carveme/universe_bacteria.xml.gz from github.com/cdanielmachado/carveme (carveme/data/generated/)
python scripts/gapfill_embl_models.py                       # -> models/gapfilled/<org>.xml.gz + <org>_gapfill.json
python scripts/run_carbon_fitness_generic.py --variant gapfilled --orgs Btheta,Putida,MR1,Smeli
python scripts/run_carbon_fitness_generic.py --variant shipped   --orgs Btheta,Putida,MR1,Smeli,Keio
python scripts/verify_carbon_fitness_multi.py               # independent recomputation of every card
python scripts/run_carbon_fitness_generic.py --orgs Putida --variant curated     # iJN1463 (models/bigg/)
python scripts/compare_models_paired.py Putida <dir_A> <dir_B> draft iJN1463     # paired comparison on common genes/conditions
python scripts/find_bypasses.py                              # which reactions/gene rules protect persistently important genes
python scripts/run_carbon_fitness_generic.py --variant gapfilled --patch data/reference/universe_patches_v0.1.json
python scripts/evaluate_patch.py v0.1                        # verifier step of the fix-at-source loop
python scripts/run_carbon_fitness_generic.py --variant gapfilled --no-drop-rich --patch data/reference/universe_patches_v0.1.json --gpr-patch data/reference/gpr_patches_v0.2.json
python scripts/evaluate_patch.py "v0.1+gpr-v0.2__nodroprich"  # gene-rule patches, all genes with fitness data
python scripts/check_gene_rules.py Putida                    # annotation-based screen of OR gene rules (recall 36/51 on the hand-adjudicated set)
python scripts/run_wbm_iem.py Harvey_1_03d --iems GA2 OXOP --out-suffix _shard2   # IEM protocol in disjoint shards
```
Gene maps need data/genpept/<org>_genpept_map.tsv (NCBI efetch GenPept, parsed by tools/genpept_parse.js; see data/genpept/PROVENANCE.md).

## Running

```
python scripts/run_ecoli_fitness.py iJR904 iAF1260 iJO1366 iML1515 iML1515_corrected   # ~2 min per model
python scripts/run_yeast_essential.py                                                  # ~30 s
python scripts/run_frog_cross_solver.py --solvers glpk hybrid                          # hours (HiGHS via optlang is slow)
python scripts/parse_iem_ground_truth.py && python scripts/link_iem_ground_truth.py    # IEM ground truth v0.1
```

See docs/sprints/ for what each run found. Paths in the scripts assume this repository root as the working directory (external/, models/, data/, results/ alongside gembench/).

## Quinone development experiment and cofactor checks

The [fixed experiment definition](docs/studies/quinone-repair-v1.md) separates source-model chemistry, sequence-informed gene hypotheses, growth demand and benchmark coverage. Its [results](docs/sprints/2026-09-06-quinone-repair.md) remain a development evaluation. Use the pinned audit environment in `requirements-audit.txt` and fresh output paths:

```sh
.venv/bin/python scripts/run_quinone_repair.py --out results/quinone_repair_new
.venv/bin/python scripts/verify_quinone_repair.py --run-dir results/quinone_repair_new
.venv/bin/python scripts/verify_quinone_repair_solvers.py --study results/quinone_repair_new --out results/quinone_repair_new/solver_verification
```

The runner freezes its inputs before preparation, saves six candidate models and physical reports, and computes five benchmark arms. `--physical-only` omits numeric fitness-data loading. The source transfer is provisional; historical patch v0.4 is preserved. Candidate-model source terms and attribution are recorded separately from software in [MODEL_SOURCES.md](results/quinone_repair_2026_09_06/MODEL_SOURCES.md).

For other explicitly configured COBRApy models, `gembench.cofactor.pool_balance_certificate` checks a declared weighted pool and `probe_metabolite_production` maximizes separate outward demands under the existing constraints. See the [method and limitations](docs/studies/cofactor-audit-method.md); these helpers do not infer biological cofactor requirements.
