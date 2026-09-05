# Benchmark card — yeast_deletion_viability_v0 (yeast-GEM essentialGenes protocol)

Created: 2026-09-04T15:31:25Z

## Model

- model_id: yeast-GEM
- file: models/yeast-GEM.xml
- source: github.com/SysBioChalmers/yeast-GEM (main branch, model/yeast-GEM.xml)
- version_note: 
- n_reactions: 4105
- n_metabolites: 2748
- n_genes: 1143
- sha256: 30842b15eefb0ef7e36cbdea86a9efddfacf69a871c8b054165faa9af9f6c8eb

## Dataset provenance

- dataset: Stanford yeast deletion project inviable ORFs (14 Aug 2011 snapshot) and SGD verified ORFs (27 Aug 2013)
- copy_obtained_from: github.com/SysBioChalmers/yeast-GEM data/essentialGenes (yeast-GEM version 9.1.1)
- copy_commit: 2d594ae1c4a2d550ccef120d96a58c7bbf586255
- n_inviable_unique: 1122
- n_verified: 5061
- caveat: Screened in complex media supplemented with auxotrophic markers; the yeast-GEM test compares against a synthetic complete medium, so the reference is imperfect but stable (yeast-GEM README).

## Protocol

- ko_tol: 1e-06
- constrained_uptake: -0.5
- glucose_uptake: -20.0
- unconstrained_uptake: -1000.0
- restrict_to_verified: true
- processes: 1
- solver: glpk

## Leakage

- ground_truth_used_in_model_curation: yes — yeast-GEM has used this exact essential-gene test as a curation regression check since yeast7/8; treat as in-sample
- ground_truth_public_since: 2002 (Giaever et al.) / Stanford deletion project downloads
- frontier_model_training_exposure: almost certainly in training corpora
- held_out_recommendation: use condition-specific deletion phenotypes not used by the yeast-GEM test suite (e.g. Nichols-style chemical genomics for yeast, or newer Tn-seq/CRISPRi screens)
- notes: []

## Results

- mcc: {"point": 0.5322864939389729, "ci95": [0.450798412783923, 0.603826454462355]}
- accuracy: {"point": 0.9015356820234869, "ci95": [0.8825654923215899, 0.9177958446251129]}
- balanced_accuracy: {"point": 0.6964911233182072, "ci95": [0.6555835338534514, 0.7363957250513791]}
- aucpr_inviable_as_positive: {"point": 0.48581510670620215, "ci95": [0.40581237363721057, 0.5666235769338451]}
- counts: {"model_genes": 1143, "evaluated_genes": 1107, "exp_inviable": 159, "exp_viable": 948, "TP": 933, "TN": 65, "FP": 94, "FN": 15}
- wt_growth: 2.813938388921867
- missing_exchanges: []
- reference_yeastGEM_9.1.1_testResults: {"TP": 933, "TN": 65, "FP": 94, "FN": 15, "source": "external/yeast-GEM/data/testResults/essentialGenes.tsv"}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
