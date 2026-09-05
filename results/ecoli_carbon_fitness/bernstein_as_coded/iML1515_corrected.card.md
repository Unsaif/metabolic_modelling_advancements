# Benchmark card — ecoli_carbon_fitness_v0 (Bernstein et al. 2023 protocol)

Created: 2026-09-04T15:38:15Z

## Model

- model_id: iML1515_corrected
- file: external/E_coli_GEM_validation/Analysis/iML1515_model_adjusted_all_corrections.xml
- source: Bernstein et al. 2023 iML1515 with all corrections (Analysis/iML1515_model_adjusted_all_corrections.xml)
- version_note: 
- n_reactions: 2714
- n_metabolites: 1877
- n_genes: 1516
- sha256: df10ea412e87327d08dda737934d8be64cfb0378e9faa3441b40f857f18e9707

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, E. coli BW25113 (orgId 'Keio')
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- copy_obtained_from: github.com/dbernste/E_coli_GEM_validation (MIT-licensed code; data redistributed there)
- copy_commit: 7e34a7df9bc16b613ac9776e2cee27133528e4c2
- selection: bernstein2023
- n_genes: 3789
- n_carbon_sources: 27
- n_experiments: 54
- units: log2 gene fitness (Fitness Browser convention); replicates averaged per carbon source
- licence_note: Fitness Browser data are publicly available; check fit.genomics.lbl.gov terms before redistribution.

## Protocol

- variant: bernstein_as_coded
- base_medium: {"medium": "M9_minimal_noCarbon_BiGG", "n_components": 23}
- carbon_uptake: -10.0
- growth_threshold: 0.001
- fitness_threshold: -2.0
- strain_adjustment: true
- drop_rich_medium_essentials: true
- rich_medium_uptake: -1000.0
- carbon_sources_to_remove: ["man", "sucr"]
- keep_atpm: true
- processes: 2
- solver: glpk

## Leakage

- ground_truth_used_in_model_curation: yes — this model was corrected against this dataset (in-sample)
- ground_truth_public_since: 2018 (Price et al., Nature; Fitness Browser)
- frontier_model_training_exposure: almost certainly present in pre-2026 training corpora (public TSVs, GitHub copies)
- held_out_recommendation: for any claim about model-driven curation, hold out organisms/conditions absent from the Fitness Browser at training cut-off, or use newly generated phenotypes
- notes: ["Sucrose/'man' exclusion: notebook comment says mannitol, code removes D-mannose ('man'); variant flag records which was used."]

## Results

- aucpr_bernstein: {"point": 0.7642940162026393, "ci95": [0.6853601561829478, 0.8309291000380492]}
- aucpr_standard: {"point": 0.5668779208878092, "ci95": [0.4859853814537711, 0.6462389998491415]}
- auroc_standard: {"point": 0.8186639065845487, "ci95": [0.7818875767351581, 0.8525744666833056]}
- mcc: {"point": 0.6642491670361856, "ci95": [0.596903229644515, 0.7257940574114938]}
- balanced_accuracy: {"point": 0.7917015817076392, "ci95": [0.7545253179587057, 0.8276013789401655]}
- accuracy: {"point": 0.9490196078431372, "ci95": [0.9394871794871795, 0.957317496229261]}
- confusion_growth_vs_important: {"tp": 29635, "tn": 1825, "fp": 1222, "fn": 468}
- counts: {"model_genes": 1516, "genes_matched": 1339, "carbon_sources_in_data": 27, "genes_after_adjustment": 1326, "carbon_sources_after_adjustment": 25}
- dropped_strain_genes: ["b0061", "b0344", "b3902"]
- n_dropped_rich_medium_essentials: 10
- dropped_carbon_sources: ["man", "sucr"]
- missing_medium_components: []
- missing_carbon_exchanges: []
- per_carbon_source: {"ac": {"aucpr_bernstein": 0.6841594117648498, "mcc": 0.5875705139504301, "wt_growth": 0.2103301248920697, "n_pred_no_growth": 90, "n_exp_important": 145}, "acgam": {"aucpr_bernstein": 0.821905529396446, "mcc": 0.696231872757012, "wt_growth": 1.132495021667263, "n_pred_no_growth": 92, "n_exp_important": 115}, "akg": {"aucpr_bernstein": 0.6989860315780024, "mcc": 0.5881771016267502, "wt_growth": 0. …
- timings_s: {"rich_medium_essentials_s": 4.684260129928589, "knockout_simulation_s": 105.56063747406006, "total_s": 110.49606323242188}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
