# Benchmark card — ecoli_carbon_fitness_v0 (Bernstein et al. 2023 protocol)

Created: 2026-09-04T15:50:43Z

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

- variant: mannitol_excluded
- base_medium: {"medium": "M9_minimal_noCarbon_BiGG", "n_components": 23}
- carbon_uptake: -10.0
- growth_threshold: 0.001
- fitness_threshold: -2.0
- strain_adjustment: true
- drop_rich_medium_essentials: true
- rich_medium_uptake: -1000.0
- carbon_sources_to_remove: ["mnl", "sucr"]
- keep_atpm: true
- processes: 1
- solver: glpk

## Leakage

- ground_truth_used_in_model_curation: yes — this model was corrected against this dataset (in-sample)
- ground_truth_public_since: 2018 (Price et al., Nature; Fitness Browser)
- frontier_model_training_exposure: almost certainly present in pre-2026 training corpora (public TSVs, GitHub copies)
- held_out_recommendation: for any claim about model-driven curation, hold out organisms/conditions absent from the Fitness Browser at training cut-off, or use newly generated phenotypes
- notes: ["Sucrose/'man' exclusion: notebook comment says mannitol, code removes D-mannose ('man'); variant flag records which was used."]

## Results

- aucpr_bernstein: {"point": 0.7633293623908308, "ci95": [0.6834870253426781, 0.8308202914667444]}
- aucpr_standard: {"point": 0.5745661258834901, "ci95": [0.495804898035123, 0.6545259446380847]}
- auroc_standard: {"point": 0.8208962586743344, "ci95": [0.7855082489466291, 0.8542358804390205]}
- mcc: {"point": 0.6658550183920241, "ci95": [0.59844690961827, 0.7271190442279881]}
- balanced_accuracy: {"point": 0.7928863902532819, "ci95": [0.7560706216155036, 0.8283856783870077]}
- accuracy: {"point": 0.9493212669683257, "ci95": [0.9398461538461538, 0.9575301659125188]}
- confusion_growth_vs_important: {"tp": 29642, "tn": 1828, "fp": 1212, "fn": 468}
- counts: {"model_genes": 1516, "genes_matched": 1339, "carbon_sources_in_data": 27, "genes_after_adjustment": 1326, "carbon_sources_after_adjustment": 25}
- dropped_strain_genes: ["b0061", "b0344", "b3902"]
- n_dropped_rich_medium_essentials: 10
- dropped_carbon_sources: ["mnl", "sucr"]
- missing_medium_components: []
- missing_carbon_exchanges: []
- per_carbon_source: {"ac": {"aucpr_bernstein": 0.6841594117648498, "mcc": 0.5875705139504301, "wt_growth": 0.21033012489208694, "n_pred_no_growth": 90, "n_exp_important": 145}, "acgam": {"aucpr_bernstein": 0.821905529396446, "mcc": 0.696231872757012, "wt_growth": 1.1324950216674115, "n_pred_no_growth": 92, "n_exp_important": 115}, "akg": {"aucpr_bernstein": 0.6989860315780024, "mcc": 0.5881771016267502, "wt_growth":  …
- timings_s: {"rich_medium_essentials_s": 8.098396301269531, "knockout_simulation_s": 195.957524061203, "total_s": 204.37143635749817}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
