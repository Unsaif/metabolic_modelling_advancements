# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Putida

Created: 2026-09-06T09:19:49Z

## Model

- model_id: iJN1463
- file: models/bigg/iJN1463.xml
- source: BiGG iJN1463 (Nogales et al. 2020), downloaded from bigg.ucsd.edu/static/models/iJN1463.xml on 5 Sept 2026
- version_note: 
- n_reactions: 2927
- n_metabolites: 2153
- n_genes: 1462
- sha256: d573833328ffae0dfa752a1fa3262ed939ed5862288beab287fca30d0fefb4a1

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Putida'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 4778
- n_experiments: 314
- conditions_unmapped: 1,3-Butandiol; 1,4-Butanediol; 1,5-Pentanediol; 1-Pentanol; 2-Piperidinone; 2-methyl-1-butanol; 3-methyl-3-butenol; 4-Hydroxyvalerate; Butyl stearate; Heptanoic acid; Levulinic Acid; Nonanoic acid; Tween 20; Valeric acid

## Protocol

- variant: curated
- patches: null
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": true, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": false, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "identity (BiGG gene ids are locus tags)"}
- condition_selection: expGroup == 'carbon source'; condition_2 empty or DMSO; media with a BiGG mapping; replicates averaged per condition x medium
- scoring: gene-level metrics pooled over conditions where the wild-type model grows (>= growth_threshold); condition-level recall = fraction of experimentally growing carbon sources on which the wild-type model grows

## Leakage

- ground_truth_used_in_model_curation: no for EMBL draft models (automated reconstruction from genome annotation); for iML1515: partly (E. coli curation used phenotype data); for iJN1463: partly (Nogales et al. 2020 validated against growth phenotypes and gene essentiality data)
- ground_truth_public_since: Fitness Browser releases 2015-2018 (Price et al. 2018)
- frontier_model_training_exposure: Fitness Browser tables are public and partly in training corpora; the mapping tables here are new
- held_out_recommendation: unpublished RB-TnSeq experiments, or organisms added to the Browser after the model's training cut-off
- notes: ['Draft models are untouched by any phenotype data, so this is a true prospective test of automated reconstruction.']

## Results

- condition_level: {"n_conditions_mapped": 43, "n_conditions_wt_grows": 36, "wt_growth_recall": 0.8372093023255814, "conditions_with_absent_exchange": 7}
- gene_level_conditions_where_wt_grows: {"n_genes": 1121, "n_conditions": 36, "n_gene_condition_pairs": 40356, "aucpr_bernstein": {"point": 0.47233015857784966, "ci95": [0.37152535193036834, 0.5580266208264155]}, "aucpr_standard": {"point": 0.3428738970115837, "ci95": [0.256437698613031, 0.42386762656500826]}, "auroc_standard": {"point": 0.8016704913965896, "ci95": [0.7481998373508131, 0.8446200541197451]}, "mcc": {"point": 0.4897286943 …
- gene_level_all_mapped_conditions: {"n_genes": 1121, "n_conditions": 43, "n_gene_condition_pairs": 48203, "aucpr_bernstein": {"point": 0.4083617621554833, "ci95": [0.3604269783477923, 0.45293244712058683]}, "aucpr_standard": {"point": 0.20963199482274977, "ci95": [0.1549064280564526, 0.2589917545452547]}, "auroc_standard": {"point": 0.7553018584536457, "ci95": [0.7121008220069435, 0.7915387073189045]}, "mcc": {"point": 0.2697495178 …
- gene_map: {"model_genes": 1462, "mapped": 1440, "mapped_with_fitness_data": 1148}
- counts: {"model_genes": 1462, "model_genes_mapped": 1440, "genes_with_fitness": 1148, "genes_after_adjustment": 1121, "conditions_total": 57, "conditions_mapped": 43, "conditions_wt_grows": 36, "medium_completion_exchanges_added": 0}
- timings_s: {"rich_medium_essentials_s": 8.103351831436157, "knockout_simulation_s": 338.03265929222107, "total_s": 347.13403058052063}
- dropped_rich_medium_essentials: 27

## Warnings

- 14 of 57 conditions have no BiGG mapping
- medium 'MOPS minimal media_noCarbon' components absent from the model: ['EX_slnt_e']
- medium 'RCH2_defined_noCarbon' components absent from the model: ['EX_slnt_e', 'EX_pydxn_e', 'EX_4abz_e', 'EX_lipoate_e', 'EX_ribflv_e', 'EX_thm_e', 'EX_fol_e']

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
