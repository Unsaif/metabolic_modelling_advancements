# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Keio

Created: 2026-09-06T09:42:54Z

## Model

- model_id: iML1515
- file: models/iML1515.xml
- source: BiGG iML1515 (Monk et al. 2017) via github.com/dbernste/E_coli_GEM_validation Models/
- version_note: 
- n_reactions: 2714
- n_metabolites: 1878
- n_genes: 1516
- sha256: 9c772d44ca43350e40dc7ee86c7aa148796856be1eea45e5406c6df8f7dcde28

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Keio'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 3789
- n_experiments: 168
- conditions_unmapped: 

## Protocol

- variant: shipped
- patches: {"file": "data/reference/universe_patches_v0.1.json", "model_file": "data/reference/model_patches_v0.2.json", "gpr_file": null, "applied": [{"reaction": "CBMKr", "bounds_before": [-1000.0, 1000.0], "bounds_after": [-1000.0, 0.0]}, {"reaction": "URIC", "bounds_before": [0.0, 1000.0], "bounds_after": [0.0, 0.0]}, {"reaction": "GCALDt", "added": "gcald_e <=> gcald_c", "gpr": "", "bounds": [-1000.0, 1 …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": true, "rich_medium_uptake": -1000.0, "knockout_genes": ["b0062", "b0063", "b3903", "b3904", "b0061", "b0344", "b3902"], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": false, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": nul …
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

- condition_level: {"n_conditions_mapped": 32, "n_conditions_wt_grows": 32, "wt_growth_recall": 1.0, "conditions_with_absent_exchange": 0}
- gene_level_conditions_where_wt_grows: {"n_genes": 1329, "n_conditions": 32, "n_gene_condition_pairs": 42528, "aucpr_bernstein": {"point": 0.6288231372296793, "ci95": [0.5166895301837996, 0.7168712114360852]}, "aucpr_standard": {"point": 0.49236375576195507, "ci95": [0.41014090290428723, 0.5648463734001795]}, "auroc_standard": {"point": 0.8033447407252235, "ci95": [0.7670054917220306, 0.8348084398501101]}, "mcc": {"point": 0.6108788863 …
- gene_level_all_mapped_conditions: {"n_genes": 1329, "n_conditions": 32, "n_gene_condition_pairs": 42528, "aucpr_bernstein": {"point": 0.6288231372296793, "ci95": [0.5166895301837996, 0.7168712114360852]}, "aucpr_standard": {"point": 0.49236375576195507, "ci95": [0.41014090290428723, 0.5648463734001795]}, "auroc_standard": {"point": 0.8033447407252235, "ci95": [0.7670054917220306, 0.8348084398501101]}, "mcc": {"point": 0.6108788863 …
- gene_map: {"model_genes": 1516, "mapped": 1515, "mapped_with_fitness_data": 1339}
- counts: {"model_genes": 1516, "model_genes_mapped": 1515, "genes_with_fitness": 1339, "genes_after_adjustment": 1329, "conditions_total": 32, "conditions_mapped": 32, "conditions_wt_grows": 32, "medium_completion_exchanges_added": 0}
- timings_s: {"rich_medium_essentials_s": 8.607654571533203, "knockout_simulation_s": 270.67939949035645, "total_s": 280.10890769958496}
- dropped_rich_medium_essentials: 10

## Warnings

- 0 of 32 conditions have no BiGG mapping
- medium 'MOPS Rich Defined media_noCarbon' components absent from the model: ['EX_4abz_e', 'EX_4hbz_e', 'EX_23dhb_e']
- medium 'RCH2_defined_noCarbon' components absent from the model: ['EX_4abz_e', 'EX_ribflv_e', 'EX_fol_e']

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
