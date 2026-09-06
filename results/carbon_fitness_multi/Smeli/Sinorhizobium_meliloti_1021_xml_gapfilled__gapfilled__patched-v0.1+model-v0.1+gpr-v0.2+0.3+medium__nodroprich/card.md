# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Smeli

Created: 2026-09-06T07:45:13Z

## Model

- model_id: Sinorhizobium_meliloti_1021_xml_gapfilled
- file: models/gapfilled/Smeli.xml.gz
- source: EMBL GEMs (CarveMe draft, Machado et al. 2018) minimally gap-filled by gembench.gapfill for growth on RCH2_defined_noCarbon + EX_glc__D_e: added ['DXPS']
- version_note: 
- n_reactions: 1998
- n_metabolites: 1376
- n_genes: 1202
- sha256: aed3446b54ecac8a7a48e22c94e75fc57d32c35666a6431ae41678183cf85dfa

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Smeli'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 5133
- n_experiments: 90
- conditions_unmapped: 

## Protocol

- variant: gapfilled
- patches: {"file": "data/reference/universe_patches_v0.1.json", "model_file": "data/reference/model_patches_v0.1.json", "gpr_file": "data/reference/gpr_patches_v0.2.json,data/reference/gpr_patches_v0.3.json", "applied": [{"reaction": "CBMKr", "bounds_before": [-1000.0, 1000.0], "bounds_after": [-1000.0, 0.0]}, {"reaction": "CBPS", "gpr_before": "NP_385682_1 or NP_386426_1 or (NP_385682_1 and NP_386426_1)",  …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept/Smeli_genpept_map.tsv"}
- condition_selection: expGroup == 'carbon source'; condition_2 empty or DMSO; media with a BiGG mapping; replicates averaged per condition x medium
- scoring: gene-level metrics pooled over conditions where the wild-type model grows (>= growth_threshold); condition-level recall = fraction of experimentally growing carbon sources on which the wild-type model grows

## Leakage

- ground_truth_used_in_model_curation: no for EMBL draft models (automated reconstruction from genome annotation); for iML1515: partly (E. coli curation used phenotype data); for iJN1463: partly (Nogales et al. 2020 validated against growth phenotypes and gene essentiality data)
- ground_truth_public_since: Fitness Browser releases 2015-2018 (Price et al. 2018)
- frontier_model_training_exposure: Fitness Browser tables are public and partly in training corpora; the mapping tables here are new
- held_out_recommendation: unpublished RB-TnSeq experiments, or organisms added to the Browser after the model's training cut-off
- notes: ['Draft models are untouched by any phenotype data, so this is a true prospective test of automated reconstruction.']

## Results

- condition_level: {"n_conditions_mapped": 33, "n_conditions_wt_grows": 21, "wt_growth_recall": 0.6363636363636364, "conditions_with_absent_exchange": 6}
- gene_level_conditions_where_wt_grows: {"n_genes": 984, "n_conditions": 21, "n_gene_condition_pairs": 20664, "aucpr_bernstein": {"point": 0.6567302733210937, "ci95": [0.5350793195598229, 0.7888479474935663]}, "aucpr_standard": {"point": 0.583091203343994, "ci95": [0.4900090318107861, 0.676024153849841]}, "auroc_standard": {"point": 0.7762508820857809, "ci95": [0.7272495728338032, 0.8242887023166676]}, "mcc": {"point": 0.676422445488014 …
- gene_level_all_mapped_conditions: {"n_genes": 984, "n_conditions": 33, "n_gene_condition_pairs": 32472, "aucpr_bernstein": {"point": 0.5037833346896707, "ci95": [0.47030597708369976, 0.5374614369547243]}, "aucpr_standard": {"point": 0.27566796830667056, "ci95": [0.23033986606646656, 0.3187801826899798]}, "auroc_standard": {"point": 0.6864944518343227, "ci95": [0.6526330005567008, 0.7198042628748831]}, "mcc": {"point": 0.1994262032 …
- gene_map: {"model_genes": 1202, "mapped": 1201, "matched_by_version": 1195, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 0, "browser_genes_hit": 1201, "mapped_with_fitness_data": 984}
- counts: {"model_genes": 1202, "model_genes_mapped": 1201, "genes_with_fitness": 984, "genes_after_adjustment": 984, "conditions_total": 33, "conditions_mapped": 33, "conditions_wt_grows": 21, "medium_completion_exchanges_added": 6}
- timings_s: {"rich_medium_essentials_s": 7.152557373046875e-07, "knockout_simulation_s": 46.2934205532074, "total_s": 46.69479823112488}
- dropped_rich_medium_essentials: 0

## Warnings

- 0 of 33 conditions have no BiGG mapping
- medium 'RCH2_defined_noCarbon' components absent from the model: ['EX_sel_e', 'EX_slnt_e', 'EX_tungs_e', 'EX_lipoate_e', 'EX_pnto__R_e', 'EX_fol_e', 'EX_cbl1_e', 'EX_adocbl_e']
- medium completion added exchange+uptake for: ['EX_4abz_e', 'EX_btn_e', 'EX_na1_e', 'EX_nac_e', 'EX_ni2_e', 'EX_ribflv_e']

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
