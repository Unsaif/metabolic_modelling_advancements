# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Smeli

Created: 2026-09-06T09:11:53Z

## Model

- model_id: Sinorhizobium_meliloti_1021_xml_gapfilled
- file: models/gapfilled/Smeli.xml.gz
- source: EMBL GEMs (CarveMe draft, Machado et al. 2018) minimally gap-filled by gembench.gapfill for growth on RCH2_defined_noCarbon + EX_glc__D_e: added ['DXPS']
- version_note: 
- n_reactions: 1993
- n_metabolites: 1373
- n_genes: 1196
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
- patches: {"file": "data/reference/universe_patches_v0.1.json", "model_file": null, "gpr_file": "data/reference/gpr_patches_v0.2.json", "applied": [{"reaction": "CBMKr", "bounds_before": [-1000.0, 1000.0], "bounds_after": [-1000.0, 0.0]}, {"reaction": "CBPS", "gpr_before": "NP_385682_1 or NP_386426_1 or (NP_385682_1 and NP_386426_1)", "gpr_after": "NP_385682_1 and NP_386426_1"}, {"reaction": "UMPK", "gpr_be …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": true, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": false, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
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
- gene_level_conditions_where_wt_grows: {"n_genes": 954, "n_conditions": 21, "n_gene_condition_pairs": 20034, "aucpr_bernstein": {"point": 0.5891940373660133, "ci95": [0.4217242408289087, 0.7250426871713316]}, "aucpr_standard": {"point": 0.5087386152911402, "ci95": [0.3899723588205477, 0.601243277047567]}, "auroc_standard": {"point": 0.7325862407816986, "ci95": [0.6727597791816604, 0.7811172009664493]}, "mcc": {"point": 0.62803810936659 …
- gene_level_all_mapped_conditions: {"n_genes": 954, "n_conditions": 33, "n_gene_condition_pairs": 31482, "aucpr_bernstein": {"point": 0.46694223614978725, "ci95": [0.435717367291273, 0.49472300239469097]}, "aucpr_standard": {"point": 0.21070146179123556, "ci95": [0.15620814384983628, 0.2567213630426936]}, "auroc_standard": {"point": 0.6549174482665888, "ci95": [0.6142770046696513, 0.6865000879447881]}, "mcc": {"point": 0.1573597012 …
- gene_map: {"model_genes": 1196, "mapped": 1195, "matched_by_version": 1195, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 0, "browser_genes_hit": 1195, "mapped_with_fitness_data": 979}
- counts: {"model_genes": 1196, "model_genes_mapped": 1195, "genes_with_fitness": 979, "genes_after_adjustment": 954, "conditions_total": 33, "conditions_mapped": 33, "conditions_wt_grows": 21, "medium_completion_exchanges_added": 0}
- timings_s: {"rich_medium_essentials_s": 5.266183614730835, "knockout_simulation_s": 88.56570768356323, "total_s": 94.5303041934967}
- dropped_rich_medium_essentials: 25

## Warnings

- 0 of 33 conditions have no BiGG mapping
- medium 'RCH2_defined_noCarbon' components absent from the model: ['EX_na1_e', 'EX_ni2_e', 'EX_sel_e', 'EX_slnt_e', 'EX_tungs_e', 'EX_4abz_e', 'EX_lipoate_e', 'EX_nac_e', 'EX_ribflv_e', 'EX_thm_e', 'EX_pnto__R_e', 'EX_btn_e', 'EX_fol_e', 'EX_cbl1_e', 'EX_adocbl_e']

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
