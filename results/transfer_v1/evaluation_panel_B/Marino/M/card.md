# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Marino

Created: 2026-10-03T13:29:21Z

## Model

- model_id: Marino_transfer_v1_base
- file: models/transfer_v1/Marino_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/Marino.xml.gz)
- version_note: 
- n_reactions: 1649
- n_metabolites: 1164
- n_genes: 959
- sha256: 7e95c275de7f729761dc874cdec7e70edc099dbf404e7f2b68c5bd373e67fd06

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Marino'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 3650
- n_experiments: 255

## Protocol

- study: transfer_v1
- arm: M
- arm_definition: {"transforms": ["universal_reaction_patches", "normalize_conjunction_rules", "universal_model_additions", "add_atp_synthase_if_annotated", "remove_menaquinol_if_no_pathway", "decisions:data/studies/transfer_v1/decisions/{org}.json"]}
- applied: [{"transform": "universal_reaction_patches", "records": [{"reaction": "CBMKr", "bounds_before": [-1000.0, 1000.0], "bounds_after": [-1000.0, 0.0]}, {"reaction": "CBPS", "gpr_before": "WP_014578297_1 or WP_014578298_1 or (WP_014578297_1 and WP_014578298_1)", "gpr_after": "WP_014578297_1 and WP_014578298_1"}]}, {"transform": "normalize_conjunction_rules", "records": [{"reaction": "ACOATA", "gpr_befo …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 1, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/studies/transfer_v1_replication/media_bigg_panel_B.tsv
- carbon_source_mapping: data/studies/transfer_v1_replication/carbon_sources_bigg_panel_B.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept_panel/Marino_genpept_map.tsv"}
- role: evaluation_panel_B
- base_model: {"org": "Marino", "role": "evaluation_panel_B", "draft": "models/embl_pinned/Marino.xml.gz", "draft_sha256": "13c5e03e67934b36dda53d1dd6ec44762155541c49e9c8424b7f7de7904ee617", "universe": "external/carveme/universe_bacteria.xml.gz", "universe_sha256": "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19", "reference": {"medium": "DinoMM_noCarbon_HighNutrient", "carbon_exchange": "EX_ …

## Leakage

- ground_truth_used_in_model_curation: no for the correction rules: frozen before this organism's fitness tables were downloaded (see the study freeze manifests); the base model's gap-fill uses observed wild-type growth on one reference condition
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: see docs/studies/transfer-method-v1.md and the panel selection record
- notes: ['Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 29, "n_conditions_wt_grows": 16, "conditions_with_absent_exchange": 15}
- gene_level_conditions_where_wt_grows: {"n_genes": 615, "n_conditions": 16, "n_gene_condition_pairs": 9840, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.3114028723436781, "ci95": [0.1587039627494976, 0.4730213186214669]}, "aucpr_standard": {"point": 0.21950664628606706, "ci95": [0.11799862757186601, 0.34164333220440685]}, "auroc_standard": {"point": 0.7315541367285432, "ci95": [0.6603125382203692, …
- gene_map: {"model_genes": 959, "mapped": 876, "matched_by_version": 958, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 82, "browser_genes_hit": 876, "mapped_with_fitness_data": 615}
- counts: {"model_genes": 959, "model_genes_mapped": 876, "genes_with_fitness": 615, "genes_after_adjustment": 615, "conditions_total": 30, "conditions_mapped": 29, "conditions_wt_grows": 16, "medium_completion_exchanges_added": 7}
- timings_s: {"rich_medium_essentials_s": 1.6689300537109375e-06, "knockout_simulation_s": 29.825432538986206, "total_s": 30.234694242477417}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
