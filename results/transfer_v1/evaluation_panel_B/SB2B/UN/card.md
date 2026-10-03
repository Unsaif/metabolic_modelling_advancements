# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism SB2B

Created: 2026-10-03T13:39:49Z

## Model

- model_id: SB2B_transfer_v1_base
- file: models/transfer_v1/SB2B_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/SB2B.xml.gz)
- version_note: 
- n_reactions: 2020
- n_metabolites: 1393
- n_genes: 880
- sha256: 951611b5768be36ad0ad10b4a643e5d07e6f2be0ad8c3a9aed797e669b7102cd

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'SB2B'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 3121
- n_experiments: 190

## Protocol

- study: transfer_v1
- arm: UN
- arm_definition: {"transforms": ["universal_reaction_patches", "normalize_conjunction_rules", "universal_model_additions", "add_atp_synthase_if_annotated"]}
- applied: [{"transform": "universal_reaction_patches", "records": [{"reaction": "CBPS", "gpr_before": "WP_011758957_1 or WP_011760612_1 or WP_011760613_1 or (WP_011760612_1 and WP_011760613_1)", "gpr_after": "WP_011760612_1 and WP_011760613_1"}]}, {"transform": "normalize_conjunction_rules", "records": [{"reaction": "ACOATA", "gpr_before": "WP_011759642_1 or WP_011760095_1 or WP_011760267_1 or WP_011761498_ …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 1, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/studies/transfer_v1_replication/media_bigg_panel_B.tsv
- carbon_source_mapping: data/studies/transfer_v1_replication/carbon_sources_bigg_panel_B.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept_panel/SB2B_genpept_map.tsv"}
- role: evaluation_panel_B
- base_model: {"org": "SB2B", "role": "evaluation_panel_B", "draft": "models/embl_pinned/SB2B.xml.gz", "draft_sha256": "6e43306f5050d2a7f66a7b296b8c355cff4fd011fd81f5d52a1098495afdf651", "universe": "external/carveme/universe_bacteria.xml.gz", "universe_sha256": "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19", "reference": {"medium": "ShewMM_noCarbon", "carbon_exchange": "EX_glc__D_e", "condi …

## Leakage

- ground_truth_used_in_model_curation: no for the correction rules: frozen before this organism's fitness tables were downloaded (see the study freeze manifests); the base model's gap-fill uses observed wild-type growth on one reference condition
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: see docs/studies/transfer-method-v1.md and the panel selection record
- notes: ['Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 20, "n_conditions_wt_grows": 15, "conditions_with_absent_exchange": 3}
- gene_level_conditions_where_wt_grows: {"n_genes": 648, "n_conditions": 15, "n_gene_condition_pairs": 9720, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.4463202730428332, "ci95": [0.3498800542581148, 0.5480256223001676]}, "aucpr_standard": {"point": 0.4895398478411719, "ci95": [0.405242528371406, 0.5662516909484572]}, "auroc_standard": {"point": 0.7553675251044514, "ci95": [0.7174702370474628, 0.7 …
- gene_map: {"model_genes": 880, "mapped": 863, "matched_by_version": 871, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 16, "browser_genes_hit": 863, "mapped_with_fitness_data": 648}
- counts: {"model_genes": 880, "model_genes_mapped": 863, "genes_with_fitness": 648, "genes_after_adjustment": 648, "conditions_total": 23, "conditions_mapped": 20, "conditions_wt_grows": 15, "medium_completion_exchanges_added": 8}
- timings_s: {"rich_medium_essentials_s": 2.384185791015625e-06, "knockout_simulation_s": 39.459418296813965, "total_s": 40.10746192932129}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
