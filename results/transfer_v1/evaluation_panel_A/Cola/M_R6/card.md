# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Cola

Created: 2026-10-03T10:15:51Z

## Model

- model_id: Cola_transfer_v1_base
- file: models/transfer_v1/Cola_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/Cola.xml.gz)
- version_note: 
- n_reactions: 1513
- n_metabolites: 1100
- n_genes: 730
- sha256: 2940de8c6e0497e446bf1d220ad34ee2351831e0de04b7ab466c8160b76a7af4

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Cola'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 3954
- n_experiments: 202

## Protocol

- study: transfer_v1
- arm: M_R6
- arm_definition: {"transforms": ["universal_reaction_patches", "normalize_conjunction_rules", "universal_model_additions", "add_atp_synthase_if_annotated", "remove_menaquinol_if_no_pathway", "decisions:data/studies/transfer_v1/decisions/{org}.json#R6"]}
- applied: [{"transform": "universal_reaction_patches", "records": [{"reaction": "CBPS", "gpr_before": "WP_015264089_1 or WP_015264609_1 or (WP_015264089_1 and WP_015264609_1)", "gpr_after": "WP_015264089_1 and WP_015264609_1"}]}, {"transform": "normalize_conjunction_rules", "records": [{"reaction": "ACOATA", "gpr_before": "WP_015263960_1 or WP_015264053_1 or WP_015266120_1 or WP_015267642_1 or (WP_010855324 …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/studies/transfer_v1/media_bigg.tsv
- carbon_source_mapping: data/studies/transfer_v1/carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept_panel/Cola_genpept_map.tsv"}
- role: evaluation_panel_A
- base_model: {"org": "Cola", "role": "evaluation_panel_A", "draft": "models/embl_pinned/Cola.xml.gz", "draft_sha256": "482301c55627b4186e9c305074ad03622b7abd9f7abb4af69b90c77b2ae6ee7f", "universe": "external/carveme/universe_bacteria.xml.gz", "universe_sha256": "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19", "reference": {"medium": "DinoMM_noCarbon_HighNutrient", "carbon_exchange": "EX_glc_ …

## Leakage

- ground_truth_used_in_model_curation: no for the correction rules: frozen before this organism's fitness tables were downloaded (see the study freeze manifests); the base model's gap-fill uses observed wild-type growth on one reference condition
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: see docs/studies/transfer-method-v1.md and the panel selection record
- notes: ['Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 14, "n_conditions_wt_grows": 7, "conditions_with_absent_exchange": 5}
- gene_level_conditions_where_wt_grows: {"n_genes": 500, "n_conditions": 7, "n_gene_condition_pairs": 3500, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.3735627729358546, "ci95": [0.22734457973355715, 0.5284680572908766]}, "aucpr_standard": {"point": 0.4181191870876605, "ci95": [0.29984606297832783, 0.5413835230175903]}, "auroc_standard": {"point": 0.6970620619398009, "ci95": [0.6329787271299243, 0 …
- gene_map: {"model_genes": 730, "mapped": 718, "matched_by_version": 720, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 11, "browser_genes_hit": 718, "mapped_with_fitness_data": 500}
- counts: {"model_genes": 730, "model_genes_mapped": 718, "genes_with_fitness": 500, "genes_after_adjustment": 500, "conditions_total": 19, "conditions_mapped": 14, "conditions_wt_grows": 7, "medium_completion_exchanges_added": 8}
- timings_s: {"rich_medium_essentials_s": 9.5367431640625e-07, "knockout_simulation_s": 11.009409666061401, "total_s": 11.32816767692566}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
