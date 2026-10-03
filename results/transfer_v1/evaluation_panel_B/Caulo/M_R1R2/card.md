# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Caulo

Created: 2026-10-03T13:10:14Z

## Model

- model_id: Caulo_transfer_v1_base
- file: models/transfer_v1/Caulo_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/Caulo.xml.gz)
- version_note: 
- n_reactions: 1693
- n_metabolites: 1259
- n_genes: 750
- sha256: 78f1ceee721a644a8845ca3dc8974f89f3eba536835a7695a143d7777f69180f

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Caulo'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 3312
- n_experiments: 198

## Protocol

- study: transfer_v1
- arm: M_R1R2
- arm_definition: {"transforms": ["universal_reaction_patches", "normalize_conjunction_rules", "universal_model_additions", "add_atp_synthase_if_annotated", "remove_menaquinol_if_no_pathway", "decisions:data/studies/transfer_v1/decisions/{org}.json#R1,R2,R2+R1"]}
- applied: [{"transform": "universal_reaction_patches", "records": [{"reaction": "CBPS", "gpr_before": "YP_002518298_1 or YP_002518367_1 or (YP_002518298_1 and YP_002518367_1)", "gpr_after": "YP_002518298_1 and YP_002518367_1"}]}, {"transform": "normalize_conjunction_rules", "records": [{"reaction": "AACPS1", "gpr_before": "YP_002517180_1 or (YP_002517122_1 and YP_002517180_1)", "gpr_after": "YP_002517122_1  …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 1, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/studies/transfer_v1_replication/media_bigg_panel_B.tsv
- carbon_source_mapping: data/studies/transfer_v1_replication/carbon_sources_bigg_panel_B.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept_panel/Caulo_genpept_map.tsv"}
- role: evaluation_panel_B
- base_model: {"org": "Caulo", "role": "evaluation_panel_B", "draft": "models/embl_pinned/Caulo.xml.gz", "draft_sha256": "a10f3d4a2da6629ff1f37d8e344e618bf7724c3f9c1ad23e682e8fb5ca5acb6b", "universe": "external/carveme/universe_bacteria.xml.gz", "universe_sha256": "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19", "reference": {"medium": "M2_noCarbon", "carbon_exchange": "EX_glc__D_e", "conditi …

## Leakage

- ground_truth_used_in_model_curation: no for the correction rules: frozen before this organism's fitness tables were downloaded (see the study freeze manifests); the base model's gap-fill uses observed wild-type growth on one reference condition
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: see docs/studies/transfer-method-v1.md and the panel selection record
- notes: ['Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 13, "n_conditions_wt_grows": 5, "conditions_with_absent_exchange": 6}
- gene_level_conditions_where_wt_grows: {"n_genes": 547, "n_conditions": 5, "n_gene_condition_pairs": 2735, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.4210000945957048, "ci95": [0.3212736761429678, 0.5278579303519428]}, "aucpr_standard": {"point": 0.42036464683255237, "ci95": [0.3463509032925566, 0.504792344964385]}, "auroc_standard": {"point": 0.7389721282944883, "ci95": [0.696459271113091, 0.78 …
- gene_map: {"model_genes": 750, "mapped": 748, "matched_by_version": 741, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 1, "browser_genes_hit": 748, "mapped_with_fitness_data": 547}
- counts: {"model_genes": 750, "model_genes_mapped": 748, "genes_with_fitness": 547, "genes_after_adjustment": 547, "conditions_total": 14, "conditions_mapped": 13, "conditions_wt_grows": 5, "medium_completion_exchanges_added": 2}
- timings_s: {"rich_medium_essentials_s": 1.1920928955078125e-06, "knockout_simulation_s": 8.855463743209839, "total_s": 9.258094072341919}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
