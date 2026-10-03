# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism PV4

Created: 2026-10-03T13:36:54Z

## Model

- model_id: PV4_transfer_v1_base
- file: models/transfer_v1/PV4_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/PV4.xml.gz)
- version_note: 
- n_reactions: 2016
- n_metabolites: 1394
- n_genes: 911
- sha256: ef0dbe3b2d6351cff66d1ac11e5459a86e364c5980febfbc8f1e9f8fa52ecc77

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'PV4'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 3009
- n_experiments: 160

## Protocol

- study: transfer_v1
- arm: M_R1R2
- arm_definition: {"transforms": ["universal_reaction_patches", "normalize_conjunction_rules", "universal_model_additions", "add_atp_synthase_if_annotated", "remove_menaquinol_if_no_pathway", "decisions:data/studies/transfer_v1/decisions/{org}.json#R1,R2,R2+R1"]}
- applied: [{"transform": "universal_reaction_patches", "records": [{"reaction": "CBPS", "gpr_before": "WP_011864817_1 or WP_011866635_1 or WP_011866636_1 or (WP_011866635_1 and WP_011866636_1)", "gpr_after": "WP_011866635_1 and WP_011866636_1"}]}, {"transform": "normalize_conjunction_rules", "records": [{"reaction": "ACOATA", "gpr_before": "WP_011865134_1 or WP_011865400_1 or WP_011866163_1 or WP_011866207_ …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 1, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/studies/transfer_v1_replication/media_bigg_panel_B.tsv
- carbon_source_mapping: data/studies/transfer_v1_replication/carbon_sources_bigg_panel_B.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept_panel/PV4_genpept_map.tsv"}
- role: evaluation_panel_B
- base_model: {"org": "PV4", "role": "evaluation_panel_B", "draft": "models/embl_pinned/PV4.xml.gz", "draft_sha256": "dfec4d30777c80b4215aa1958c29ef377e7df37d9671b1e2981b1c9e762a78b7", "universe": "external/carveme/universe_bacteria.xml.gz", "universe_sha256": "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19", "reference": {"medium": "ShewMM_noCarbon", "carbon_exchange": "EX_lac__L_e", "conditi …

## Leakage

- ground_truth_used_in_model_curation: no for the correction rules: frozen before this organism's fitness tables were downloaded (see the study freeze manifests); the base model's gap-fill uses observed wild-type growth on one reference condition
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: see docs/studies/transfer-method-v1.md and the panel selection record
- notes: ['Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 12, "n_conditions_wt_grows": 10, "conditions_with_absent_exchange": 2}
- gene_level_conditions_where_wt_grows: {"n_genes": 645, "n_conditions": 10, "n_gene_condition_pairs": 6450, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.4922591129301061, "ci95": [0.35707633809386985, 0.6195847084480428]}, "aucpr_standard": {"point": 0.4881998699313506, "ci95": [0.3986809614890483, 0.5799745343170389]}, "auroc_standard": {"point": 0.7375295343209312, "ci95": [0.6929039358946693, 0 …
- gene_map: {"model_genes": 911, "mapped": 897, "matched_by_version": 902, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 13, "browser_genes_hit": 897, "mapped_with_fitness_data": 645}
- counts: {"model_genes": 911, "model_genes_mapped": 897, "genes_with_fitness": 645, "genes_after_adjustment": 645, "conditions_total": 13, "conditions_mapped": 12, "conditions_wt_grows": 10, "medium_completion_exchanges_added": 8}
- timings_s: {"rich_medium_essentials_s": 1.430511474609375e-06, "knockout_simulation_s": 26.376359224319458, "total_s": 26.86467742919922}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
