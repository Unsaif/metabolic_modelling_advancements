# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Cup4G11

Created: 2026-10-03T13:21:44Z

## Model

- model_id: Cup4G11_transfer_v1_base
- file: models/transfer_v1/Cup4G11_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/Cup4G11.xml.gz)
- version_note: 
- n_reactions: 2218
- n_metabolites: 1493
- n_genes: 1305
- sha256: 37c47a0dc229b340eebd0a44fd6b5ac270b399905f7e182921968784c077da90

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Cup4G11'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 6384
- n_experiments: 132

## Protocol

- study: transfer_v1
- arm: M_R6
- arm_definition: {"transforms": ["universal_reaction_patches", "normalize_conjunction_rules", "universal_model_additions", "add_atp_synthase_if_annotated", "remove_menaquinol_if_no_pathway", "decisions:data/studies/transfer_v1/decisions/{org}.json#R6"]}
- applied: [{"transform": "universal_reaction_patches", "records": [{"reaction": "CBPS", "gpr_before": "WP_043347612_1 or WP_043347616_1 or (WP_043347612_1 and WP_043347616_1)", "gpr_after": "WP_043347612_1 and WP_043347616_1"}]}, {"transform": "normalize_conjunction_rules", "records": [{"reaction": "ACOATA", "gpr_before": "WP_043347879_1 or WP_043347892_1 or WP_043355695_1 or WP_043355703_1 or (WP_009522116 …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 1, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/studies/transfer_v1_replication/media_bigg_panel_B.tsv
- carbon_source_mapping: data/studies/transfer_v1_replication/carbon_sources_bigg_panel_B.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept_panel/Cup4G11_genpept_map.tsv"}
- role: evaluation_panel_B
- base_model: {"org": "Cup4G11", "role": "evaluation_panel_B", "draft": "models/embl_pinned/Cup4G11.xml.gz", "draft_sha256": "bc51082574bcc225b06adc2e7dc2f2a142678a244b68607725c0fb86640502d6", "universe": "external/carveme/universe_bacteria.xml.gz", "universe_sha256": "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19", "reference": {"medium": "RCH2_defined_noCarbon", "carbon_exchange": "EX_fru_e …

## Leakage

- ground_truth_used_in_model_curation: no for the correction rules: frozen before this organism's fitness tables were downloaded (see the study freeze manifests); the base model's gap-fill uses observed wild-type growth on one reference condition
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: see docs/studies/transfer-method-v1.md and the panel selection record
- notes: ['Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 36, "n_conditions_wt_grows": 16, "conditions_with_absent_exchange": 6}
- gene_level_conditions_where_wt_grows: {"n_genes": 961, "n_conditions": 16, "n_gene_condition_pairs": 15376, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.29328111560959075, "ci95": [0.1319393471702188, 0.5003080122660692]}, "aucpr_standard": {"point": 0.24673791280875823, "ci95": [0.12240182047374117, 0.39728459238225694]}, "auroc_standard": {"point": 0.6542986908521148, "ci95": [0.587797019297676 …
- gene_map: {"model_genes": 1305, "mapped": 1219, "matched_by_version": 1293, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 85, "browser_genes_hit": 1219, "mapped_with_fitness_data": 961}
- counts: {"model_genes": 1305, "model_genes_mapped": 1219, "genes_with_fitness": 961, "genes_after_adjustment": 961, "conditions_total": 38, "conditions_mapped": 36, "conditions_wt_grows": 16, "medium_completion_exchanges_added": 5}
- timings_s: {"rich_medium_essentials_s": 1.430511474609375e-06, "knockout_simulation_s": 66.33268451690674, "total_s": 66.80836987495422}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
