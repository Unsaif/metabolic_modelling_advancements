# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Cola

Created: 2026-10-03T10:15:18Z

## Model

- model_id: Cola_transfer_v1_base
- file: models/transfer_v1/Cola_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/Cola.xml.gz)
- version_note: 
- n_reactions: 1511
- n_metabolites: 1099
- n_genes: 727
- sha256: 2940de8c6e0497e446bf1d220ad34ee2351831e0de04b7ab466c8160b76a7af4

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Cola'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 3954
- n_experiments: 202

## Protocol

- study: transfer_v1
- arm: Uonly_atp
- arm_definition: {"transforms": ["add_atp_synthase_if_annotated"]}
- applied: [{"transform": "add_atp_synthase_if_annotated", "records": [{"reaction": "ATPS4rpp", "added": "adp_c + 4.0 h_p + pi_c <=> atp_c + h2o_c + 3.0 h_c", "gpr": "Echvi_0607 and Echvi_0608 and Echvi_0609 and Echvi_0610 and Echvi_0611 and Echvi_0612", "subunit_types": ["a", "alpha", "b", "beta", "c", "delta", "epsilon", "gamma"], "clusters": [["Echvi_0607", "Echvi_0608", "Echvi_0609", "Echvi_0610", "Echvi …
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
- gene_level_conditions_where_wt_grows: {"n_genes": 497, "n_conditions": 7, "n_gene_condition_pairs": 3479, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.3376206333815432, "ci95": [0.16675599091956317, 0.509010705963627]}, "aucpr_standard": {"point": 0.3622101953662247, "ci95": [0.24657255762903674, 0.480302495329251]}, "auroc_standard": {"point": 0.6745757761210638, "ci95": [0.6112987126532886, 0.7 …
- gene_map: {"model_genes": 727, "mapped": 715, "matched_by_version": 720, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 11, "browser_genes_hit": 715, "mapped_with_fitness_data": 497}
- counts: {"model_genes": 727, "model_genes_mapped": 715, "genes_with_fitness": 497, "genes_after_adjustment": 497, "conditions_total": 19, "conditions_mapped": 14, "conditions_wt_grows": 7, "medium_completion_exchanges_added": 8}
- timings_s: {"rich_medium_essentials_s": 1.9073486328125e-06, "knockout_simulation_s": 10.197486162185669, "total_s": 10.510364532470703}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
