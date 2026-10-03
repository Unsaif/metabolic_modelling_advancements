# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism RPal_CGA009

Created: 2026-10-03T10:28:37Z

## Model

- model_id: RPal_CGA009_transfer_v1_base
- file: models/transfer_v1/RPal_CGA009_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/RPal_CGA009.xml.gz)
- version_note: 
- n_reactions: 1843
- n_metabolites: 1324
- n_genes: 955
- sha256: 98523291d5d082d65d580cfdd57565b171a983664756c30d68030cefa56ca1e9

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'RPal_CGA009'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 4322
- n_experiments: 208

## Protocol

- study: transfer_v1
- arm: Uonly_atp
- arm_definition: {"transforms": ["add_atp_synthase_if_annotated"]}
- applied: [{"transform": "add_atp_synthase_if_annotated", "records": [{"reaction": "ATPS4rpp", "added": "adp_c + 4.0 h_p + pi_c <=> atp_c + h2o_c + 3.0 h_c", "gpr": "TX73_000910 and TX73_000915 and TX73_000920 and TX73_000925 and TX73_000930 and TX73_004345 and TX73_004350 and TX73_004355 and TX73_004360", "subunit_types": ["a", "alpha", "b", "beta", "c", "delta", "epsilon", "gamma"], "clusters": [["TX73_00 …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/studies/transfer_v1/media_bigg.tsv
- carbon_source_mapping: data/studies/transfer_v1/carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept_panel/RPal_CGA009_genpept_map.tsv"}
- role: evaluation_panel_A
- base_model: {"org": "RPal_CGA009", "role": "evaluation_panel_A", "draft": "models/embl_pinned/RPal_CGA009.xml.gz", "draft_sha256": "1cb49f1a06b244929b4f83943d433ef641546fa36c86e54477c73a9d3ad7f32f", "universe": "external/carveme/universe_bacteria.xml.gz", "universe_sha256": "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19", "reference": {"medium": "PM", "carbon_exchange": "EX_succ_e", "condit …

## Leakage

- ground_truth_used_in_model_curation: no for the correction rules: frozen before this organism's fitness tables were downloaded (see the study freeze manifests); the base model's gap-fill uses observed wild-type growth on one reference condition
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: see docs/studies/transfer-method-v1.md and the panel selection record
- notes: ['Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 7, "n_conditions_wt_grows": 3, "conditions_with_absent_exchange": 0}
- gene_level_conditions_where_wt_grows: {"n_genes": 679, "n_conditions": 3, "n_gene_condition_pairs": 2037, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.06112900631475161, "ci95": [0.023381038526289196, 0.13650841225104743]}, "aucpr_standard": {"point": 0.05097892610843541, "ci95": [0.010240994011708209, 0.1801314575068744]}, "auroc_standard": {"point": 0.653521019010825, "ci95": [0.463091395247973 …
- gene_map: {"model_genes": 955, "mapped": 933, "matched_by_version": 945, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 21, "browser_genes_hit": 933, "mapped_with_fitness_data": 679}
- counts: {"model_genes": 955, "model_genes_mapped": 933, "genes_with_fitness": 679, "genes_after_adjustment": 679, "conditions_total": 7, "conditions_mapped": 7, "conditions_wt_grows": 3, "medium_completion_exchanges_added": 2}
- timings_s: {"rich_medium_essentials_s": 1.9073486328125e-06, "knockout_simulation_s": 7.472981691360474, "total_s": 7.844852447509766}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
