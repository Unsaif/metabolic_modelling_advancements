# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Dino

Created: 2026-10-03T10:18:30Z

## Model

- model_id: Dino_transfer_v1_base
- file: models/transfer_v1/Dino_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/Dino.xml.gz)
- version_note: 
- n_reactions: 1663
- n_metabolites: 1196
- n_genes: 867
- sha256: b409621713f7dc2cb1e31ce4db773426e04ae1609ab5ae4a4f034ff4d226888c

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Dino'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 3187
- n_experiments: 186

## Protocol

- study: transfer_v1
- arm: Uonly_atp
- arm_definition: {"transforms": ["add_atp_synthase_if_annotated"]}
- applied: [{"transform": "add_atp_synthase_if_annotated", "records": [{"reaction": "ATPS4rpp", "added": "adp_c + 4.0 h_p + pi_c <=> atp_c + h2o_c + 3.0 h_c", "gpr": "Dshi_0435 and Dshi_0436 and Dshi_0439 and Dshi_0440 and Dshi_0442 and Dshi_2933 and Dshi_2934 and Dshi_2935 and Dshi_2936 and Dshi_2937 and Dshi_3030", "subunit_types": ["a", "alpha", "beta", "c", "delta", "epsilon", "gamma"], "clusters": [["Ds …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/studies/transfer_v1/media_bigg.tsv
- carbon_source_mapping: data/studies/transfer_v1/carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept_panel/Dino_genpept_map.tsv"}
- role: evaluation_panel_A
- base_model: {"org": "Dino", "role": "evaluation_panel_A", "draft": "models/embl_pinned/Dino.xml.gz", "draft_sha256": "991bf8bcb2c8cc12c2c4857b9cddaeb17b003165065a7fe0b1e439745d46ce5c", "universe": "external/carveme/universe_bacteria.xml.gz", "universe_sha256": "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19", "reference": {"medium": "DinoMM_noCarbon_HighNutrient", "carbon_exchange": "EX_glc_ …

## Leakage

- ground_truth_used_in_model_curation: no for the correction rules: frozen before this organism's fitness tables were downloaded (see the study freeze manifests); the base model's gap-fill uses observed wild-type growth on one reference condition
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: see docs/studies/transfer-method-v1.md and the panel selection record
- notes: ['Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 21, "n_conditions_wt_grows": 7, "conditions_with_absent_exchange": 11}
- gene_level_conditions_where_wt_grows: {"n_genes": 584, "n_conditions": 7, "n_gene_condition_pairs": 4088, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.28972226111362503, "ci95": [0.12104031820174718, 0.46642819267026897]}, "aucpr_standard": {"point": 0.25063420778079354, "ci95": [0.10820519081692037, 0.42461259209336266]}, "auroc_standard": {"point": 0.666422452091032, "ci95": [0.5727575065277887 …
- gene_map: {"model_genes": 867, "mapped": 829, "matched_by_version": 855, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 37, "browser_genes_hit": 829, "mapped_with_fitness_data": 584}
- counts: {"model_genes": 867, "model_genes_mapped": 829, "genes_with_fitness": 584, "genes_after_adjustment": 584, "conditions_total": 21, "conditions_mapped": 21, "conditions_wt_grows": 7, "medium_completion_exchanges_added": 5}
- timings_s: {"rich_medium_essentials_s": 1.1920928955078125e-06, "knockout_simulation_s": 14.683675527572632, "total_s": 15.09714150428772}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
