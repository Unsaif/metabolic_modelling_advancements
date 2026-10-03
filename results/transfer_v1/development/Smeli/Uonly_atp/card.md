# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Smeli

Created: 2026-10-03T08:14:28Z

## Model

- model_id: Smeli_transfer_v1_base
- file: models/transfer_v1/Smeli_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/Smeli.xml.gz)
- version_note: 
- n_reactions: 1994
- n_metabolites: 1373
- n_genes: 1205
- sha256: f2709520e35ba1dd1991b80938dfbe7c34688ac6e95178f1c314c32ab609b671

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Smeli'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 5133
- n_experiments: 90

## Protocol

- study: transfer_v1
- arm: Uonly_atp
- arm_definition: {"transforms": ["add_atp_synthase_if_annotated"]}
- applied: [{"transform": "add_atp_synthase_if_annotated", "records": [{"reaction": "ATPS4rpp", "added": "adp_c + 4.0 h_p + pi_c <=> atp_c + h2o_c + 3.0 h_c", "gpr": "SMc00871 and SMc00870 and SMc00869 and SMc00868 and SMc02502 and SMc02501 and SMc02500 and SMc02499 and SMc02498", "subunit_types": ["a", "alpha", "b", "beta", "c", "delta", "epsilon", "gamma"], "clusters": [["SMc00871", "SMc00870", "SMc00869", …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept/Smeli_genpept_map.tsv"}
- role: development
- base_model: {"org": "Smeli", "role": "development", "draft": "models/embl_pinned/Smeli.xml.gz", "draft_sha256": "9881183f0104a2966fa81e80e471c3ab26c948ac2206911e9363da0ec2341362", "universe": "external/carveme/universe_bacteria.xml.gz", "universe_sha256": "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19", "reference": {"medium": "RCH2_defined_noCarbon", "carbon_exchange": "EX_glc__D_e"}, "min …

## Leakage

- ground_truth_used_in_model_curation: yes: development organism; its phenotypes guided the correction rules under test (development cycles 1-7 and this study's rule selection); scores are retrospective
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: see docs/studies/transfer-method-v1.md and the panel selection record
- notes: ['Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 33, "n_conditions_wt_grows": 22, "conditions_with_absent_exchange": 6}
- gene_level_conditions_where_wt_grows: {"n_genes": 979, "n_conditions": 22, "n_gene_condition_pairs": 21538, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.5309188726170851, "ci95": [0.38475441784046616, 0.6756023506505255]}, "aucpr_standard": {"point": 0.45683625046898113, "ci95": [0.36015372284466585, 0.5628756161398464]}, "auroc_standard": {"point": 0.7118491954053999, "ci95": [0.6628917485704255 …
- gene_map: {"model_genes": 1205, "mapped": 1204, "matched_by_version": 1195, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 0, "browser_genes_hit": 1204, "mapped_with_fitness_data": 979}
- counts: {"model_genes": 1205, "model_genes_mapped": 1204, "genes_with_fitness": 979, "genes_after_adjustment": 979, "conditions_total": 33, "conditions_mapped": 33, "conditions_wt_grows": 22, "medium_completion_exchanges_added": 6}
- timings_s: {"rich_medium_essentials_s": 1.9073486328125e-06, "knockout_simulation_s": 89.66858291625977, "total_s": 90.3474338054657}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
