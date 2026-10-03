# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism MR1

Created: 2026-10-03T09:05:35Z

## Model

- model_id: MR1_transfer_v1_base
- file: models/transfer_v1/MR1_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/MR1.xml.gz)
- version_note: 
- n_reactions: 1978
- n_metabolites: 1363
- n_genes: 903
- sha256: 405850d4656d42f4c367ec1a8a247b4bdea549505974389bc0478a41eec5a69d

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'MR1'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 3782
- n_experiments: 176

## Protocol

- study: transfer_v1
- arm: M_R1R2
- arm_definition: {"transforms": ["universal_reaction_patches", "normalize_conjunction_rules", "universal_model_additions", "add_atp_synthase_if_annotated", "remove_menaquinol_if_no_pathway", "decisions:data/studies/transfer_v1/decisions/{org}.json#R1,R2,R2+R1"]}
- applied: [{"transform": "universal_reaction_patches", "records": [{"reaction": "CBPS", "gpr_before": "NP_716766_1 or NP_716767_1 or NP_716921_1 or (NP_716766_1 and NP_716767_1)", "gpr_after": "NP_716766_1 and NP_716767_1"}]}, {"transform": "normalize_conjunction_rules", "records": [{"reaction": "ACOATA", "gpr_before": "NP_717352_1 or NP_718359_1 or NP_718425_1 or NP_718636_2 or NP_719904_1 or NP_719907_1 o …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept/MR1_genpept_map.tsv"}
- role: development
- base_model: {"org": "MR1", "role": "development", "draft": "models/embl_pinned/MR1.xml.gz", "draft_sha256": "591d8e5b5fdd1124a4e78a947b0ad9907398c9a64f1c4028c6dee92452f1220e", "universe": "external/carveme/universe_bacteria.xml.gz", "universe_sha256": "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19", "reference": {"medium": "ShewMM_noCarbon", "carbon_exchange": "EX_lac__L_e"}, "min_growth":  …

## Leakage

- ground_truth_used_in_model_curation: yes: development organism; its phenotypes guided the correction rules under test (development cycles 1-7 and this study's rule selection); scores are retrospective
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: see docs/studies/transfer-method-v1.md and the panel selection record
- notes: ['Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 12, "n_conditions_wt_grows": 9, "conditions_with_absent_exchange": 1}
- gene_level_conditions_where_wt_grows: {"n_genes": 713, "n_conditions": 9, "n_gene_condition_pairs": 6417, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.46503627316551627, "ci95": [0.36717050088481773, 0.5711129471311789]}, "aucpr_standard": {"point": 0.48038991119321706, "ci95": [0.40912757360431967, 0.558535739935594]}, "auroc_standard": {"point": 0.7080900825982301, "ci95": [0.6739806355653707,  …
- gene_map: {"model_genes": 903, "mapped": 899, "matched_by_version": 894, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 3, "browser_genes_hit": 899, "mapped_with_fitness_data": 713}
- counts: {"model_genes": 903, "model_genes_mapped": 899, "genes_with_fitness": 713, "genes_after_adjustment": 713, "conditions_total": 16, "conditions_mapped": 12, "conditions_wt_grows": 9, "medium_completion_exchanges_added": 7}
- timings_s: {"rich_medium_essentials_s": 2.1457672119140625e-06, "knockout_simulation_s": 26.799330949783325, "total_s": 27.345039129257202}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
