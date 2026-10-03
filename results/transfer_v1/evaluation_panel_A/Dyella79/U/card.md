# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Dyella79

Created: 2026-10-03T10:20:10Z

## Model

- model_id: Dyella79_transfer_v1_base
- file: models/transfer_v1/Dyella79_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/Dyella79.xml.gz)
- version_note: 
- n_reactions: 1715
- n_metabolites: 1211
- n_genes: 782
- sha256: 2b44f8c3f2ac7e46bbe757fa432b76ae187e705940e5823adaa288b1bdaa271f

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Dyella79'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 3700
- n_experiments: 69

## Protocol

- study: transfer_v1
- arm: U
- arm_definition: {"transforms": ["universal_reaction_patches", "universal_model_additions", "add_atp_synthase_if_annotated"]}
- applied: [{"transform": "universal_reaction_patches", "records": [{"reaction": "CBPS", "gpr_before": "WP_026634352_1 or WP_026634353_1 or (WP_026634352_1 and WP_026634353_1)", "gpr_after": "WP_026634352_1 and WP_026634353_1"}]}, {"transform": "universal_model_additions", "records": [{"reaction": "GCALDt", "added": "gcald_e <=> gcald_c", "gpr": "", "bounds": [-1000.0, 1000.0]}, {"reaction": "EX_gcald_e", "a …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/studies/transfer_v1/media_bigg.tsv
- carbon_source_mapping: data/studies/transfer_v1/carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept_panel/Dyella79_genpept_map.tsv"}
- role: evaluation_panel_A
- base_model: {"org": "Dyella79", "role": "evaluation_panel_A", "draft": "models/embl_pinned/Dyella79.xml.gz", "draft_sha256": "7a52ae165d3cd90609f4bb6cc6eab10a0da3a2e388c279cd1d687afb48af2250", "universe": "external/carveme/universe_bacteria.xml.gz", "universe_sha256": "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19", "reference": {"medium": "RCH2_defined_noCarbon", "carbon_exchange": "EX_glc …

## Leakage

- ground_truth_used_in_model_curation: no for the correction rules: frozen before this organism's fitness tables were downloaded (see the study freeze manifests); the base model's gap-fill uses observed wild-type growth on one reference condition
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: see docs/studies/transfer-method-v1.md and the panel selection record
- notes: ['Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 9, "n_conditions_wt_grows": 6, "conditions_with_absent_exchange": 1}
- gene_level_conditions_where_wt_grows: {"n_genes": 542, "n_conditions": 6, "n_gene_condition_pairs": 3252, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.46777346871871817, "ci95": [0.35552818715439993, 0.5872509802597755]}, "aucpr_standard": {"point": 0.45436186176368315, "ci95": [0.35531994331125133, 0.5567442948970653]}, "auroc_standard": {"point": 0.7210195226998681, "ci95": [0.6672072109004502, …
- gene_map: {"model_genes": 782, "mapped": 747, "matched_by_version": 781, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 34, "browser_genes_hit": 747, "mapped_with_fitness_data": 542}
- counts: {"model_genes": 782, "model_genes_mapped": 747, "genes_with_fitness": 542, "genes_after_adjustment": 542, "conditions_total": 10, "conditions_mapped": 9, "conditions_wt_grows": 6, "medium_completion_exchanges_added": 6}
- timings_s: {"rich_medium_essentials_s": 9.5367431640625e-07, "knockout_simulation_s": 13.186872482299805, "total_s": 13.713182210922241}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
