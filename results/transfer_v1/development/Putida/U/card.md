# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Putida

Created: 2026-10-03T08:18:34Z

## Model

- model_id: Putida_transfer_v1_base
- file: models/transfer_v1/Putida_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/Putida.xml.gz)
- version_note: 
- n_reactions: 1899
- n_metabolites: 1317
- n_genes: 1292
- sha256: d35b0827b289142b79d3fd28ce75f7df7f4d032a6d8e964e2020270f6f5ad9e2

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Putida'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 4778
- n_experiments: 314

## Protocol

- study: transfer_v1
- arm: U
- arm_definition: {"transforms": ["universal_reaction_patches", "universal_model_additions", "add_atp_synthase_if_annotated"]}
- applied: [{"transform": "universal_reaction_patches", "records": [{"reaction": "CBMKr", "bounds_before": [-1000.0, 1000.0], "bounds_after": [-1000.0, 0.0]}]}, {"transform": "universal_model_additions", "records": [{"reaction": "GCALDt", "added": "gcald_e <=> gcald_c", "gpr": "", "bounds": [-1000.0, 1000.0]}, {"reaction": "EX_gcald_e", "added": "gcald_e --> ", "gpr": "", "bounds": [0.0, 1000.0]}]}, {"transf …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept/Putida_genpept_map.tsv"}
- role: development
- base_model: {"org": "Putida", "role": "development", "draft": "models/embl_pinned/Putida.xml.gz", "draft_sha256": "5b69bf07cf15fdecb7cf5a6c87e35271cd94e97c7d16747aee62f4c7594fd7c1", "universe": "external/carveme/universe_bacteria.xml.gz", "universe_sha256": "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19", "reference": {"medium": "MOPS minimal media_noCarbon", "carbon_exchange": "EX_glc__D_e …

## Leakage

- ground_truth_used_in_model_curation: yes: development organism; its phenotypes guided the correction rules under test (development cycles 1-7 and this study's rule selection); scores are retrospective
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: see docs/studies/transfer-method-v1.md and the panel selection record
- notes: ['Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 43, "n_conditions_wt_grows": 28, "conditions_with_absent_exchange": 10}
- gene_level_conditions_where_wt_grows: {"n_genes": 1038, "n_conditions": 28, "n_gene_condition_pairs": 29064, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.45152335075547984, "ci95": [0.3328097261799563, 0.5670530852023752]}, "aucpr_standard": {"point": 0.3755540266127471, "ci95": [0.27280345954097884, 0.47598807351990335]}, "auroc_standard": {"point": 0.7279130017985852, "ci95": [0.674938340154112 …
- gene_map: {"model_genes": 1292, "mapped": 1291, "matched_by_version": 1291, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 0, "browser_genes_hit": 1291, "mapped_with_fitness_data": 1038}
- counts: {"model_genes": 1292, "model_genes_mapped": 1291, "genes_with_fitness": 1038, "genes_after_adjustment": 1038, "conditions_total": 57, "conditions_mapped": 43, "conditions_wt_grows": 28, "medium_completion_exchanges_added": 7}
- timings_s: {"rich_medium_essentials_s": 1.9073486328125e-06, "knockout_simulation_s": 106.34268879890442, "total_s": 106.88873767852783}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
