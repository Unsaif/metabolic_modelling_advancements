# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Smeli

Created: 2026-10-03T08:04:24Z

## Model

- model_id: Smeli_transfer_v1_base
- file: models/transfer_v1/Smeli_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/Smeli.xml.gz)
- version_note: 
- n_reactions: 1993
- n_metabolites: 1373
- n_genes: 1196
- sha256: f2709520e35ba1dd1991b80938dfbe7c34688ac6e95178f1c314c32ab609b671

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Smeli'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 5133
- n_experiments: 90

## Protocol

- study: transfer_v1
- arm: B0
- arm_definition: {"transforms": []}
- applied: []
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

- condition_level: {"n_conditions_mapped": 33, "n_conditions_wt_grows": 21, "conditions_with_absent_exchange": 6}
- gene_level_conditions_where_wt_grows: {"n_genes": 979, "n_conditions": 21, "n_gene_condition_pairs": 20559, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.5245603965015921, "ci95": [0.3810956287237439, 0.6708572450960442]}, "aucpr_standard": {"point": 0.45055233791160165, "ci95": [0.3536312329702701, 0.5552005930956341]}, "auroc_standard": {"point": 0.7117256694215942, "ci95": [0.6634342089272258,  …
- gene_map: {"model_genes": 1196, "mapped": 1195, "matched_by_version": 1195, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 0, "browser_genes_hit": 1195, "mapped_with_fitness_data": 979}
- counts: {"model_genes": 1196, "model_genes_mapped": 1195, "genes_with_fitness": 979, "genes_after_adjustment": 979, "conditions_total": 33, "conditions_mapped": 33, "conditions_wt_grows": 21, "medium_completion_exchanges_added": 6}
- timings_s: {"rich_medium_essentials_s": 3.0994415283203125e-06, "knockout_simulation_s": 80.52700924873352, "total_s": 81.1760368347168}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
