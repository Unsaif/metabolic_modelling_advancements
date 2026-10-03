# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Btheta

Created: 2026-10-03T08:02:48Z

## Model

- model_id: Btheta_transfer_v1_base
- file: models/transfer_v1/Btheta_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/Btheta.xml.gz)
- version_note: 
- n_reactions: 1493
- n_metabolites: 1066
- n_genes: 683
- sha256: f4c4fd384ebd27ff390c51a2b084ca2a4662b46e9d9bdd68846f25128948a5a8

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Btheta'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 4055
- n_experiments: 542

## Protocol

- study: transfer_v1
- arm: Uonly_atp
- arm_definition: {"transforms": ["add_atp_synthase_if_annotated"]}
- applied: [{"transform": "add_atp_synthase_if_annotated", "records": [{"reaction": "ATPS4rpp", "added": "adp_c + 4.0 h_p + pi_c <=> atp_c + h2o_c + 3.0 h_c", "gpr": "BT0711 and BT0712 and BT0714 and BT0715 and BT0716 and BT0717 and BT0718 and BT0719", "subunit_types": ["a", "alpha", "b", "beta", "c", "delta", "epsilon", "gamma"], "clusters": [["BT0711", "BT0712", "BT0714", "BT0715", "BT0716", "BT0717", "BT0 …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept/Btheta_genpept_map.tsv"}
- role: development
- base_model: {"org": "Btheta", "role": "development", "draft": "models/embl_pinned/Btheta.xml.gz", "draft_sha256": "0d38cb19d7e2deba4224c12eb114ec4b9f812678284091dd6a30cb5f394cd985", "universe": "external/carveme/universe_bacteria.xml.gz", "universe_sha256": "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19", "reference": {"medium": "Varel_Bryant_medium", "carbon_exchange": "EX_glc__D_e"}, "min …

## Leakage

- ground_truth_used_in_model_curation: yes: development organism; its phenotypes guided the correction rules under test (development cycles 1-7 and this study's rule selection); scores are retrospective
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: see docs/studies/transfer-method-v1.md and the panel selection record
- notes: ['Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 25, "n_conditions_wt_grows": 14, "conditions_with_absent_exchange": 9}
- gene_level_conditions_where_wt_grows: {"n_genes": 507, "n_conditions": 14, "n_gene_condition_pairs": 7098, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.4706600775859737, "ci95": [0.34893343573431246, 0.5969653689780478]}, "aucpr_standard": {"point": 0.4171687485468481, "ci95": [0.31331848174097476, 0.5359505500254621]}, "auroc_standard": {"point": 0.7438855121424418, "ci95": [0.6900710570758748,  …
- gene_map: {"model_genes": 683, "mapped": 682, "matched_by_version": 674, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 0, "browser_genes_hit": 682, "mapped_with_fitness_data": 507}
- counts: {"model_genes": 683, "model_genes_mapped": 682, "genes_with_fitness": 507, "genes_after_adjustment": 507, "conditions_total": 47, "conditions_mapped": 25, "conditions_wt_grows": 14, "medium_completion_exchanges_added": 4}
- timings_s: {"rich_medium_essentials_s": 1.9073486328125e-06, "knockout_simulation_s": 20.861642360687256, "total_s": 21.15462350845337}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
