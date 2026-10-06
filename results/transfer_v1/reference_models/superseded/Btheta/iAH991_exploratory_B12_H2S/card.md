# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Btheta

Created: 2026-10-06T12:19:07Z

## Model

- model_id: iAH991_bigg_view
- file: models/curated/iAH991/iAH991_bigg_view.xml.gz
- source: Rebuilt iAH991 (as label iAH991); EXPLORATORY run with vitamin B12 and sulfide added to the medium
- version_note: 
- n_reactions: 1488
- n_metabolites: 1152
- n_genes: 991
- sha256: bafc8663ea5c6a30f16697be88add4abdfc7ec32ae7e3de9103324cce4cf9f8c

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Btheta'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 4055
- n_experiments: 542

## Protocol

- study: transfer_v1 reference models (post hoc)
- arm: REF_iAH991_exploratory_B12_H2S
- applied: []
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "locus tag with the underscore removed (BT_0554 -> BT0554 = Fitness Browser sysName)"}
- role: development
- medium_supplement: {"cbl1": -0.001, "h2s": -1000.0}
- exploratory: true

## Leakage

- ground_truth_used_in_model_curation: unknown for the published curated model; its authors may have used these or related phenotypes
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018
- frontier_model_training_exposure: not applicable (no AI step)
- held_out_recommendation: reference point only; not a held-out test
- notes: ["Scored post hoc with the transfer study's fixed protocol, for context in Paper 1.", "EXPLORATORY deviation from the protocol: every base medium supplemented with {'cbl1': -0.001, 'h2s': -1000.0}"]

## Results

- condition_level: {"n_conditions_mapped": 25, "n_conditions_wt_grows": 24, "conditions_with_absent_exchange": 1}
- gene_level_conditions_where_wt_grows: {"n_genes": 818, "n_conditions": 24, "n_gene_condition_pairs": 19632, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.47093505753383264, "ci95": [0.3729653706110869, 0.5737326660215637]}, "aucpr_standard": {"point": 0.3512042648850781, "ci95": [0.26715517443711645, 0.4530018868818377]}, "auroc_standard": {"point": 0.7845206469918031, "ci95": [0.7304021008246486, …
- gene_map: {"model_genes": 991, "mapped": 986, "mapped_with_fitness_data": 818}
- counts: {"model_genes": 991, "model_genes_mapped": 986, "genes_with_fitness": 818, "genes_after_adjustment": 818, "conditions_total": 47, "conditions_mapped": 25, "conditions_wt_grows": 24, "medium_completion_exchanges_added": 0}
- timings_s: {"rich_medium_essentials_s": 1.9073486328125e-06, "knockout_simulation_s": 38.400030851364136, "total_s": 38.893004417419434}

## Software

- python: 3.13.16
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.5.3
- scipy: 1.18.1
- scikit-learn: 1.9.1
