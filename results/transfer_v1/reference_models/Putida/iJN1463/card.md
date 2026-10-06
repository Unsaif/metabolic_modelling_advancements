# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Putida

Created: 2026-10-06T09:11:24Z

## Model

- model_id: iJN1463
- file: models/bigg/iJN1463.xml
- source: BiGG iJN1463 (Nogales et al. 2020), downloaded 5 Sept 2026
- version_note: 
- n_reactions: 2927
- n_metabolites: 2153
- n_genes: 1462
- sha256: d573833328ffae0dfa752a1fa3262ed939ed5862288beab287fca30d0fefb4a1

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Putida'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 4778
- n_experiments: 314

## Protocol

- study: transfer_v1 reference models (post hoc)
- arm: REF_iJN1463
- applied: []
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "identity (BiGG gene ids are locus tags = Fitness Browser sysName)"}
- role: development

## Leakage

- ground_truth_used_in_model_curation: unknown for the published curated model; its authors may have used these or related phenotypes
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018
- frontier_model_training_exposure: not applicable (no AI step)
- held_out_recommendation: reference point only; not a held-out test
- notes: ["Scored post hoc with the transfer study's fixed protocol, for context in Paper 1."]

## Results

- condition_level: {"n_conditions_mapped": 43, "n_conditions_wt_grows": 36, "conditions_with_absent_exchange": 7}
- gene_level_conditions_where_wt_grows: {"n_genes": 1148, "n_conditions": 36, "n_gene_condition_pairs": 41328, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.4491030986431062, "ci95": [0.3633416915536046, 0.5272681052675446]}, "aucpr_standard": {"point": 0.29222026290890024, "ci95": [0.2170246008392563, 0.37065545955199236]}, "auroc_standard": {"point": 0.7961049096813735, "ci95": [0.7440518731201029 …
- gene_map: {"model_genes": 1462, "mapped": 1440, "mapped_with_fitness_data": 1148}
- counts: {"model_genes": 1462, "model_genes_mapped": 1440, "genes_with_fitness": 1148, "genes_after_adjustment": 1148, "conditions_total": 57, "conditions_mapped": 43, "conditions_wt_grows": 36, "medium_completion_exchanges_added": 4}
- timings_s: {"rich_medium_essentials_s": 2.86102294921875e-06, "knockout_simulation_s": 195.71993350982666, "total_s": 196.62029719352722}

## Software

- python: 3.13.16
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.5.3
- scipy: 1.18.1
- scikit-learn: 1.9.1
