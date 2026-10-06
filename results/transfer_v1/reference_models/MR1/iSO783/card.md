# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism MR1

Created: 2026-10-06T12:13:10Z

## Model

- model_id: iSO783_bigg_view
- file: models/curated/iSO783/iSO783_bigg_view.xml.gz
- source: BioModels MODEL1507180036 (Pinchuk et al. 2010), BiGG identifiers via scripts/translate_curated_model.py
- version_note: 
- n_reactions: 870
- n_metabolites: 713
- n_genes: 783
- sha256: 67dafd323dcd0d6e08cf79fc9c29634d91afd1ba5d5f70b8131e5d9c9f431cf7

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'MR1'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 3782
- n_experiments: 176

## Protocol

- study: transfer_v1 reference models (post hoc)
- arm: REF_iSO783
- applied: []
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "identity, else the same locus tag with an underscore after the letter prefix (SO0419 -> SO_0419: the Fitness Browser lists a few MR-1 loci in the RefSeq form)"}
- role: development

## Leakage

- ground_truth_used_in_model_curation: unknown for the published curated model; its authors may have used these or related phenotypes
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018
- frontier_model_training_exposure: not applicable (no AI step)
- held_out_recommendation: reference point only; not a held-out test
- notes: ["Scored post hoc with the transfer study's fixed protocol, for context in Paper 1."]

## Results

- condition_level: {"n_conditions_mapped": 12, "n_conditions_wt_grows": 9, "conditions_with_absent_exchange": 1}
- gene_level_conditions_where_wt_grows: {"n_genes": 570, "n_conditions": 9, "n_gene_condition_pairs": 5130, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.5770073157055962, "ci95": [0.47661884038666286, 0.678089300389326]}, "aucpr_standard": {"point": 0.5577501118907346, "ci95": [0.4895980186514825, 0.6292215351570161]}, "auroc_standard": {"point": 0.7542761716311869, "ci95": [0.7163025080545286, 0.7 …
- gene_map: {"model_genes": 783, "mapped": 782, "mapped_with_fitness_data": 570}
- counts: {"model_genes": 783, "model_genes_mapped": 782, "genes_with_fitness": 570, "genes_after_adjustment": 570, "conditions_total": 16, "conditions_mapped": 12, "conditions_wt_grows": 9, "medium_completion_exchanges_added": 5}
- timings_s: {"rich_medium_essentials_s": 2.1457672119140625e-06, "knockout_simulation_s": 6.63133978843689, "total_s": 6.8523406982421875}

## Software

- python: 3.13.16
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.5.3
- scipy: 1.18.1
- scikit-learn: 1.9.1
