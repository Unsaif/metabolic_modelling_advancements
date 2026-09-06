# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Smeli

Created: 2026-09-06T10:39:14Z

## Model

- model_id: Sinorhizobium_meliloti_1021_xml_gapfilled
- file: models/gapfilled/Smeli.xml.gz
- source: EMBL GEMs (CarveMe draft, Machado et al. 2018) minimally gap-filled by gembench.gapfill for growth on RCH2_defined_noCarbon + EX_glc__D_e: added ['DXPS']
- version_note: 
- n_reactions: 2003
- n_metabolites: 1377
- n_genes: 1220
- sha256: aed3446b54ecac8a7a48e22c94e75fc57d32c35666a6431ae41678183cf85dfa

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Smeli'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 5133
- n_experiments: 90
- conditions_unmapped: 

## Protocol

- variant: gapfilled
- patches: {"file": "data/reference/universe_patches_v0.1.json", "model_file": "data/reference/model_patches_v0.4.json", "gpr_file": "data/reference/gpr_patches_v0.2.json,data/reference/gpr_patches_v0.3.json,data/reference/gpr_patches_v0.4.json", "applied": [{"reaction": "CBMKr", "bounds_before": [-1000.0, 1000.0], "bounds_after": [-1000.0, 0.0]}, {"reaction": "CBPS", "gpr_before": "NP_385682_1 or NP_386426_ …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept/Smeli_genpept_map.tsv"}
- condition_selection: expGroup == 'carbon source'; condition_2 empty or DMSO; media with a BiGG mapping; replicates averaged per condition x medium
- scoring: gene-level metrics pooled over conditions where the wild-type model grows (>= growth_threshold); condition-level recall = fraction of experimentally growing carbon sources on which the wild-type model grows
- input_sha256: {"data/fitness_browser/Smeli/experiments.tsv": "3c09d633f61e18afa480e69547cc2791fb2aaa93b708058689455ae295c20ea2", "data/fitness_browser/Smeli/fit_logratios.tsv": "f16dbe9b13da71a9583e4ca0a7345bc6abfb114a9aab0b3f9f4cf33a345d640e", "data/fitness_browser/Smeli/genes.tsv": "4f587c69af42b451f7d99519398ab0444cf989c7f9371087fb3691338f07798e", "data/genpept/Smeli_genpept_map.tsv": "74b32adcc980040e38cd2c …
- evaluation_role: retrospective_development

## Leakage

- ground_truth_used_in_model_curation: yes: these patch sets were developed and accepted/held/rejected using Fitness Browser phenotypes also used for scoring; annotation support does not remove this reuse
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: Freeze code, media rules and patches before selecting independent data. Use organisms/experiments not inspected during development, or nested evaluation that repeats every data-guided selection step. A retrospective split of already-inspected data does not restore independence.
- notes: ['Retrospective development evaluation. No independent held-out or prospective test is established.', 'Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare gene-level scores on matched genes, conditions and finite observations; report wild-type growth coverage separately.', 'Medium completion adds assumed gene-less transport; recipe presence alone does not establish uptake or its energy cost.']

## Results

- condition_level: {"n_conditions_mapped": 33, "n_conditions_wt_grows": 27, "wt_growth_recall": 0.8181818181818182, "conditions_with_absent_exchange": 5}
- gene_level_conditions_where_wt_grows: {"n_genes": 989, "n_conditions": 27, "n_gene_condition_pairs": 26703, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.6692229118454969, "ci95": [0.535502423483713, 0.7792823134951639]}, "aucpr_standard": {"point": 0.5899675970804311, "ci95": [0.49003748615283166, 0.6743556399342774]}, "auroc_standard": {"point": 0.7719803978749248, "ci95": [0.7185619877363413, 0 …
- gene_level_all_mapped_conditions: {"n_genes": 989, "n_conditions": 33, "n_gene_condition_pairs": 32637, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.39743010632364945, "ci95": [0.33929938253526304, 0.45213174730833655]}, "aucpr_standard": {"point": 0.321578970996812, "ci95": [0.2633349803527581, 0.3763885043450084]}, "auroc_standard": {"point": 0.7199300799511452, "ci95": [0.6771847368891343, …
- gene_map: {"model_genes": 1220, "mapped": 1219, "matched_by_version": 1195, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 0, "browser_genes_hit": 1219, "mapped_with_fitness_data": 989}
- counts: {"model_genes": 1220, "model_genes_mapped": 1219, "genes_with_fitness": 989, "genes_after_adjustment": 989, "conditions_total": 33, "conditions_mapped": 33, "conditions_wt_grows": 27, "medium_completion_exchanges_added": 6}
- timings_s: {"rich_medium_essentials_s": 9.5367431640625e-07, "knockout_simulation_s": 143.0988359451294, "total_s": 143.57371377944946}
- dropped_rich_medium_essentials: 0

## Warnings

- 0 of 33 conditions have no BiGG mapping
- medium 'RCH2_defined_noCarbon' components absent from the model: ['EX_sel_e', 'EX_slnt_e', 'EX_tungs_e', 'EX_lipoate_e', 'EX_pnto__R_e', 'EX_fol_e', 'EX_cbl1_e', 'EX_adocbl_e']
- medium completion added exchange+uptake for: ['EX_4abz_e', 'EX_btn_e', 'EX_na1_e', 'EX_nac_e', 'EX_ni2_e', 'EX_ribflv_e']

## Software

- python: 3.14.3
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.5.2
- scipy: 1.18.1
- scikit-learn: 1.9.0
