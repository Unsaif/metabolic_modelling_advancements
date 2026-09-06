# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Keio

Created: 2026-09-06T07:50:37Z

## Model

- model_id: iML1515
- file: models/iML1515.xml
- source: BiGG iML1515 (Monk et al. 2017) via github.com/dbernste/E_coli_GEM_validation Models/
- version_note: 
- n_reactions: 2714
- n_metabolites: 1878
- n_genes: 1516
- sha256: 9c772d44ca43350e40dc7ee86c7aa148796856be1eea45e5406c6df8f7dcde28

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Keio'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 3789
- n_experiments: 168
- conditions_unmapped: 

## Protocol

- variant: shipped
- patches: {"file": "data/reference/universe_patches_v0.1.json", "model_file": "data/reference/model_patches_v0.1.json", "gpr_file": null, "applied": [{"reaction": "CBMKr", "bounds_before": [-1000.0, 1000.0], "bounds_after": [-1000.0, 0.0]}, {"reaction": "URIC", "bounds_before": [0.0, 1000.0], "bounds_after": [0.0, 0.0]}, {"reaction": "GCALDt", "added": "gcald_e <=> gcald_c", "gpr": "", "bounds": [-1000.0, 1 …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": true, "rich_medium_uptake": -1000.0, "knockout_genes": ["b0062", "b0063", "b3903", "b3904", "b0061", "b0344", "b3902"], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": false, "medium_completion_exclude": ["pnto__R", "fol"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "identity (BiGG gene ids are locus tags)"}
- condition_selection: expGroup == 'carbon source'; condition_2 empty or DMSO; media with a BiGG mapping; replicates averaged per condition x medium
- scoring: gene-level metrics pooled over conditions where the wild-type model grows (>= growth_threshold); condition-level recall = fraction of experimentally growing carbon sources on which the wild-type model grows
- reporting_audit: {"date": "2026-09-06", "base_commit": "b1997d3", "note": "Leakage metadata corrected after review. Historical simulations and numeric results were not rerun or altered."}
- evaluation_role: retrospective_development

## Leakage

- ground_truth_used_in_model_curation: yes: these patch sets were developed and accepted/held/rejected using Fitness Browser phenotypes also used for scoring; annotation support does not remove this reuse
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: Freeze code, media rules and patches before selecting independent data. Use organisms/experiments not inspected during development, or nested evaluation that repeats every data-guided selection step. A retrospective split of already-inspected data does not restore independence.
- notes: ['Retrospective development evaluation. No independent held-out or prospective test is established.', 'Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare gene-level scores on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 32, "n_conditions_wt_grows": 32, "wt_growth_recall": 1.0, "conditions_with_absent_exchange": 0}
- gene_level_conditions_where_wt_grows: {"n_genes": 1329, "n_conditions": 32, "n_gene_condition_pairs": 42528, "aucpr_bernstein": {"point": 0.6288231372296793, "ci95": [0.5166895301837996, 0.7168712114360852]}, "aucpr_standard": {"point": 0.49969683766365103, "ci95": [0.41819981351318797, 0.5739199524397489]}, "auroc_standard": {"point": 0.8035248186308137, "ci95": [0.7670020296109888, 0.8348994188701666]}, "mcc": {"point": 0.6108788863 …
- gene_level_all_mapped_conditions: {"n_genes": 1329, "n_conditions": 32, "n_gene_condition_pairs": 42528, "aucpr_bernstein": {"point": 0.6288231372296793, "ci95": [0.5166895301837996, 0.7168712114360852]}, "aucpr_standard": {"point": 0.49969683766365103, "ci95": [0.41819981351318797, 0.5739199524397489]}, "auroc_standard": {"point": 0.8035248186308137, "ci95": [0.7670020296109888, 0.8348994188701666]}, "mcc": {"point": 0.6108788863 …
- gene_map: {"model_genes": 1516, "mapped": 1515, "mapped_with_fitness_data": 1339}
- counts: {"model_genes": 1516, "model_genes_mapped": 1515, "genes_with_fitness": 1339, "genes_after_adjustment": 1329, "conditions_total": 32, "conditions_mapped": 32, "conditions_wt_grows": 32, "medium_completion_exchanges_added": 0}
- timings_s: {"rich_medium_essentials_s": 6.66391396522522, "knockout_simulation_s": 177.42419147491455, "total_s": 184.6200442314148}
- dropped_rich_medium_essentials: 10

## Warnings

- 0 of 32 conditions have no BiGG mapping
- medium 'MOPS Rich Defined media_noCarbon' components absent from the model: ['EX_4abz_e', 'EX_4hbz_e', 'EX_23dhb_e']
- medium 'RCH2_defined_noCarbon' components absent from the model: ['EX_4abz_e', 'EX_ribflv_e', 'EX_fol_e']

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
