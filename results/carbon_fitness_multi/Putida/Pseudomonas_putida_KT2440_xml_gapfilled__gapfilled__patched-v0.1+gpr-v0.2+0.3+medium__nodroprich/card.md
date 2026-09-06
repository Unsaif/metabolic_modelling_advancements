# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Putida

Created: 2026-09-06T08:32:38Z

## Model

- model_id: Pseudomonas_putida_KT2440_xml_gapfilled
- file: models/gapfilled/Putida.xml.gz
- source: EMBL GEMs (CarveMe draft, Machado et al. 2018) minimally gap-filled by gembench.gapfill for growth on MOPS minimal media_noCarbon + EX_glc__D_e: added ['OXCDC', 'PHPYROX']
- version_note: 
- n_reactions: 1897
- n_metabolites: 1316
- n_genes: 1292
- sha256: 378f67c3c8d4642c135942917cd46f3b6109c53c16f7987c9751fb643ade5739

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Putida'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 4778
- n_experiments: 314
- conditions_unmapped: 1,3-Butandiol; 1,4-Butanediol; 1,5-Pentanediol; 1-Pentanol; 2-Piperidinone; 2-methyl-1-butanol; 3-methyl-3-butenol; 4-Hydroxyvalerate; Butyl stearate; Heptanoic acid; Levulinic Acid; Nonanoic acid; Tween 20; Valeric acid

## Protocol

- variant: gapfilled
- patches: {"file": "data/reference/universe_patches_v0.1.json", "model_file": null, "gpr_file": "data/reference/gpr_patches_v0.2.json,data/reference/gpr_patches_v0.3.json", "applied": [{"reaction": "CBMKr", "bounds_before": [-1000.0, 1000.0], "bounds_after": [-1000.0, 0.0]}, {"reaction": "ANPRT", "gpr_before": "NP_742587_1 or NP_746127_1", "gpr_after": "NP_742587_1", "rule": "R1", "status": "accepted", "gen …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept/Putida_genpept_map.tsv"}
- condition_selection: expGroup == 'carbon source'; condition_2 empty or DMSO; media with a BiGG mapping; replicates averaged per condition x medium
- scoring: gene-level metrics pooled over conditions where the wild-type model grows (>= growth_threshold); condition-level recall = fraction of experimentally growing carbon sources on which the wild-type model grows
- reporting_audit: {"date": "2026-09-06", "base_commit": "b1997d3", "note": "Leakage metadata corrected after review. Historical simulations and numeric results were not rerun or altered."}
- evaluation_role: retrospective_development

## Leakage

- ground_truth_used_in_model_curation: yes: these patch sets were developed and accepted/held/rejected using Fitness Browser phenotypes also used for scoring; annotation support does not remove this reuse
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: Freeze code, media rules and patches before selecting independent data. Use organisms/experiments not inspected during development, or nested evaluation that repeats every data-guided selection step. A retrospective split of already-inspected data does not restore independence.
- notes: ['Retrospective development evaluation. No independent held-out or prospective test is established.', 'Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare gene-level scores on matched genes, conditions and finite observations; report wild-type growth coverage separately.', 'Medium completion adds assumed gene-less transport; recipe presence alone does not establish uptake or its energy cost.']

## Results

- condition_level: {"n_conditions_mapped": 43, "n_conditions_wt_grows": 28, "wt_growth_recall": 0.6511627906976745, "conditions_with_absent_exchange": 10}
- gene_level_conditions_where_wt_grows: {"n_genes": 1038, "n_conditions": 28, "n_gene_condition_pairs": 29064, "aucpr_bernstein": {"point": 0.46705789812570053, "ci95": [0.3379307735671281, 0.5834587382004236]}, "aucpr_standard": {"point": 0.38415835003148724, "ci95": [0.2829147599426083, 0.4860465598393833]}, "auroc_standard": {"point": 0.7572200600965889, "ci95": [0.6947766980583319, 0.8060510663896613]}, "mcc": {"point": 0.5458884077 …
- gene_level_all_mapped_conditions: {"n_genes": 1038, "n_conditions": 43, "n_gene_condition_pairs": 44634, "aucpr_bernstein": {"point": 0.47191248504729916, "ci95": [0.44113937629653543, 0.5040144988234748]}, "aucpr_standard": {"point": 0.16074425761588018, "ci95": [0.11728586019656445, 0.20342307488950312]}, "auroc_standard": {"point": 0.662626317047222, "ci95": [0.6215050237070057, 0.6963193272633845]}, "mcc": {"point": 0.15720719 …
- gene_map: {"model_genes": 1292, "mapped": 1291, "matched_by_version": 1291, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 0, "browser_genes_hit": 1291, "mapped_with_fitness_data": 1038}
- counts: {"model_genes": 1292, "model_genes_mapped": 1291, "genes_with_fitness": 1038, "genes_after_adjustment": 1038, "conditions_total": 57, "conditions_mapped": 43, "conditions_wt_grows": 28, "medium_completion_exchanges_added": 7}
- timings_s: {"rich_medium_essentials_s": 1.430511474609375e-06, "knockout_simulation_s": 128.9878499507904, "total_s": 129.65734791755676}
- dropped_rich_medium_essentials: 0

## Warnings

- 14 of 57 conditions have no BiGG mapping
- medium 'MOPS minimal media_noCarbon' components absent from the model: ['EX_sel_e', 'EX_slnt_e', 'EX_tungs_e']
- medium 'RCH2_defined_noCarbon' components absent from the model: ['EX_sel_e', 'EX_slnt_e', 'EX_tungs_e', 'EX_lipoate_e', 'EX_fol_e', 'EX_cbl1_e', 'EX_adocbl_e']
- medium completion added exchange+uptake for: ['EX_4abz_e', 'EX_btn_e', 'EX_na1_e', 'EX_nac_e', 'EX_ni2_e', 'EX_ribflv_e', 'EX_thm_e']

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
