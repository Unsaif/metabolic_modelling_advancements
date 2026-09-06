# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Smeli

Created: 2026-09-06T08:41:39Z

## Model

- model_id: Sinorhizobium_meliloti_1021_xml_gapfilled
- file: models/gapfilled/Smeli.xml.gz
- source: EMBL GEMs (CarveMe draft, Machado et al. 2018) minimally gap-filled by gembench.gapfill for growth on RCH2_defined_noCarbon + EX_glc__D_e: added ['DXPS']
- version_note: 
- n_reactions: 1993
- n_metabolites: 1373
- n_genes: 1198
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
- patches: {"file": "data/reference/universe_patches_v0.1.json", "model_file": null, "gpr_file": "data/reference/gpr_patches_v0.2.json,data/reference/gpr_patches_v0.3.json", "applied": [{"reaction": "CBMKr", "bounds_before": [-1000.0, 1000.0], "bounds_after": [-1000.0, 0.0]}, {"reaction": "CBPS", "gpr_before": "NP_385682_1 or NP_386426_1 or (NP_385682_1 and NP_386426_1)", "gpr_after": "NP_385682_1 and NP_386 …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": false, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept/Smeli_genpept_map.tsv"}
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

- condition_level: {"n_conditions_mapped": 33, "n_conditions_wt_grows": 21, "wt_growth_recall": 0.6363636363636364, "conditions_with_absent_exchange": 6}
- gene_level_conditions_where_wt_grows: {"n_genes": 981, "n_conditions": 21, "n_gene_condition_pairs": 20601, "aucpr_bernstein": {"point": 0.5791913111628567, "ci95": [0.45542323244189425, 0.6980912839651681]}, "aucpr_standard": {"point": 0.5158566562598321, "ci95": [0.40622207880624656, 0.6085400352425657]}, "auroc_standard": {"point": 0.7721567879589714, "ci95": [0.7193456148282951, 0.8207204642590549]}, "mcc": {"point": 0.61551527798 …
- gene_level_all_mapped_conditions: {"n_genes": 981, "n_conditions": 33, "n_gene_condition_pairs": 32373, "aucpr_bernstein": {"point": 0.5072946556127782, "ci95": [0.47361982534299124, 0.5407059533571629]}, "aucpr_standard": {"point": 0.24315504851642952, "ci95": [0.19336625102065744, 0.28997767152528714]}, "auroc_standard": {"point": 0.6792609863994414, "ci95": [0.6452995811590071, 0.7110719008651741]}, "mcc": {"point": 0.195592772 …
- gene_map: {"model_genes": 1198, "mapped": 1197, "matched_by_version": 1195, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 0, "browser_genes_hit": 1197, "mapped_with_fitness_data": 981}
- counts: {"model_genes": 1198, "model_genes_mapped": 1197, "genes_with_fitness": 981, "genes_after_adjustment": 981, "conditions_total": 33, "conditions_mapped": 33, "conditions_wt_grows": 21, "medium_completion_exchanges_added": 0}
- timings_s: {"rich_medium_essentials_s": 1.9073486328125e-06, "knockout_simulation_s": 92.23414444923401, "total_s": 92.97604060173035}
- dropped_rich_medium_essentials: 0

## Warnings

- 0 of 33 conditions have no BiGG mapping
- medium 'RCH2_defined_noCarbon' components absent from the model: ['EX_na1_e', 'EX_ni2_e', 'EX_sel_e', 'EX_slnt_e', 'EX_tungs_e', 'EX_4abz_e', 'EX_lipoate_e', 'EX_nac_e', 'EX_ribflv_e', 'EX_thm_e', 'EX_pnto__R_e', 'EX_btn_e', 'EX_fol_e', 'EX_cbl1_e', 'EX_adocbl_e']

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
