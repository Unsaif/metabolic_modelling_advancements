# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism MR1

Created: 2026-09-06T09:03:40Z

## Model

- model_id: Shewanella_oneidensis_MR_1_xml_gapfilled
- file: models/gapfilled/MR1.xml.gz
- source: EMBL GEMs (CarveMe draft, Machado et al. 2018) minimally gap-filled by gembench.gapfill for growth on ShewMM_noCarbon + EX_lac__L_e: added ['CITACt', 'LDH_L', 'ORNCD']
- version_note: 
- n_reactions: 1975
- n_metabolites: 1362
- n_genes: 895
- sha256: 7aae2a3ea5b644241d66bd719fd8a5d5a20d6e2d612ff9706f689a91fa75bdc5

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'MR1'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 3782
- n_experiments: 176
- conditions_unmapped: Gelatin; Gly-DL-Asp; Gly-Glu; Tween 20

## Protocol

- variant: gapfilled
- patches: {"file": "data/reference/universe_patches_v0.1.json", "model_file": null, "gpr_file": null, "applied": [{"reaction": "CBPS", "gpr_before": "NP_716766_1 or NP_716767_1 or NP_716921_1 or (NP_716766_1 and NP_716767_1)", "gpr_after": "NP_716766_1 and NP_716767_1"}]}
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": true, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": false, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept/MR1_genpept_map.tsv"}
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

- condition_level: {"n_conditions_mapped": 12, "n_conditions_wt_grows": 8, "wt_growth_recall": 0.6666666666666666, "conditions_with_absent_exchange": 1}
- gene_level_conditions_where_wt_grows: {"n_genes": 674, "n_conditions": 8, "n_gene_condition_pairs": 5392, "aucpr_bernstein": {"point": 0.3077000587961048, "ci95": [0.1870902501660819, 0.4406838903986626]}, "aucpr_standard": {"point": 0.37424611671276875, "ci95": [0.2984428482820715, 0.4450226487985565]}, "auroc_standard": {"point": 0.6343347261381111, "ci95": [0.5995111223826717, 0.669910701311228]}, "mcc": {"point": 0.416457899681449 …
- gene_level_all_mapped_conditions: {"n_genes": 674, "n_conditions": 12, "n_gene_condition_pairs": 8088, "aucpr_bernstein": {"point": 0.425525005536688, "ci95": [0.39728021872284974, 0.45516607856469077]}, "aucpr_standard": {"point": 0.23502899681178752, "ci95": [0.19980347274004465, 0.2728163963450868]}, "auroc_standard": {"point": 0.6106949158654801, "ci95": [0.5891820804158101, 0.6333441083601414]}, "mcc": {"point": 0.14217164943 …
- gene_map: {"model_genes": 895, "mapped": 891, "matched_by_version": 894, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 3, "browser_genes_hit": 891, "mapped_with_fitness_data": 713}
- counts: {"model_genes": 895, "model_genes_mapped": 891, "genes_with_fitness": 713, "genes_after_adjustment": 674, "conditions_total": 16, "conditions_mapped": 12, "conditions_wt_grows": 8, "medium_completion_exchanges_added": 0}
- timings_s: {"rich_medium_essentials_s": 3.6866228580474854, "knockout_simulation_s": 23.294018983840942, "total_s": 27.656171321868896}
- dropped_rich_medium_essentials: 39

## Warnings

- 4 of 16 conditions have no BiGG mapping
- medium 'ShewMM_noCarbon' components absent from the model: ['EX_na1_e', 'EX_ni2_e', 'EX_sel_e', 'EX_slnt_e', 'EX_tungs_e', 'EX_4abz_e', 'EX_lipoate_e', 'EX_ribflv_e', 'EX_thm_e', 'EX_pnto__R_e', 'EX_btn_e', 'EX_fol_e', 'EX_cbl1_e', 'EX_adocbl_e']

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
