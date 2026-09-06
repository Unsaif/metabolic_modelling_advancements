# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism MR1

Created: 2026-09-06T08:26:49Z

## Model

- model_id: Shewanella_oneidensis_MR_1_xml_gapfilled
- file: models/gapfilled/MR1.xml.gz
- source: EMBL GEMs (CarveMe draft, Machado et al. 2018) minimally gap-filled by gembench.gapfill for growth on ShewMM_noCarbon + EX_lac__L_e: added ['CITACt', 'LDH_L', 'ORNCD']
- version_note: 
- n_reactions: 1977
- n_metabolites: 1363
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
- patches: {"file": "data/reference/universe_patches_v0.1.json", "model_file": "data/reference/model_patches_v0.1.json", "gpr_file": "data/reference/gpr_patches_v0.2.json,data/reference/gpr_patches_v0.3.json", "applied": [{"reaction": "CBPS", "gpr_before": "NP_716766_1 or NP_716767_1 or NP_716921_1 or (NP_716766_1 and NP_716767_1)", "gpr_after": "NP_716766_1 and NP_716767_1"}, {"reaction": "GCALDt", "added": …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
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
- notes: ['Retrospective development evaluation. No independent held-out or prospective test is established.', 'Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare gene-level scores on matched genes, conditions and finite observations; report wild-type growth coverage separately.', 'Medium completion adds assumed gene-less transport; recipe presence alone does not establish uptake or its energy cost.']

## Results

- condition_level: {"n_conditions_mapped": 12, "n_conditions_wt_grows": 8, "wt_growth_recall": 0.6666666666666666, "conditions_with_absent_exchange": 1}
- gene_level_conditions_where_wt_grows: {"n_genes": 713, "n_conditions": 8, "n_gene_condition_pairs": 5704, "aucpr_bernstein": {"point": 0.5236184053778908, "ci95": [0.4215280684279829, 0.6314691004766045]}, "aucpr_standard": {"point": 0.5054801568537552, "ci95": [0.42718233200907185, 0.5801871436946239]}, "auroc_standard": {"point": 0.7295560005139383, "ci95": [0.6915490386248643, 0.7661546875102443]}, "mcc": {"point": 0.53508078138018 …
- gene_level_all_mapped_conditions: {"n_genes": 713, "n_conditions": 12, "n_gene_condition_pairs": 8556, "aucpr_bernstein": {"point": 0.5397084290873828, "ci95": [0.5012174047716558, 0.5796439452084359]}, "aucpr_standard": {"point": 0.3067449366432, "ci95": [0.26269293978275166, 0.3501940199069169]}, "auroc_standard": {"point": 0.6635336624767058, "ci95": [0.6415157257024027, 0.6852775129898715]}, "mcc": {"point": 0.2442678633199769 …
- gene_map: {"model_genes": 895, "mapped": 891, "matched_by_version": 894, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 3, "browser_genes_hit": 891, "mapped_with_fitness_data": 713}
- counts: {"model_genes": 895, "model_genes_mapped": 891, "genes_with_fitness": 713, "genes_after_adjustment": 713, "conditions_total": 16, "conditions_mapped": 12, "conditions_wt_grows": 8, "medium_completion_exchanges_added": 7}
- timings_s: {"rich_medium_essentials_s": 1.430511474609375e-06, "knockout_simulation_s": 26.214118480682373, "total_s": 26.725131034851074}
- dropped_rich_medium_essentials: 0

## Warnings

- 4 of 16 conditions have no BiGG mapping
- medium 'ShewMM_noCarbon' components absent from the model: ['EX_sel_e', 'EX_slnt_e', 'EX_tungs_e', 'EX_lipoate_e', 'EX_thm_e', 'EX_pnto__R_e', 'EX_fol_e']
- medium completion added exchange+uptake for: ['EX_4abz_e', 'EX_adocbl_e', 'EX_btn_e', 'EX_cbl1_e', 'EX_na1_e', 'EX_ni2_e', 'EX_ribflv_e']

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
