# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism MR1

Created: 2026-09-06T10:02:53Z

## Model

- model_id: Shewanella_oneidensis_MR_1_xml_gapfilled
- file: models/gapfilled/MR1.xml.gz
- source: EMBL GEMs (CarveMe draft, Machado et al. 2018) minimally gap-filled by gembench.gapfill for growth on ShewMM_noCarbon + EX_lac__L_e: added ['CITACt', 'LDH_L', 'ORNCD']
- version_note:
- n_reactions: 1980
- n_metabolites: 1363
- n_genes: 903
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
- patches: {"file": "data/reference/universe_patches_v0.1.json", "model_file": "data/reference/model_patches_v0.3.json", "gpr_file": "data/reference/gpr_patches_v0.2.json,data/reference/gpr_patches_v0.3.json,data/reference/gpr_patches_v0.4.json", "applied": [{"reaction": "CBPS", "gpr_before": "NP_716766_1 or NP_716767_1 or NP_716921_1 or (NP_716766_1 and NP_716767_1)", "gpr_after": "NP_716766_1 and NP_716767 …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept/MR1_genpept_map.tsv"}
- condition_selection: expGroup == 'carbon source'; condition_2 empty or DMSO; media with a BiGG mapping; replicates averaged per condition x medium
- scoring: gene-level metrics pooled over conditions where the wild-type model grows (>= growth_threshold); condition-level recall = fraction of experimentally growing carbon sources on which the wild-type model grows
- input_sha256: {"data/fitness_browser/MR1/experiments.tsv": "5c16b10c64bbafa9cdcfa914e7683c505ecd22785f97f4aa56e67eead1955ef0", "data/fitness_browser/MR1/fit_logratios.tsv": "3d6677312ff34d6dcc2cb3966181fe10bfc3806f59b0fe42353bead8f584d2d8", "data/fitness_browser/MR1/genes.tsv": "c4ee5eb901accc6459bf50ba2f55a3c4f9772928f9c663cadf9a513fcc733d1d", "data/genpept/MR1_genpept_map.tsv": "e380ac438948d5c56bea457ee27cd3 …
- evaluation_role: retrospective_development

## Leakage

- ground_truth_used_in_model_curation: yes: these patch sets were developed and accepted/held/rejected using Fitness Browser phenotypes also used for scoring; annotation support does not remove this reuse
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: Freeze code, media rules and patches before selecting independent data. Use organisms/experiments not inspected during development, or nested evaluation that repeats every data-guided selection step. A retrospective split of already-inspected data does not restore independence.
- notes: ['Retrospective development evaluation. No independent held-out or prospective test is established.', 'Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare gene-level scores on matched genes, conditions and finite observations; report wild-type growth coverage separately.', 'Medium completion adds assumed gene-less transport; recipe presence alone does not establish uptake or its energy cost.']

## Results

- condition_level: {"n_conditions_mapped": 12, "n_conditions_wt_grows": 11, "wt_growth_recall": 0.9166666666666666, "conditions_with_absent_exchange": 1}
- gene_level_conditions_where_wt_grows: {"n_genes": 713, "n_conditions": 11, "n_gene_condition_pairs": 7843, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.5319703445081962, "ci95": [0.43273080730656244, 0.6392777514958424]}, "aucpr_standard": {"point": 0.5119530773659523, "ci95": [0.4372421264502152, 0.5848258252437087]}, "auroc_standard": {"point": 0.7273485261500038, "ci95": [0.6929196127126783, 0 …
- gene_level_all_mapped_conditions: {"n_genes": 713, "n_conditions": 12, "n_gene_condition_pairs": 8556, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.449802477286639, "ci95": [0.3763066918081864, 0.5304547193922844]}, "aucpr_standard": {"point": 0.4369917759936935, "ci95": [0.3760125909514763, 0.5012522540336545]}, "auroc_standard": {"point": 0.7205667275182983, "ci95": [0.6890393555233982, 0.7 …
- gene_map: {"model_genes": 903, "mapped": 899, "matched_by_version": 894, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 3, "browser_genes_hit": 899, "mapped_with_fitness_data": 713}
- counts: {"model_genes": 903, "model_genes_mapped": 899, "genes_with_fitness": 713, "genes_after_adjustment": 713, "conditions_total": 16, "conditions_mapped": 12, "conditions_wt_grows": 11, "medium_completion_exchanges_added": 7}
- timings_s: {"rich_medium_essentials_s": 2.1457672119140625e-06, "knockout_simulation_s": 53.79000926017761, "total_s": 54.52835297584534}
- dropped_rich_medium_essentials: 0

## Warnings

- 4 of 16 conditions have no BiGG mapping
- medium 'ShewMM_noCarbon' components absent from the model: ['EX_sel_e', 'EX_slnt_e', 'EX_tungs_e', 'EX_lipoate_e', 'EX_thm_e', 'EX_pnto__R_e', 'EX_fol_e']
- medium completion added exchange+uptake for: ['EX_4abz_e', 'EX_adocbl_e', 'EX_btn_e', 'EX_cbl1_e', 'EX_na1_e', 'EX_ni2_e', 'EX_ribflv_e']

## Software

- python: 3.14.3
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.5.2
- scipy: 1.18.1
- scikit-learn: 1.9.0
