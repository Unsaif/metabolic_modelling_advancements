# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism MR1

Created: 2026-09-06T07:38:24Z

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
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept/MR1_genpept_map.tsv"}
- condition_selection: expGroup == 'carbon source'; condition_2 empty or DMSO; media with a BiGG mapping; replicates averaged per condition x medium
- scoring: gene-level metrics pooled over conditions where the wild-type model grows (>= growth_threshold); condition-level recall = fraction of experimentally growing carbon sources on which the wild-type model grows

## Leakage

- ground_truth_used_in_model_curation: no for EMBL draft models (automated reconstruction from genome annotation); for iML1515: partly (E. coli curation used phenotype data); for iJN1463: partly (Nogales et al. 2020 validated against growth phenotypes and gene essentiality data)
- ground_truth_public_since: Fitness Browser releases 2015-2018 (Price et al. 2018)
- frontier_model_training_exposure: Fitness Browser tables are public and partly in training corpora; the mapping tables here are new
- held_out_recommendation: unpublished RB-TnSeq experiments, or organisms added to the Browser after the model's training cut-off
- notes: ['Draft models are untouched by any phenotype data, so this is a true prospective test of automated reconstruction.']

## Results

- condition_level: {"n_conditions_mapped": 12, "n_conditions_wt_grows": 8, "wt_growth_recall": 0.6666666666666666, "conditions_with_absent_exchange": 1}
- gene_level_conditions_where_wt_grows: {"n_genes": 713, "n_conditions": 8, "n_gene_condition_pairs": 5704, "aucpr_bernstein": {"point": 0.4765894298950845, "ci95": [0.3776764953548839, 0.584901639499142]}, "aucpr_standard": {"point": 0.5092440266238073, "ci95": [0.42903287672042933, 0.5842391396465685]}, "auroc_standard": {"point": 0.7301786390413458, "ci95": [0.6930834497481702, 0.7657546344239645]}, "mcc": {"point": 0.515327362213365 …
- gene_level_all_mapped_conditions: {"n_genes": 713, "n_conditions": 12, "n_gene_condition_pairs": 8556, "aucpr_bernstein": {"point": 0.5254158840403106, "ci95": [0.48682037180654514, 0.5678448530352725]}, "aucpr_standard": {"point": 0.3095870931786962, "ci95": [0.26319752776774463, 0.35150673106039126]}, "auroc_standard": {"point": 0.6676940524337169, "ci95": [0.6440285616865175, 0.6915945353162462]}, "mcc": {"point": 0.23358379529 …
- gene_map: {"model_genes": 895, "mapped": 891, "matched_by_version": 894, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 3, "browser_genes_hit": 891, "mapped_with_fitness_data": 713}
- counts: {"model_genes": 895, "model_genes_mapped": 891, "genes_with_fitness": 713, "genes_after_adjustment": 713, "conditions_total": 16, "conditions_mapped": 12, "conditions_wt_grows": 8, "medium_completion_exchanges_added": 7}
- timings_s: {"rich_medium_essentials_s": 4.76837158203125e-07, "knockout_simulation_s": 13.05344033241272, "total_s": 13.416937351226807}
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
