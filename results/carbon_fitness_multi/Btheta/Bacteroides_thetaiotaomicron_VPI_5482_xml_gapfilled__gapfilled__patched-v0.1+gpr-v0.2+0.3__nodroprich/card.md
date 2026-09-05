# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Btheta

Created: 2026-09-05T20:02:40Z

## Model

- model_id: Bacteroides_thetaiotaomicron_VPI_5482_xml_gapfilled
- file: models/gapfilled/Btheta.xml.gz
- source: EMBL GEMs (CarveMe draft, Machado et al. 2018) minimally gap-filled by gembench.gapfill for growth on Varel_Bryant_medium + EX_glc__D_e: added ['ASPO2y', 'DHORD6', 'FE3t', 'OCBT_2', 'PDX5PO2']
- version_note: 
- n_reactions: 1492
- n_metabolites: 1066
- n_genes: 676
- sha256: ce97c05e8bed2f682cb9ef744c8f5a447d536f6555fff2a5ca6bfa62d00436b6

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Btheta'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 4055
- n_experiments: 542
- conditions_unmapped: (+)-Arabinogalactan; 1,4-B-D-Galactobiose; 6-O-Acetyl-D-glucose; Agro_defined_trehalose; Chondroitin sulfate A sodium salt from bovine trachea; D-Leucrose; Fructooligosaccharides (FOS); Heparin sodium salt from porcine intestinal mucosa; Hyaluronic acid sodium salt from Streptococcus equi; Isomaltose; Lactulose; Levan - from Erwinia herbicola; Maltitol; Mannan from Saccharomyces cerevisiae; Pectin; Rhamnogalacturonan - from potato; Supernatant; Agrobacterium rhizogenes K599 grown in Agro_defined_trehalose (~7.5 mM 3-keto-trehalose); a-Cyclodextrin; amylopectin from maize; dextran, Mw ~200,000; palatinose hydrate; polygalacturonic acid

## Protocol

- variant: gapfilled
- patches: {"file": "data/reference/universe_patches_v0.1.json", "gpr_file": "data/reference/gpr_patches_v0.2.json,data/reference/gpr_patches_v0.3.json", "applied": [{"reaction": "CBPS", "gpr_before": "NP_809470_1 or (NP_809469_1 and NP_809470_1)", "gpr_after": "NP_809469_1 and NP_809470_1"}, {"reaction": "ACLS", "gpr_before": "NP_810990_1 or (NP_809485_1 and NP_810990_1) or (NP_810989_1 and NP_810990_1)", " …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept/Btheta_genpept_map.tsv"}
- condition_selection: expGroup == 'carbon source'; condition_2 empty or DMSO; media with a BiGG mapping; replicates averaged per condition x medium
- scoring: gene-level metrics pooled over conditions where the wild-type model grows (>= growth_threshold); condition-level recall = fraction of experimentally growing carbon sources on which the wild-type model grows

## Leakage

- ground_truth_used_in_model_curation: no for EMBL draft models (automated reconstruction from genome annotation); for iML1515: partly (E. coli curation used phenotype data); for iJN1463: partly (Nogales et al. 2020 validated against growth phenotypes and gene essentiality data)
- ground_truth_public_since: Fitness Browser releases 2015-2018 (Price et al. 2018)
- frontier_model_training_exposure: Fitness Browser tables are public and partly in training corpora; the mapping tables here are new
- held_out_recommendation: unpublished RB-TnSeq experiments, or organisms added to the Browser after the model's training cut-off
- notes: ['Draft models are untouched by any phenotype data, so this is a true prospective test of automated reconstruction.']

## Results

- condition_level: {"n_conditions_mapped": 25, "n_conditions_wt_grows": 14, "wt_growth_recall": 0.56, "conditions_with_absent_exchange": 9}
- gene_level_conditions_where_wt_grows: {"n_genes": 508, "n_conditions": 14, "n_gene_condition_pairs": 7112, "aucpr_bernstein": {"point": 0.5490979085625025, "ci95": [0.4300416812523127, 0.6584502557775443]}, "aucpr_standard": {"point": 0.48682003074364594, "ci95": [0.38402118686255304, 0.5903506307960295]}, "auroc_standard": {"point": 0.777920955104316, "ci95": [0.7238457405515145, 0.8291089096522909]}, "mcc": {"point": 0.5603726909866 …
- gene_level_all_mapped_conditions: {"n_genes": 508, "n_conditions": 25, "n_gene_condition_pairs": 12700, "aucpr_bernstein": {"point": 0.6079641060518925, "ci95": [0.5700824454954703, 0.6497235021789756]}, "aucpr_standard": {"point": 0.24065894286177691, "ci95": [0.1917643212251781, 0.29255260560122176]}, "auroc_standard": {"point": 0.6456805513940511, "ci95": [0.6114663716667125, 0.6738217417597475]}, "mcc": {"point": 0.19903126403 …
- gene_map: {"model_genes": 676, "mapped": 675, "matched_by_version": 674, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 0, "browser_genes_hit": 675, "mapped_with_fitness_data": 508}
- counts: {"model_genes": 676, "model_genes_mapped": 675, "genes_with_fitness": 508, "genes_after_adjustment": 508, "conditions_total": 47, "conditions_mapped": 25, "conditions_wt_grows": 14}
- timings_s: {"rich_medium_essentials_s": 7.152557373046875e-07, "knockout_simulation_s": 19.808263063430786, "total_s": 20.133124589920044}
- dropped_rich_medium_essentials: 0

## Warnings

- 22 of 47 conditions have no BiGG mapping
- medium 'Varel_Bryant_medium' components absent from the model: ['EX_hco3_e', 'EX_na1_e', 'EX_mobd_e', 'EX_ni2_e', 'EX_sel_e', 'EX_slnt_e', 'EX_tungs_e', 'EX_pheme_e', 'EX_met__L_e', 'EX_cys__L_e']

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
