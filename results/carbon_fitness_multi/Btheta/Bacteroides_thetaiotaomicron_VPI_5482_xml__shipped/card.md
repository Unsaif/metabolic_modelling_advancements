# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Btheta

Created: 2026-09-05T12:15:48Z

## Model

- model_id: Bacteroides_thetaiotaomicron_VPI_5482_xml
- file: models/embl/Bacteroides_thetaiotaomicron_VPI_5482.xml.gz
- source: EMBL GEMs (CarveMe draft, Machado et al. 2018; github.com/cdanielmachado/embl_gems)
- version_note: 
- n_reactions: 1487
- n_metabolites: 1066
- n_genes: 675
- sha256: 0d38cb19d7e2deba4224c12eb114ec4b9f812678284091dd6a30cb5f394cd985

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Btheta'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 4055
- n_experiments: 542
- conditions_unmapped: (+)-Arabinogalactan; 1,4-B-D-Galactobiose; 6-O-Acetyl-D-glucose; Agro_defined_trehalose; Chondroitin sulfate A sodium salt from bovine trachea; D-Leucrose; Fructooligosaccharides (FOS); Heparin sodium salt from porcine intestinal mucosa; Hyaluronic acid sodium salt from Streptococcus equi; Isomaltose; Lactulose; Levan - from Erwinia herbicola; Maltitol; Mannan from Saccharomyces cerevisiae; Pectin; Rhamnogalacturonan - from potato; Supernatant; Agrobacterium rhizogenes K599 grown in Agro_defined_trehalose (~7.5 mM 3-keto-trehalose); a-Cyclodextrin; amylopectin from maize; dextran, Mw ~200,000; palatinose hydrate; polygalacturonic acid

## Protocol

- variant: shipped
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": true, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept/Btheta_genpept_map.tsv"}
- condition_selection: expGroup == 'carbon source'; condition_2 empty or DMSO; media with a BiGG mapping; replicates averaged per condition x medium
- scoring: gene-level metrics pooled over conditions where the wild-type model grows (>= growth_threshold); condition-level recall = fraction of experimentally growing carbon sources on which the wild-type model grows

## Leakage

- ground_truth_used_in_model_curation: no for EMBL draft models (automated reconstruction from genome annotation); for iML1515: partly (E. coli curation used phenotype data)
- ground_truth_public_since: Fitness Browser releases 2015-2018 (Price et al. 2018)
- frontier_model_training_exposure: Fitness Browser tables are public and partly in training corpora; the mapping tables here are new
- held_out_recommendation: unpublished RB-TnSeq experiments, or organisms added to the Browser after the model's training cut-off
- notes: ['Draft models are untouched by any phenotype data, so this is a true prospective test of automated reconstruction.']

## Results

- condition_level: {"n_conditions_mapped": 25, "n_conditions_wt_grows": 0, "wt_growth_recall": 0.0, "conditions_with_absent_exchange": 9}
- gene_level_conditions_where_wt_grows: {"n_genes": 0, "n_conditions": 0, "n_gene_condition_pairs": 0}
- gene_level_all_mapped_conditions: {"n_genes": 469, "n_conditions": 25, "n_gene_condition_pairs": 11725, "aucpr_bernstein": {"point": NaN, "ci95": [NaN, NaN]}, "aucpr_standard": {"point": 0.10311300639658849, "ci95": [0.08108315565031983, 0.12649040511727078]}, "auroc_standard": {"point": 0.5, "ci95": [0.5, 0.5]}, "mcc": {"point": 0.0, "ci95": [0.0, 0.0]}, "balanced_accuracy": {"point": 0.5, "ci95": [0.5, 0.5]}, "accuracy": {"point …
- gene_map: {"model_genes": 675, "mapped": 674, "matched_by_version": 674, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 0, "browser_genes_hit": 674, "mapped_with_fitness_data": 507}
- counts: {"model_genes": 675, "model_genes_mapped": 674, "genes_with_fitness": 507, "genes_after_adjustment": 469, "conditions_total": 47, "conditions_mapped": 25, "conditions_wt_grows": 0}
- timings_s: {"rich_medium_essentials_s": 2.201889753341675, "knockout_simulation_s": 0.3440845012664795, "total_s": 2.9485666751861572}
- dropped_rich_medium_essentials: 38

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
