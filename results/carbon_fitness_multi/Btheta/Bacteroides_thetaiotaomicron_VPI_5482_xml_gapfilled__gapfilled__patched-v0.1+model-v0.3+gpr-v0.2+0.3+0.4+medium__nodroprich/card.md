# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Btheta

Created: 2026-09-06T09:44:44Z

## Model

- model_id: Bacteroides_thetaiotaomicron_VPI_5482_xml_gapfilled
- file: models/gapfilled/Btheta.xml.gz
- source: EMBL GEMs (CarveMe draft, Machado et al. 2018) minimally gap-filled by gembench.gapfill for growth on Varel_Bryant_medium + EX_glc__D_e: added ['ASPO2y', 'DHORD6', 'FE3t', 'OCBT_2', 'PDX5PO2']
- version_note: 
- n_reactions: 1498
- n_metabolites: 1068
- n_genes: 686
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
- patches: {"file": "data/reference/universe_patches_v0.1.json", "model_file": "data/reference/model_patches_v0.3.json", "gpr_file": "data/reference/gpr_patches_v0.2.json,data/reference/gpr_patches_v0.3.json,data/reference/gpr_patches_v0.4.json", "applied": [{"reaction": "CBPS", "gpr_before": "NP_809470_1 or (NP_809469_1 and NP_809470_1)", "gpr_after": "NP_809469_1 and NP_809470_1"}, {"reaction": "PPCK", "bo …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/reference/fitness_browser_media_bigg.tsv
- carbon_source_mapping: data/reference/fitness_browser_carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept/Btheta_genpept_map.tsv"}
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

- condition_level: {"n_conditions_mapped": 25, "n_conditions_wt_grows": 16, "wt_growth_recall": 0.64, "conditions_with_absent_exchange": 8}
- gene_level_conditions_where_wt_grows: {"n_genes": 510, "n_conditions": 16, "n_gene_condition_pairs": 8160, "aucpr_bernstein": {"point": 0.5694308457034167, "ci95": [0.46527704767374484, 0.6878599107433634]}, "aucpr_standard": {"point": 0.4710885457578758, "ci95": [0.37995855643982096, 0.5830634935143262]}, "auroc_standard": {"point": 0.7875396302607278, "ci95": [0.7406040576413654, 0.8448573614330568]}, "mcc": {"point": 0.573097226655 …
- gene_level_all_mapped_conditions: {"n_genes": 510, "n_conditions": 25, "n_gene_condition_pairs": 12750, "aucpr_bernstein": {"point": 0.5641729613121677, "ci95": [0.5227525783641876, 0.6106862335059436]}, "aucpr_standard": {"point": 0.22849617620077173, "ci95": [0.18439589240204954, 0.28306670735236256]}, "auroc_standard": {"point": 0.6652552970541702, "ci95": [0.6341741894783957, 0.7017628471659146]}, "mcc": {"point": 0.2341181784 …
- gene_map: {"model_genes": 686, "mapped": 685, "matched_by_version": 674, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 0, "browser_genes_hit": 685, "mapped_with_fitness_data": 510}
- counts: {"model_genes": 686, "model_genes_mapped": 685, "genes_with_fitness": 510, "genes_after_adjustment": 510, "conditions_total": 47, "conditions_mapped": 25, "conditions_wt_grows": 16, "medium_completion_exchanges_added": 4}
- timings_s: {"rich_medium_essentials_s": 3.5762786865234375e-06, "knockout_simulation_s": 26.311267614364624, "total_s": 26.815990209579468}
- dropped_rich_medium_essentials: 0

## Warnings

- 22 of 47 conditions have no BiGG mapping
- medium 'Varel_Bryant_medium' components absent from the model: ['EX_hco3_e', 'EX_mobd_e', 'EX_sel_e', 'EX_slnt_e', 'EX_tungs_e', 'EX_pheme_e']
- medium completion added exchange+uptake for: ['EX_cys__L_e', 'EX_met__L_e', 'EX_na1_e', 'EX_ni2_e']

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
