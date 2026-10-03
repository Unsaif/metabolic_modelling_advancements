# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism MycoTube

Created: 2026-10-03T10:26:12Z

## Model

- model_id: MycoTube_transfer_v1_base
- file: models/transfer_v1/MycoTube_base.xml.gz
- source: transfer_v1 base (models/embl_pinned/MycoTube.xml.gz)
- version_note: 
- n_reactions: 1270
- n_metabolites: 937
- n_genes: 803
- sha256: 4bdd558f8e542aa476183236cd0d16e8c84a6cb1a30659dd0bd1813769f731f4

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'MycoTube'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 2881
- n_experiments: 212

## Protocol

- study: transfer_v1
- arm: M
- arm_definition: {"transforms": ["universal_reaction_patches", "normalize_conjunction_rules", "universal_model_additions", "add_atp_synthase_if_annotated", "remove_menaquinol_if_no_pathway", "decisions:data/studies/transfer_v1/decisions/{org}.json"]}
- applied: [{"transform": "universal_reaction_patches", "records": []}, {"transform": "normalize_conjunction_rules", "records": [{"reaction": "SSALy", "gpr_before": "NP_214737_1 or NP_214748_2 or NP_216247_2 or (NP_214748_2 and NP_216247_2)", "gpr_after": "NP_214737_1 or (NP_214748_2 and NP_216247_2)"}]}, {"transform": "universal_model_additions", "records": [{"reaction": "GCALDt", "added": "gcald_e <=> gcal …
- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- media_mapping: data/studies/transfer_v1/media_bigg.tsv
- carbon_source_mapping: data/studies/transfer_v1/carbon_sources_bigg.tsv
- gene_mapping: {"method": "RefSeq protein accession -> /locus_tag from NCBI GenPept records (efetch, db=protein, rettype=gp), 5 September 2026; locus tags normalised to Fitness Browser sysName by removing underscores where needed", "genpept_table": "data/genpept_panel/MycoTube_genpept_map.tsv"}
- role: evaluation_panel_A
- base_model: {"org": "MycoTube", "role": "evaluation_panel_A", "draft": "models/embl_pinned/MycoTube.xml.gz", "draft_sha256": "c0c8fb4a6ba0233d013f8e9be6adf54a7c41983a4d0881b068140364ebc6c236", "universe": "external/carveme/universe_bacteria.xml.gz", "universe_sha256": "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19", "reference": {"medium": "Sautons minimal media with no carbon", "carbon_exc …

## Leakage

- ground_truth_used_in_model_curation: no for the correction rules: frozen before this organism's fitness tables were downloaded (see the study freeze manifests); the base model's gap-fill uses observed wild-type growth on one reference condition
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: see docs/studies/transfer-method-v1.md and the panel selection record
- notes: ['Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately.']

## Results

- condition_level: {"n_conditions_mapped": 19, "n_conditions_wt_grows": 19, "conditions_with_absent_exchange": 10}
- gene_level_conditions_where_wt_grows: {"n_genes": 463, "n_conditions": 19, "n_gene_condition_pairs": 8797, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.04396491558495129, "ci95": [0.026105539588222076, 0.06452804072269061]}, "aucpr_standard": {"point": 0.005725528589071665, "ci95": [0.003674886997287462, 0.008951202512723825]}, "auroc_standard": {"point": 0.4689112965872016, "ci95": [0.3955573400 …
- gene_map: {"model_genes": 803, "mapped": 802, "matched_by_version": 794, "matched_by_accession_only": 0, "no_genpept_record": 1, "locus_tag_not_in_browser": 0, "browser_genes_hit": 802, "mapped_with_fitness_data": 463}
- counts: {"model_genes": 803, "model_genes_mapped": 802, "genes_with_fitness": 463, "genes_after_adjustment": 463, "conditions_total": 20, "conditions_mapped": 19, "conditions_wt_grows": 19, "medium_completion_exchanges_added": 2}
- timings_s: {"rich_medium_essentials_s": 1.6689300537109375e-06, "knockout_simulation_s": 19.803802967071533, "total_s": 20.073816537857056}

## Software

- python: 3.11.15
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.4.4
- scipy: 1.17.1
- scikit-learn: 1.8.0
