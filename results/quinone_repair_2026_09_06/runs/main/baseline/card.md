# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Putida

Created: 2026-09-06T11:33:05Z

## Model

- model_id: Pseudomonas_putida_KT2440_xml_gapfilled__baseline
- file: /Users/timhulshof/.codex/.chatgpt-projects/g-p-6a9d1896e28881919f06af51e1e3c459/metabolic_modelling_advancements/results/quinone_repair_2026_09_06/runs/main/baseline/model.xml.gz
- source: Frozen v0.4 prepared model plus explicitly declared provisional intervention
- version_note: 
- n_reactions: 1907
- n_metabolites: 1321
- n_genes: 1301
- sha256: c147610f9f2063af7517c5f25810beeec109effa6176b52ccde8c71f1d963b83

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Putida'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 4778
- n_experiments: 314

## Protocol

- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- evaluation_role: development
- intervention: {"configuration": {"id": "baseline", "add_pathway": false, "biomass_coefficient_source": "none", "gpr_overrides": {}, "run_benchmark": true}, "transfer": null, "gpr_changes": [], "applied_biomass_coefficient": 0.0}
- study_fingerprint: d4a72ef731a9fa0451e0b1c9a0463837872e07dc1007772ce7ed2d1b53820734

## Leakage

- ground_truth_used_in_model_curation: yes: these patch sets were developed and accepted/held/rejected using Fitness Browser phenotypes also used for scoring; annotation support does not remove this reuse
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: Freeze code, media rules and patches before selecting independent data. Use organisms/experiments not inspected during development, or nested evaluation that repeats every data-guided selection step. A retrospective split of already-inspected data does not restore independence.
- notes: ['Retrospective development evaluation. No independent held-out or prospective test is established.', 'Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare gene-level scores on matched genes, conditions and finite observations; report wild-type growth coverage separately.', 'Medium completion adds assumed gene-less transport; recipe presence alone does not establish uptake or its energy cost.']

## Results

- condition_level: {"n_conditions_mapped": 43, "n_conditions_wt_grows": 34}
- gene_level_conditions_where_wt_grows: {"n_genes": 1047, "n_conditions": 34, "n_gene_condition_pairs": 35598, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.5138583932136408, "ci95": [0.39082686445499215, 0.6304530465670097]}, "aucpr_standard": {"point": 0.42879974899088785, "ci95": [0.319304765160858, 0.5388529313154024]}, "auroc_standard": {"point": 0.7625839995107795, "ci95": [0.7080815514904726, …
- counts: {"model_genes": 1301, "model_genes_mapped": 1300, "genes_with_fitness": 1047, "genes_after_adjustment": 1047, "conditions_total": 57, "conditions_mapped": 43, "conditions_wt_grows": 34, "medium_completion_exchanges_added": 7}

## Warnings

- These are exposed development data; no independent validation set exists.
- Source-model reaction transfer retains its bounds, oxygenase abstractions, Q8 representation and unmodified upstream assumptions.
- The sequence-informed GPR arm encodes an inferred UbiD/UbiX cofactor dependency, not a heteromeric catalyst; UbiB accessory requirements are not represented.
- No exported fitness row does not establish mutant-library absence, essentiality or dispensability.
- Native KT2440 quinone tail length, pool size and reducing-donor specificity remain unresolved.

## Software

- python: 3.14.3
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.5.2
- scipy: 1.18.1
- scikit-learn: 1.9.0
