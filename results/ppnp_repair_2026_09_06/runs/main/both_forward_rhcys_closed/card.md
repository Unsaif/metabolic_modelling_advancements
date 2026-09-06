# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Putida

Created: 2026-09-06T12:52:45Z

## Model

- model_id: Pseudomonas_putida_KT2440_xml_gapfilled__curated_path_template__ppnp_both_forward_rhcys_closed
- file: results/ppnp_repair_2026_09_06/runs/main/both_forward_rhcys_closed/model.xml.gz
- source: Exact saved provisional quinone-demand parent plus declared PpnP/direction/RHCYS intervention
- version_note:
- n_reactions: 1914
- n_metabolites: 1325
- n_genes: 1306
- sha256: 8510a97c7c4dbab112c9f1249a35f4f1c4ea846e2c6a0394c3c820d51100a2fe

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Putida'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 4778
- n_experiments: 314

## Protocol

- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- evaluation_role: development
- intervention: {"configuration": {"id": "both_forward_rhcys_closed", "substrates": ["adenosine", "inosine"], "close_rhcys": true, "direction": "forward", "run_benchmark": true}, "source_reactions": ["PUNP1", "PUNP5"], "added_reaction_records": [{"id": "PUNP1", "name": "Purine-nucleoside phosphorylase (Adenosine)", "equation": "adn_c + pi_c --> ade_c + r1p_c", "metabolites": {"ade_c": 1.0, "adn_c": -1.0, "pi_c":  …
- study_fingerprint: ec7dc1dfdb7eed71276dd3c84d33984340c44d68203bf18b5aeb46efecdb2ede

## Leakage

- ground_truth_used_in_model_curation: yes: these patch sets were developed and accepted/held/rejected using Fitness Browser phenotypes also used for scoring; annotation support does not remove this reuse
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: Freeze code, media rules and patches before selecting independent data. Use organisms/experiments not inspected during development, or nested evaluation that repeats every data-guided selection step. A retrospective split of already-inspected data does not restore independence.
- notes: ['Retrospective development evaluation. No independent held-out or prospective test is established.', 'Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare gene-level scores on matched genes, conditions and finite observations; report wild-type growth coverage separately.', 'Medium completion adds assumed gene-less transport; recipe presence alone does not establish uptake or its energy cost.']

## Results

- condition_level: {"n_conditions_mapped": 43, "n_conditions_wt_grows": 34}
- gene_level_conditions_where_wt_grows: {"n_genes": 1047, "n_conditions": 34, "n_gene_condition_pairs": 35598, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.5081152768302855, "ci95": [0.3865200114755869, 0.6266372910169251]}, "aucpr_standard": {"point": 0.41107900535196257, "ci95": [0.3089711146953494, 0.5268684473600408]}, "auroc_standard": {"point": 0.7617802046752788, "ci95": [0.7075729503003891, …
- counts: {"model_genes": 1306, "model_genes_mapped": 1305, "genes_with_fitness": 1047, "genes_after_adjustment": 1047, "conditions_total": 57, "conditions_mapped": 43, "conditions_wt_grows": 34, "medium_completion_exchanges_added": 7}

## Warnings

- Hypothesis selected after inspecting previous development errors; all phenotype data remain exposed development data, not independent validation.
- PpnP substrate assays concern E. coli; KT2440 PP_4248 is an exact sequence match to a homology-only annotation, not direct catalytic validation.
- Primary PUNP1/PUNP5 directions permit phosphorolysis only; source reversibility is a separate sensitivity because the original supplement contradicts itself about reverse catalysis.
- Existing PPM PP_1777 substrate assignment is provisional; another pre-existing ATP-dependent salvage route may bypass it.
- RHCYS closure is a sensitivity to questionable historical chemistry, not proof that KT2440 lacks all relevant enzymatic activity.
- Parent quinone tail length, demand amount, hydroxylase reducing donors, upstream assumptions and zero ATP-maintenance lower bound are inherited unresolved limitations.
- Reaction capacities of 1000 are modelling conventions, not measured enzyme capacities.
- Missing exported fitness rows do not establish mutant-library absence, essentiality or dispensability.

## Software

- python: 3.14.3
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.5.2
- scipy: 1.18.1
- scikit-learn: 1.9.0
