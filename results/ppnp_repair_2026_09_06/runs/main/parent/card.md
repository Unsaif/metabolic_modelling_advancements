# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Putida

Created: 2026-09-06T12:46:34Z

## Model

- model_id: Pseudomonas_putida_KT2440_xml_gapfilled__curated_path_template__ppnp_parent
- file: results/ppnp_repair_2026_09_06/runs/main/parent/model.xml.gz
- source: Exact saved provisional quinone-demand parent plus declared PpnP/direction/RHCYS intervention
- version_note:
- n_reactions: 1912
- n_metabolites: 1325
- n_genes: 1305
- sha256: 4f51b83e9812ab9060714ed85d7ead79a0fd1496c856c5c7d10a32bcfc389ed7

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Putida'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 4778
- n_experiments: 314

## Protocol

- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- evaluation_role: development
- intervention: {"configuration": {"id": "parent", "substrates": [], "close_rhcys": false, "direction": "forward", "run_benchmark": true}, "source_reactions": [], "added_reaction_records": [], "closed_reactions": [], "gene_aliases": {"PP_4248": "PP_4248"}, "evidence_grade": "KT2440 function by homology; E. coli substrate assays", "target_charge_note": "Preserved historical zero placeholders; source Pi and r1p cha …
- study_fingerprint: ec7dc1dfdb7eed71276dd3c84d33984340c44d68203bf18b5aeb46efecdb2ede

## Leakage

- ground_truth_used_in_model_curation: yes: these patch sets were developed and accepted/held/rejected using Fitness Browser phenotypes also used for scoring; annotation support does not remove this reuse
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: Freeze code, media rules and patches before selecting independent data. Use organisms/experiments not inspected during development, or nested evaluation that repeats every data-guided selection step. A retrospective split of already-inspected data does not restore independence.
- notes: ['Retrospective development evaluation. No independent held-out or prospective test is established.', 'Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare gene-level scores on matched genes, conditions and finite observations; report wild-type growth coverage separately.', 'Medium completion adds assumed gene-less transport; recipe presence alone does not establish uptake or its energy cost.']

## Results

- condition_level: {"n_conditions_mapped": 43, "n_conditions_wt_grows": 34}
- gene_level_conditions_where_wt_grows: {"n_genes": 1047, "n_conditions": 34, "n_gene_condition_pairs": 35598, "n_missing_fitness": 0, "n_nonfinite_simulation": 0, "aucpr_bernstein": {"point": 0.49939171817927247, "ci95": [0.37068562579727654, 0.6129579853699584]}, "aucpr_standard": {"point": 0.4015875293193779, "ci95": [0.30155842576315894, 0.5137033774925877]}, "auroc_standard": {"point": 0.7612094716058451, "ci95": [0.706722238547088 …
- counts: {"model_genes": 1305, "model_genes_mapped": 1304, "genes_with_fitness": 1047, "genes_after_adjustment": 1047, "conditions_total": 57, "conditions_mapped": 43, "conditions_wt_grows": 34, "medium_completion_exchanges_added": 7}

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
