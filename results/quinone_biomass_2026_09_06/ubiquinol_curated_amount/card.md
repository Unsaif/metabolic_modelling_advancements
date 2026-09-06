# Benchmark card — carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism Putida

Created: 2026-09-06T10:32:56Z

## Model

- model_id: Pseudomonas_putida_KT2440_xml_gapfilled__ubiquinol_curated_amount
- file: models/gapfilled/Putida.xml.gz
- source: EMBL GEMs (CarveMe draft, Machado et al. 2018) minimally gap-filled by gembench.gapfill for growth on MOPS minimal media_noCarbon + EX_glc__D_e: added ['OXCDC', 'PHPYROX']
- version_note: 
- n_reactions: 0
- n_metabolites: 0
- n_genes: 0
- sha256: 378f67c3c8d4642c135942917cd46f3b6109c53c16f7987c9751fb643ade5739

## Dataset provenance

- dataset: Fitness Browser RB-TnSeq gene fitness, orgId 'Putida'
- primary_source: Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov
- download: 5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)
- n_genes_with_fitness: 4778
- n_experiments: 314

## Protocol

- params: {"carbon_uptake": -10.0, "growth_threshold": 0.001, "fitness_threshold": -2.0, "drop_rich_medium_essentials": false, "rich_medium_uptake": -1000.0, "knockout_genes": [], "processes": 2, "solver": "glpk", "max_conditions": null, "complete_medium_transport": true, "medium_completion_exclude": ["pnto__R", "fol", "hco3"], "genes_subset": null}
- evaluation_role: retrospective_development
- intervention: {"id": "ubiquinol_curated_amount", "metabolite": "q8h2_c", "coefficient_source": "q8h2_c coefficient in models/bigg/iJN1463.xml biomass BIOMASS_KT2440_WT3", "applied_coefficient": -0.000223}
- input_sha256: {"data/studies/quinone_biomass_sensitivity_v1.json": "4d5e20ee4c5f71dd59dbfc4417d132ad65f3f8feb9949b5de1dbb9e83d422b34", "scripts/run_quinone_biomass_sensitivity.py": "8b67123ef3f31d15b50436c09170c48de2574288ccaad6d1a0e139174f95d33d", "models/gapfilled/Putida.xml.gz": "378f67c3c8d4642c135942917cd46f3b6109c53c16f7987c9751fb643ade5739", "models/bigg/iJN1463.xml": "d573833328ffae0dfa752a1fa3262ed939e …

## Leakage

- ground_truth_used_in_model_curation: yes: these patch sets were developed and accepted/held/rejected using Fitness Browser phenotypes also used for scoring; annotation support does not remove this reuse
- ground_truth_public_since: Fitness Browser releases include Price et al. 2018; exact dates differ by experiment
- frontier_model_training_exposure: unknown; public accessibility does not establish inclusion in a particular model's training data
- held_out_recommendation: Freeze code, media rules and patches before selecting independent data. Use organisms/experiments not inspected during development, or nested evaluation that repeats every data-guided selection step. A retrospective split of already-inspected data does not restore independence.
- notes: ['Retrospective development evaluation. No independent held-out or prospective test is established.', 'Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.', 'Compare gene-level scores on matched genes, conditions and finite observations; report wild-type growth coverage separately.', 'Medium completion adds assumed gene-less transport; recipe presence alone does not establish uptake or its energy cost.']

## Results

- condition_level: {"n_conditions_mapped": 43, "n_conditions_wt_grows": 0}
- gene_level_conditions_where_wt_grows: {"n_genes": 0, "n_conditions": 0, "n_gene_condition_pairs": 0, "n_missing_fitness": 0, "n_nonfinite_simulation": 0}

## Software

- python: 3.14.3
- cobra: 0.32.1
- optlang: 1.9.1
- numpy: 2.5.2
- scipy: 1.18.1
- scikit-learn: 1.9.0
