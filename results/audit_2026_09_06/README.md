# Audit evidence, 6 September 2026

Start with [the scientific review](../../docs/reviews/2026-09-06-scientific-audit.md). Historical source commit: `b1997d3`.

- `historical_card_audit.json`: recalculated point metrics and source matrix hashes for 49 saved runs. The historical arrays were not rerun by this calculation.
- `paired_development_comparisons.json`: common-gene/common-condition comparisons, shared observation masks, paired gene-bootstrap intervals. These are development-data comparisons.
- `historical_verification.txt`: independent MCC and Bernstein AUC-PR recalculation, plus sampled raw-fitness and condition checks.
- `carbon_fitness_multi/MR1/`: a complete new run of the cycle-6 MR1 arm, with input fingerprints. `MR1_cycle6_reproduction.json` compares it with the historical arm; `rerun_verification.txt` is its independent verification.
- `prepared_model_energy_checks.json`: five closed-boundary energy-dissipation tests for each organism with model patches v0.3/v0.4.
- `gene_mapping_duplicate_check.json`: no many-to-one mappings in these eight prepared models.
- `added_reaction_balance_exceptions.json`: elemental/charge checks and the inherited zero-charge annotation limitation. Charge flags do not establish wrong stoichiometry.
- `atp_survey_summary.json`: offline count of the archived 500-model screen, not a new download or full functional ATP-synthesis test.
- `iem_legacy_audit.json`: source-pinned IEM protocol discrepancies, counts and reasons legacy predictions do not support validation claims.
- `test_results.txt`: combined regression-suite output from the audit environment, recorded in `requirements-audit.txt`.

Recompute historical results and the survey count from the repository root:

```sh
.venv/bin/python scripts/audit_saved_results.py
.venv/bin/python scripts/verify_carbon_fitness_multi.py
.venv/bin/python scripts/survey_embl_atp_synthase.py --summarize results/embl_atp_synthase_survey_sample500.tsv
.venv/bin/python -m pytest -q
```

To independently rerun the MR1 arm, use a new result directory:

```sh
.venv/bin/python scripts/run_carbon_fitness_generic.py --orgs MR1 --variant gapfilled --no-drop-rich --processes 2 --complete-medium-transport --patch data/reference/universe_patches_v0.1.json --model-patch data/reference/model_patches_v0.3.json --gpr-patch data/reference/gpr_patches_v0.2.json,data/reference/gpr_patches_v0.3.json,data/reference/gpr_patches_v0.4.json --output-dir results/reproduction_mr1
.venv/bin/python scripts/verify_carbon_fitness_multi.py --results-dir results/reproduction_mr1
```

The whole-body IEM engine was tested on controlled numerical fixtures, not on a complete new Harvey/Harvetta run. No independent biological validation set was introduced during this review.
