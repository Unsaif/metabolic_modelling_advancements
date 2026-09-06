# Preserved incomplete attempt

The no-demand baseline completed. The first ubiquinol arm failed before simulation because COBRApy's `get_coefficient` raises `KeyError` for a metabolite absent from the biomass reaction. The guard was corrected to use the biomass metabolite dictionary with an explicit zero default. No hypothesis, coefficients, model patches, or protocol settings were changed in response to outcomes.

The original runner and manifest are preserved here. The rerun also excludes the unrelated, concurrently developed `gembench/study.py` utility from its input list. The completed study is in `../quinone_biomass_2026_09_06/`; its baseline is recomputed. This is a development study, with both historical outcomes and this initial baseline already exposed before the retry. Neither manifest constitutes public preregistration or external validation.
