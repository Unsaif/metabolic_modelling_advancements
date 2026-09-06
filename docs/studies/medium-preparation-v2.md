# Explicit medium preparation v2: one declared identifier exception

6 September 2026. This is an additive revision of the
[v1 definition](medium-preparation-v1.md); all v1 code and artifacts remain intact.

The v1 run completed 113 conditions for the four draft models, preserving their
226 solver vectors. It then correctly refused iJN1463's `EX_AEP_e`, whose sole
metabolite is `2ameph_e` with coefficient −1. The naming convention check was
stricter than this source model's identifiers. The partial run is retained under
`results/medium_preparation_2026_09_06/runs/main/`, with its failure and no completed
panel summary. No curated-reference solve occurred in that attempt.

An independent source audit finds exactly this one identity mismatch among 348
iJN1463 EX reactions. The metabolite is the represented extracellular
2-aminoethylphosphonate. Neither this exchange nor its canonical spelling is
requested by the 43 mapped carbon conditions or their base media. This is an
identifier-format exception, not evidence for supplying a new nutrient.

The v2 API adds an explicit `exchange_metabolites` policy mapping for **existing**
reactions only. In this study the sole exception, on iJN1463 only, is
`EX_AEP_e -> 2ameph_e`. The actual reaction must still contain exactly the declared
external metabolite at coefficient −1. The adapter changes no reaction/metabolite
identifier, stoichiometry, gene rule, chemical annotation or medium request.
Missing bindings, nonexistent reactions and incompatible shapes remain errors.

New module `gembench/medium_preparation_v2.py` and runner
`scripts/run_medium_preparation_v2.py` preserve the v1 algorithm with this explicit
identity check extension. The v2 recipe adds the exact exception, source audit,
tests and the failed v1 artifacts to its input chain. All 156 original conditions
and five withdrawal tests are rerun in a fresh directory. Exact old/new LP,
GPR and chemical equality remains mandatory immediately before each solve.
Solvers, tolerances, missing-component policy and interpretation are unchanged.

Commit and push this revision before rerunning. A source exception selected
after a format failure is a disclosed software-development revision; it is not
an independent evaluation or a biochemical correction selected from fitness.
