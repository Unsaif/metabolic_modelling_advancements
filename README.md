# Metabolic modelling advancements

Working repository for the "Metabolic modelling improvements" project: bringing frontier-model
effort to constraint-based metabolic modelling, with verification first. Tim Hulshof (Thiele lab,
University of Galway) and Claude.

Start with **`docs/roadmap/00-landscape-and-roadmap.md`** — the standing charter: where the field
stands in 2026, the diagnosis, the workstreams, the decision log, and the current sprint. The six
`docs/evidence/` briefs hold the sourced findings behind it; `docs/sprints/` holds the results of
each sprint.

## Layout

```
docs/          roadmap, evidence briefs, sprint notes
gembench/      GEM-Bench: benchmark harness (media, datasets, metrics, benchmark cards, FROG fixtures,
               whole-body-model loading/solving, IEM biomarker protocol)
scripts/       runnable entry points (see gembench-README.md)
results/       benchmark cards and result tables produced so far
data/
  fitness_browser/   RB-TnSeq gene fitness for 9 organisms (Fitness Browser; see PROVENANCE.md)
  iem/               inborn-error-of-metabolism biomarker ground truth v0.1 (linked to Orphanet/OMIM/HPO/HMDB)
tools/         small helpers
external/      (git-ignored) clones of public repositories used as data/model sources — see gembench-README.md
```

## Reproducing

`gembench-README.md` lists the Python dependencies, the external repositories to clone into
`external/`, and the commands for each benchmark. Reference data that is too large or not
redistributable (HPO, Orphanet XML, Human-GEM, Recon3D, Harvey/Harvetta) is fetched by the
commands documented there rather than stored here.

## Status (5 September 2026)

- E. coli carbon-source fitness benchmark (Bernstein et al. 2023 protocol) reproduced with bootstrap CIs.
- yeast-GEM 9.1.1 essential-gene test reproduced exactly.
- Harvey/Harvetta 1.03d solve with HiGHS (interior point, ~20 s) and GLPK in Python, coupling constraints included, solutions certified.
- IEM biomarker protocol ported to Python; full 63-IEM run in progress (shipped bounds; physiological/diet constraints not yet ported).
- IEM ground truth v0.1: 279 tuples from the lab's benchmark, linked outward; HPO cross-check done.
- Fitness Browser data for 9 organisms downloaded (5 Sept 2026).
