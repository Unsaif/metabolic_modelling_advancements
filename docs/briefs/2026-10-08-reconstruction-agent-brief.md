# Brief: evaluating an agent that does metabolic reconstruction

8 October 2026. Written for a fresh Claude Code session by the Claude session that has worked on Tim's
"Metabolic modelling improvements" project since September. It is opinion and context, not instructions from Tim.
Tim's own messages decide what you do.

## The project

Tim is a PhD student in Ines Thiele's lab (constraint-based modelling, rare-disease diagnosis). The idea is to have
AI agents do metabolic reconstruction themselves: building, curating and testing models, with tools to check their
work, such as running the COBRA Toolbox in a MATLAB Docker image (`ce-certification-agent/chatimd_interface/
matlab-dockerfile/`). An earlier attempt, Astra, worked but was slow and expensive.

## Your first job: evaluate, don't change

1. Read the code and setup: the agent loop, prompts, tool calls, the Docker image and how MATLAB is started and
   called.
2. Find where time and money went in the Astra runs, from logs if they exist: language-model calls and tokens per
   step, MATLAB start-up and per-call overhead, solver time, retries.
3. Report what to keep, what to cut, and a leaner design (below), with a cost estimate per task.
4. Write a handoff, so later sessions start from your findings.

Change nothing until Tim has read the report.

## Python first

Most of what an agent needs runs in Python. MATLAB should be kept for the few Toolbox functions that have no
equivalent.

- **COBRApy** (0.32 or later) with **HiGHS** (`highspy`) for general models. Use `gurobipy` where speed matters:
  it works on Tim's Mac with his licence. The licence bundled with pip refuses whole-body models.
- **Whole-body models in Python already work.** The repository below has a Python port of the Toolbox's whole-body
  set-up and its IEM protocol:
  - set-up bounds identical to MATLAB's for all 83,521 Harvetta reactions;
  - the same calls as MATLAB on every biomarker where MATLAB returned an optimum.
- **If MATLAB is needed,** use Tim's queue runner (`~/Documents/matlab-agent-runner/queue/`, rules in its
  `QUEUE.md`) rather than an interactive or server session.
- **Reference code:** `github.com/Unsaif/metabolic_modelling_advancements`, branch `claude/opus-continuation`.
  - `gembench/`: `wbm.py`, `wbm_constraints.py`, `wbm_iem.py`, `wbm_iem_gurobi.py`, `fitness_browser.py`,
    `patches.py`, `gapfill.py`, `checks.py`.
  - `scripts/`: runners and analyses.
  - `docs/HANDOVER.md` and `docs/studies/`.

## Lessons from this project

**A model that solves can still be wrong.** This week's audit of Harvey 1.03d found:

- pathway maxima up to 1,400–2,100 times the dietary supply, which only cycles can produce;
- biomarkers capped by fixed limits in both states;
- disease knockouts that miss part of the enzyme (FED misses 7 LCAT reactions; HYPRO1 leaves 9 PROD2 reactions
  active);
- a headline accuracy of about 85% that equals what the constant guess "always increased" scores.

So build the tests that can fail before building the agent:

- mass and charge balance;
- energy from nothing: with all uptakes closed, no ATP (Fritzemeier et al. 2017);
- blocked and unbounded reactions;
- metabolic tasks;
- growth and gene-essentiality data for microbes (Fitness Browser pipeline in `gembench/`);
- for human models, the IEM disease ranking (`scripts/iem_disease_ranking.py`).

**Keep the language model out of inner loops.** Gap-filling, flux variability, consistency checks and solver runs
are deterministic code. The model decides, cites evidence and writes the change log. Batch and cache calls, use
smaller models for routine steps, and log tokens and wall time for every step, against a budget.

**What has worked:**

- Rule-based corrections plus AI curation of draft bacterial models. Small gains, replicated exactly on a second
  panel, and on shared genes the corrected drafts showed no clear difference from hand-curated models.
- One change at a time, each with a test.
- A plan fixed before results are seen, held-out data for the final test, and an independent agent checking the
  results.

**Numerical traps on whole-body models** (about 80,000 reactions):

- Each solve takes about 8 s with Gurobi barrier and 17 s with HiGHS interior point on 2 cores.
- Simplex warm starts rarely help.
- A flux threshold of 1e-6 is below solver noise (re-solves differ by up to 4e-5).
- A bound computed with one solver can be infeasible in another.
- Always check the solver status; an infeasible or failed solve is never zero flux.

## Rules

- Never read, copy or handle licence files (`gurobi.lic`, MATLAB licences) or other credentials.
- One Gurobi job at a time, and no MATLAB queue jobs while one runs.
- MATLAB queue: a unique filename per job, never reused; the JSON sidecar first, then `.m.tmp`, renamed to `.m`.
  Do not start the old TCP MATLAB server.
- If a job times out, read its records; do not retry in a loop.
- Original model files are read-only. Save changed models elsewhere, with a record of every change.
- Numerical feasibility is not scientific acceptance.
- No paid API or agent runs (for example new Astra experiments) without Tim's explicit go-ahead and a budget.

## A first benchmark to propose

A small task with ground truth, where cost per accepted change can be measured:

- **Microbial option:** fix a draft model's known failures against Fitness Browser data. The pipeline and four
  organisms already exist.
- **Whole-body option:** repair the defects listed above in a copy of Harvey. Each test should fail before the fix
  and pass after it, with the IEM ranking as the outcome.

Success means accepted, test-backed changes at a known cost per change, not a model that merely solves.
