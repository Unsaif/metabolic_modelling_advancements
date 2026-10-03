# Transfer study v1: replication on panel B (plan)

Declared 3 October 2026 by Claude (Opus 5.5), after the panel A results ([transfer-v1-results.md](transfer-v1-results.md)) and **before** any panel B input beyond gene and experiment metadata was prepared and before any panel B fitness table was downloaded. Tim confirmed on 3 October that no external custodian is available, so this replication is again self-custodied (section 6).

## 1. Purpose

1. **Replicate** the frozen v1 method ([transfer-method-v1.md](transfer-method-v1.md)) exactly on panel B: Caulo (*Caulobacter crescentus* NA1000), Cup4G11 (*Cupriavidus basilensis* 4G11), Marino (*Marinobacter adhaerens* HP15), Miya (*Desulfovibrio vulgaris* Miyazaki F), PV4 (*Shewanella loihica* PV-4) and SB2B (*Shewanella amazonensis* SB2B). These were assigned to panel B at random on 3 October and nothing about their fitness has been downloaded.
2. Give a **pooled estimate** over panels A and B with a combination rule fixed now.
3. Test one **new secondary hypothesis** that development and panel A suggested (H3, section 4).

## 2. What is unchanged

Everything in sections 3–6 of the v1 method: gene mapping by protein sequence; the FEBA media rule and its component classes; the carbon-source rule; the reference-condition rule; the minimal gap-fill with the energy-cycle gate and six cut rounds; arms B0 to M_R1R2; blind adjudication (procedure v1, prompt v1, one fresh subagent per organism, one attempt, packets from the U′ model); the scoring protocol; the union-gene primary metric; gene and organism bootstraps; the sign test; the reading rule; and the failure handling.

The known v1 limitations are **kept on purpose**, so that this is a replication and not a new method: organic components of base media are unlimited (for example the lactate implied by the name of Miya's MoLS4 medium, as opposed to its MoLS4_no_lactate medium); counter-ions of trace salts are unlimited ions; no light input is modelled; the gap-fill stops after six energy-cycle cuts. Their effects will be reported per organism.

Implementation: `scripts/transfer_replication.py` points the frozen runner and packet builder at `data/studies/transfer_v1_replication/`, which holds `arms.json` (the v1 arms unchanged plus M_noIonR6), `organisms_panel_B.json`, the panel B media and carbon-source tables and section 3 of the component dictionary. No frozen v1 file is modified. Panel B decisions are written to `data/studies/transfer_v1/decisions/<org>.json`, as panel A's were.

## 3. Primary comparisons

- **Panel B (the replication):** H1-B, U′ (UNQ) against B0; H2-B, M against U′. Per-organism union-gene paired MCC differences; unweighted mean over evaluable organisms with a 10,000-resample organism bootstrap (seed 0); sign counts and exact sign test. Reading as in v1: supported if the mean is positive and the interval excludes zero; not supported if the mean is zero or negative; inconclusive otherwise. **H1 counts as replicated if H1-B is supported.**
- **Pooled, panels A and B:** the same statistics over all evaluable organisms of both panels (at most eleven; PS stays not evaluable). Panel A's numbers are already known, so the pooled estimate is not blind; it is reported as the combined estimate, and the panel B estimate remains the clean replication.

## 4. Secondary hypothesis H3 (new)

**H3: R6 assignments to inorganic-ion transport reactions make curation worse.** In development (MR1, a copper ABC transporter) and on panel A (*R. palustris*, a zinc ABC transporter) such assignments made the transporter genes "important" in every condition, and none was. A plausible mechanism is that models require trace metals for biomass while cells have redundant uptake routes.

- **Arm M_noIonR6:** U′ plus the organism's blind decisions, with every R6 decision on an inorganic-ion transport reaction turned into an abstention (`transfer_replication.py filter`).
- **Definition:** a reaction's transported species are the metabolites that occur in two or more of its compartments, ignoring protons and water. It is an inorganic-ion transport reaction when that set is non-empty and every member is in: zn2, cu2, cu, cobalt2, mn2, fe2, fe3, ni2, mobd, mg2, ca2, k, na1, cl, so4, so3, pi, ppi, nh4, no3, no2, tsul, slnt, sel, tungs, cd2, hg2, pb, cro4, aso3, aso4, hco3, co2, h2s, n2, o2 (`is_inorganic_ion_transport`, with tests).
- **Comparisons on panel B:** M → M_noIonR6 (H3) and U′ → M_noIonR6, with the same statistics and reading. An organism where the filter drops nothing has a difference of exactly zero and counts as unchanged.
- **Development basis (retrospective, not evidence for H3):** in the development and panel A decisions the filter drops 3 of 43 R6 decisions — MR1 CUabcpp, Cola ZN2tpp and RPal ZNabc. Union-gene M → M_noIonR6: MR1 +0.017, Cola 0.000, RPal +0.017 (`results/transfer_v1/development/MR1/M_noIonR6`, `results/transfer_v1/evaluation_panel_A/{Cola,RPal_CGA009}/M_noIonR6`).

## 5. Order of work

1. This plan, the wrapper and the replication arms are frozen (`results/study_freezes/transfer_v1_replication_plan.json`) and recorded in the Claude project.
2. Panel B inputs: pinned EMBL drafts (GitHub, commit 260d0f1); Fitness Browser and NCBI protein sequences; sequence gene maps; media and carbon-source mapping by the frozen rules (new dictionary components and carbon-source names in separately marked files); reference conditions and gap-fill; packets; blind adjudication; filtered decisions. Then an inputs freeze, also recorded in the project.
3. First download of panel B fitness tables, with hashes and times recorded.
4. Every arm run once; panel B and pooled analyses; independent recomputation by a separate agent.

## 6. Limits

- Self-custodied: the same model family designs, prepares and evaluates; separation rests on recorded order (git commits, manifests, project documents with server times).
- Two panel B organisms are *Shewanella* species, relatives of the development organism MR1. They were assigned at random and stay in the panel; results are shown per organism.
- After this replication no organism in the intersection of the Fitness Browser and the EMBL GEMs collection remains untouched; a further version of the method needs new data (for example newer Fitness Browser organisms or another fitness resource).
