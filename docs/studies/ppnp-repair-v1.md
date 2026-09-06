# PpnP recycling hypothesis: frozen development experiment

The preceding quinone repair restored synthesis and retained a growth-associated
quinone requirement. It also created 34 ribokinase (PP_2458) disagreements. An
exact stoichiometric argument showed that this dependency arose because the
model could dispose of the ribosyl part of S-adenosylhomocysteine only through
free ribose and ribokinase. This study asks whether a source-supported route
through nucleoside phosphorolysis removes that structural bottleneck.

This is an adaptive development experiment. The hypothesis was chosen after
observing those previous errors. Freezing the current proposal makes subsequent
changes auditable; it does not create independent validation or preregistration.
The existing historical patches and every prior study artifact remain unchanged.

## Candidate and evidence

The exact saved `curated_path_template` quinone model is the parent. Its biomass,
including quinone demand, stays unchanged in every arm. CarveMe universe reactions
PUNP1 and PUNP5 provide balanced adenosine/inosine plus phosphate conversion to
adenine/hypoxanthine plus ribose-1-phosphate. No metabolite, artificial sink,
transport reaction or maintenance change is introduced by this intervention.

The proposed gene is PP_4248. Exact sequence identity links the local KT2440 locus
to the PpnP annotation, but the functional assignment is by homology. The direct
substrate assays concern E. coli. The original supplement supports forward
phosphorolysis but contradicts itself about reverse catalysis. Accordingly the
primary candidate permits only forward flux; a separate arm tests reversibility.
The bounds are modelling capacities, not measured enzyme kinetics. Existing zero
charges for phosphate and ribose-1-phosphate remain disclosed placeholders;
source formulas and proper source charges must balance before transfer.

The downstream PPM reaction has an uncertain PP_1777 substrate assignment. A
pre-existing alternative route (reverse PNP, RNMK, NMNN) can convert ribose-1-P to
ribose-5-P at an ATP cost. Therefore deleting PPM is a diagnostic, not a guaranteed
blocking control. The old gene-less RHCYS chemistry is questionable; its closure
is tested separately and does not assert absence of an enzyme in the organism.

## Declared arms and checks

The recipe is `data/studies/ppnp_repair_v1.json`. Seven arms receive physical checks:
the unchanged parent; adenosine phosphorolysis alone; inosine phosphorolysis alone;
both forward reactions; parent with RHCYS closed; both forward reactions with
RHCYS closed; and both reactions reversible.

Four arms receive complete fresh fitness benchmark runs: parent, both forward,
both forward with RHCYS closed, and both reversible. The other three are physical
diagnostics only. No assertion of benchmark equivalence is made for them. Report
new gene coverage separately and compare common finite gene-condition observations
where both wild types grow; use the inherited thresholds and a 500-resample
gene bootstrap with seed 0. A parent rerun also permits comparison with the previous
saved study, without representing it as fresh biological validation.

Physical checks include wild-type glucose growth, net quinone production, all
five standard energy-from-nothing probes, eight specified single-gene deletions,
and specified double deletions involving ribokinase. In particular, removal of
the newly added reactions together with ribokinase is a rollback control. Combined
ribokinase/PpnP and ribokinase/AhcY deletion test the proposed route where both
functions are represented. Missing or orphan genes remain unrepresented rather
than being treated as dispensable. Save selected flux witnesses and their mass
and bound residuals; individual optimal fluxes are not unique or measured.

Numerical failure is an error, never zero growth. Preserve partial failures and
freeze any reviewed retry in a new directory. Independently verify benchmark
artifacts and glucose mechanistic LPs with two solvers. Source snapshots/hashes
must identify each verifier used. No reactions, gene rules, biomass coefficients,
media or thresholds may be selected using these outcomes. No candidate is promoted
to the historical model series from this experiment alone.

## Limits on interpretation

An apparent rescue establishes that the proposed chemistry can remove a model
bottleneck, not that KT2440 uses it in vivo. The inherited quinone representation,
demand amount, reducing-donor abstractions, upstream assignments, unvalidated
downstream reactions and zero positive ATP-maintenance requirement limit all
predictions. A more accurate development score cannot settle those assumptions.
Fresh independent data and suitable biological experiments remain future work.

Evidence is stored under `results/ppnp_repair_2026_09_06/evidence/`, including the
sequence, original substrate-assay, reaction-source and design audits.
