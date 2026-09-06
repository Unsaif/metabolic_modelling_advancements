# Development model compartment and boundary inventory

6 September 2026. Source structure only: no optimization, knockout simulation, energy test, numerical fitness parsing or score calculation. The CLI disables both COBRA model optimization methods as an additional guard. All work is additive; historical sources and results remain unchanged.

The [inventory JSON](medium_structure_inventory.json) records complete single-metabolite boundary inventories, inferred boundary classes, compartment assignments, medium requirements, historical reconstruction changes and exact input hashes. The [source script](inventory_medium_structure.py) exposes `build_parent(label)` and `metadata_conditions(org)` for a separately frozen comparison runner. The former returns a model before medium completion and its provenance; the latter reads experiment metadata and only the fitness table's first-line experiment names. File bytes are hashed without interpreting fitness values.

## Parent identity and availability

Putida is the exact serialized provisional PpnP `both_forward` model, SHA256 `b053283adf1257675b771439e99bb0dcaf2ee044c89e8df34096888c26d7f812`. Btheta, MR1 and Smeli have no serialized cycle-7 model in the development sprint output. They are explicitly **reconstructed in memory**, following their recorded cycle-7 invocations: saved gapfilled source, universe patch v0.1, model patch v0.4, and GPR patches v0.2/v0.3/v0.4 in order. Gene mappings use local gene and GenPept metadata. These reconstructed parents are not presented as newly recovered historical model bytes.

The already exposed Putida iJN1463 is an unmodified local structural reference using canonical compartments. Keio/iML1515 was unavailable in the clone and three expected original model locations; it is recorded as unavailable, with no new download or substitute organism.

## Observed representation differences

| Parent | Mapped conditions | Main compartments | EX IDs omitted from inferred exchanges, before → after completion | Completion exchanges added |
|---|---:|---|---:|---:|
| Putida PpnP | 43 | C_c / C_p / C_e, plus c/e | 2 → 7 | 7 |
| Btheta cycle 7 | 25 | C_c / C_p / C_e, plus e | 2 → 4 | 4 |
| MR1 cycle 7 | 12 | C_c / C_p / C_e, plus e | 1 → 6 | 7 |
| Smeli cycle 7 | 33 | C_c / C_p / C_e, plus c/e | 3 → 7 | 6 |
| Putida iJN1463 | 43 | c / p / e | 0 → 0 | 4 |

COBRA 0.32.1 infers `C_e` for all four draft-derived parents and `e` for iJN1463. The drafts already contain external-suffix metabolites in the second label `e`: Putida's glycolaldehyde and hexanoate; Btheta's fructose and glycolaldehyde; MR1's glycolaldehyde; Smeli's N-acetylglucosamine, glycolaldehyde and thiamine. Their `EX_` reactions are outside the automatically inferred exchange collection. Additional completion can create more such cases. This is an observed classification difference, not proof that the compounds occupy biologically different extracellular spaces.

Every `EX_` reaction found has exactly one metabolite with coefficient −1. No reversed-sign, multi-metabolite or nonunit exchange shape occurs in these five parents. The curated reference has 31 `DM_` and two sink-class boundaries. The drafts have respectively 7/1/9/7 sink-prefixed single-metabolite boundaries, but COBRA labels none as demand or sink: its name exclusions and reversibility rules leave these outside those inferred collections. An explicit medium helper must distinguish external exchange closure from preservation of internal demand/sink boundaries without treating the automatic collections as a complete inventory.

In every draft, sodium and nickel extracellular metabolites already exist in `C_e` even though their exchanges are missing. Inherited completion reuses those metabolites, while other missing external metabolites are created in `e`. Completion adds seven exchanges for Putida, four for Btheta, seven for MR1, six for Smeli and four for iJN1463, with a corresponding inward-only `MEDt_` carrier for each. The JSON lists every component and all policy exclusions or unavailable cytoplasmic targets.

**No known formula or charge mismatch was found between existing external and cytoplasmic metabolites for the requested medium pairs.** This negative finding is limited to those pairs and the metadata present. It is not a global chemical-consistency certificate. Existing metabolite identity, carrier collisions and declared compartment equivalence still require explicit validation.

## Proposed comparison panel

Use all **156 listed mapped conditions** across these five parents, retaining the historical union of media components for each organism. This spans Varel–Bryant, ShewMM, RCH2 and MOPS media, both compartment naming conventions, previously existing external metabolites, mixed labels, internal boundaries and different completion inventories. These remain exposed development models, not independent validation organisms.

Before real-model solves, freeze all parent/reconstruction inputs, experiment metadata, old and new preparation code, compartment declarations and completion exclusions. For every condition, compare prepared stoichiometry, GPRs, formula, charge, bounds and objective exactly. Any permitted compartment-label or display-name normalization must be declared separately. Numerical WT comparisons with GLPK and HiGHS can follow that structural gate; no mutant fitness rescore is needed for a preparation refactor.

Add synthetic tests for structures absent from the real inventory: positive-coefficient source boundaries; malformed multi-metabolite exchanges; ambiguous external-compartment ties; conflicting existing metabolite chemistry; carrier-ID collisions; and internal sink/demand preservation. Sequential medium reuse, fresh-copy equivalence, idempotence, atomic failure and input immutability should be tested explicitly. These fixtures exercise software semantics, not new biological hypotheses.
