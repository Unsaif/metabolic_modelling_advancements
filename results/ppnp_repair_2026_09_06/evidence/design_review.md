# Source-only design review: a provisional PpnP route

6 September 2026. This review specifies a testable model intervention after the quinone experiment exposed the ribokinase dependency. It uses saved model chemistry, directions, gene identities and source evidence. **No intervention, optimization or numeric fitness access was performed.** The proposed PP_4248 assignment remains a homology-based hypothesis, and any later agreement with the exposed Putida benchmark will remain a development result.

## Exact chemistry and metadata

The [candidate JSON](reaction_candidates.json) supplies the runner interface, exact source records, target metadata, complete local adjacency and input hashes. Its baseline is the unchanged saved `curated_path_template/model.xml.gz`, with the existing five quinone steps and the template quinone demand retained.

The local iJN1463 model has no exact forward or reverse equivalent of either candidate reaction. The original CarveMe universe was absent from this checkout and the user-supplied original path. The [public upstream file](https://github.com/cdanielmachado/carveme/blob/master/carveme/data/generated/universe_bacteria.xml.gz) was recovered into [source/universe_bacteria.xml.gz](source/universe_bacteria.xml.gz); its 595,770 bytes have SHA256 `b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19`, exactly matching the previously recorded Putida gap-fill universe. This establishes the recovered file's identity without silently substituting a newer universe.

| Proposed/source ID | Exact phosphorolysis stoichiometry | Source bounds | Primary hypothesis bounds |
|---|---|---|---|
| `PUNP1` | `adn_c + pi_c → ade_c + r1p_c` | −1,000 to +1,000 | 0 to +1,000 |
| `PUNP5` | `ins_c + pi_c → hxan_c + r1p_c` | −1,000 to +1,000 | 0 to +1,000 |

The universe contains exact duplicate `PUNP1_1` and `PUNP5_1` records. The canonical IDs retain richer database annotations; only those two are proposed. Neither reaction exists in the saved baseline under any exact equivalent forward or reverse stoichiometry. Both source GPRs are empty. Assigning both new reactions to `PP_4248` is an explicit new hypothesis, not a copied source-model gene assignment.

All six involved metabolites already exist. Their native/source and target IDs and `C_c` compartment labels match exactly; no new metabolite, boundary reaction or compartment alias is needed.

| Metabolite | Formula | Source charge | Preserved target charge |
|---|---|---:|---:|
| Adenosine, `adn_c` | C10H13N5O4 | 0 | 0 |
| Inosine, `ins_c` | C10H12N4O5 | 0 | 0 |
| Adenine, `ade_c` | C5H5N5 | 0 | 0 |
| Hypoxanthine, `hxan_c` | C5H4N4O | 0 | 0 |
| Phosphate, `pi_c` | HO4P | −2 | 0 |
| Ribose 1-phosphate, `r1p_c` | C5H9O8P | −2 | 0 |

Both equations are atom and charge balanced under the source metadata. The target's two zero-charge placeholders cancel as well, but that cancellation does not validate their charges. Preserve the target metadata explicitly and reject any undeclared formula/compartment difference. Inherited names and annotations differ in places and are recorded separately; these are not reasons to overwrite shared metabolite objects.

## Direction and identity are separate evidence decisions

The primary hypothesis is **forward phosphorolysis only**. The original Sévin study supplement supports phosphate-dependent cleavage of both adenosine and inosine by E. coli YaiE/PpnP. Its Figure S8 and Table 9 conflict about reverse nucleoside synthesis. Preserve that inconsistency; a separately labeled reversible sensitivity may use the source's bounds, but source reversibility does not resolve the experimental conflict. Free ribose plus phosphate is not a supported PpnP substrate pair. See [evidence brief 10](../../../docs/evidence/10-putida-ribosyl-disposal-sources.md) and the separate primary-source review for the evidence boundary. [Sévin et al., original supplement](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.4103/MediaObjects/41592_2017_BFnmeth4103_MOESM62_ESM.pdf).

Local annotation resolves `PP_4248` to current locus `PP_RS22065`, RefSeq `WP_003254278.1`, on `NC_002947.4` positions 4,834,852–4,835,136, encoding 94 residues. Neither `PP_4248`, its RefSeq-style ID, nor any mapped alias for that locus exists in the saved model. The unambiguous proposed model ID is therefore `PP_4248`. Exact sequence equality to Q88F51 and biochemical specificity are addressed by the independent sequence/source audit; locus resolution alone proves neither.

Relevant existing aliases are `PP_2458 → NP_744606_1` (RBK), `PP_4976 → YP_009237194_1` (AHCi), `PP_1777 → NP_743933_1` (PPM), and `PP_0591 → NP_742754_1` (ADA). The extraction checks uniqueness. Existing `PNP` acts on nicotinamide riboside, `PYNP2r` on uridine, and `NP1` on nicotinate riboside synthesis; their names do not make either proposed purine phosphorolysis reaction redundant.

The downstream PPM assignment is a weak link. The saved draft uses `PP_1777`; local iJN1463 uses `PP_1777 or PP_5288`. The historical PP_1777 annotation is phosphomannomutase, and species-specific ribose-phosphate activity is not established by the current evidence. Retain the exact baseline edge and gene rule for this experiment; report the uncertainty rather than silently importing the curated OR rule.

## Controls and their logical interpretation

Use a fixed glucose condition and the same growth, oxygen, maintenance, medium-completion and quinone-demand assumptions as the saved baseline. Do not change a bound after seeing a result. Genotype controls should resolve model identities before changing gene rules or knocking out genes; absent PP_4248 in the baseline must be labeled unrepresented.

| Declared comparison or control | What it can establish without presuming an outcome |
|---|---|
| Baseline wild type and baseline ΔRBK | Reference growth and the already derived exact ribokinase dependence of the saved network. |
| Forward PpnP wild type and ΔRBK | Whether the proposed route permits a feasible alternative under the fixed condition. Source equations alone do not guarantee growth rescue. |
| Forward PpnP ΔPP_4248 | Both additions must be disabled; the remaining reaction constraints should match baseline exactly. This is an implementation/equivalence control. |
| Forward PpnP ΔRBK ΔPP_4248 | Must recover the original exact no-growth implication. A positive solution requires investigating the implementation or numerical tolerance. |
| Forward PpnP ΔRBK ΔAhcY | Must block positive exact growth: with AHCi absent, the sum `ahcys_c + rhcys_c + rib__D_c` has no exit when RBK is also absent. PpnP consumes neither SAH nor free ribose. |
| Forward PpnP ΔPPM and ΔRBK ΔPPM | Diagnostic tests, **not guaranteed blocking controls**: the baseline already has an alternative r1p-to-r5p route. Compare gene deletion with explicit PPM reaction closure. |
| Separate adenosine-only and inosine-only additions, with ΔRBK and ΔRBK ΔADA | Resolve which substrate branch supports any rescue. In the inosine-only arm, ADA connects AhcY-produced adenosine to the new reaction, giving a stronger exact blocking expectation for the double deletion. No such ADA requirement follows for direct adenosine phosphorolysis. |
| Matched primary arm with RHCYS set to zero | Tests dependence on the inherited gene-less, biologically questionable ribose-producing step. It is a declared structural restriction, not an accepted correction to LuxS chemistry. |
| Separate reversible-bound sensitivity | Tests the impact of the unresolved reverse direction. Its effects must not be pooled with the forward-only primary result. |

The strong double-deletion implications concern the exact stored steady-state equations and the declared additions, not an experimentally established KT2440 phenotype. Gene deletion is also distinct from chemical inhibition: PP_0591 has DADA as another assigned activity, PP_3254 has 5DOAN, and PP_1777 participates in a PGMT OR rule as well as PPM. The artifact records every associated reaction so those differences remain reviewable.

## Why PPM closure is not a decisive negative control

The saved r1p neighborhood has four existing reaction interfaces. Besides PPM, `NP1` consumes r1p forward; `PNP` and `PYNP2r` can consume it in their allowed reverse directions. Direction compatibility is not, by itself, a full-network feasibility result. There is nevertheless a complete symbolic PPM alternative:

```text
PNP reverse:  H+ + nicotinamide + r1p → phosphate + nicotinamide riboside
RNMK forward: ATP + nicotinamide riboside → ADP + H+ + NMN
NMNN forward: H2O + NMN → H+ + nicotinamide + r5p

Net: r1p + ATP + H2O → r5p + ADP + phosphate + H+
```

All three directions are permitted. PNP and NMNN are gene-less; RNMK is assigned to PP_4218. This ATP-consuming bypass is verified by exact stoichiometric addition, without optimization. Therefore an RBK/PPM double mutant that still grows after PpnP addition would not by itself falsify the proposed ribosyl-disposal mechanism or demonstrate another PPM isozyme.

If an exact reaction-junction negative control is needed, predeclare **ΔRBK plus closure of every r1p-consuming direction in the forward-only candidate**: PPM both directions, NP1 forward, and the reverse directions of PNP and PYNP2r. Then the r1p balance forces both new phosphorolysis fluxes to zero. This is a deliberately constrained model diagnostic; it is not a biological genotype or a proposed physiological repair. Closing only the three-step nicotinamide bypass would still leave other nominal r1p outlets and would not establish exclusivity.

## Phosphate specificity and validation boundaries

Closing external phosphate is a poor whole-cell specificity control: phosphate is also required for biomass and many other reactions. A resulting loss of growth would not establish the PpnP mechanism. Keep the phosphate coefficient at −1 and verify exact atom balance. If a computational reaction-level specificity check is required, use a **separate, explicitly isolated substrate module** or biochemical assay comparing nucleoside plus phosphate with nucleoside alone, with measured base and r1p products. Such diagnostic substrate feeds do not belong in the candidate whole-cell model or benchmark. Never make the reaction phosphate-independent by deleting its phosphate term.

Before interpreting a later rescue, require consistent solver status, finite fluxes, mass/charge and bound checks, energy-from-nothing checks, a fixed positive quinone demand, and flux through the declared phosphorolysis route. Record whether r1p leaves through PPM or another represented route. A successful flux calculation would establish a property of this provisional network, not PP_4248 activity in KT2440 or independent predictive validation.

Reproduce the source extraction with `scripts/map_ppnp_candidates.py --out <new-path>`. The script rejects overwriting evidence, checks the source hash and all input hashes, and never calls an optimizer or reads numeric fitness tables.
