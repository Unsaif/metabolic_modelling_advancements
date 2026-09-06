# Local model comparison: quinone pathway

This extraction compares `models/bigg/iJN1463.xml` with Putida prepared by `scripts.run_quinone_biomass_sensitivity.prepare`, before condition-specific media/completion. It reads gene metadata, RefSeq/GenPept identities, model equations, bounds and annotations. It does not read fitness values, optimize a model, perform gene deletions, or apply a candidate.

The [exact source transfer candidate](curated_transfer_candidate.json) has SHA-256 `bee13fe5ebcc6f7ba98bf93c2ab99cccbf978f7308ddba0c8361a99a7ffcdb87`. It preserves the source file's equations, bounds, GPRs and metabolite metadata, with explicit model-gene aliases and uncertainties. It is a provisional representation, not an accepted physiological repair.

## Contiguous core route

| Source reaction | Source gene rule | Prepared draft |
|---|---|---|
| CHRPL | PP_5317 | Present |
| HBZOPT | PP_5318 | Present |
| OPHBDC | PP_5213 or PP_0548 | Present; draft uses PP_5213 only |
| OPHHX | PP_5013 | Present; historical assignment conflicts with current UbiB regulatory-kinase annotation |
| OHPHM | PP_1765 | Absent |
| OMPHHX | PP_5199 | Absent |
| OMBZLM | PP_5011 | Absent |
| OMMBLHX | PP_0427 or PP_5197 | Absent; source OR is not independently established |
| DMQMT | PP_1765 | Absent |

The missing tail runs from existing `2ohph_c` to `q8h2_c` through four absent intermediates: `2omph_c`, `2ombzl_c`, `2ommbl_c`, and `2omhmbl_c`. The five source reactions sum to:

`2ohph_c + 3 amet_c + o2_c → q8h2_c + 3 ahcys_c + 3 h_c`.

Every transferred reaction is atom/charge balanced under the curated source metadata; none is an exchange, demand, sink or artificial source. All source forward upper bounds remain `999999`. These are source-model bounds, not measured enzyme capacities. “Minimal” here means the contiguous absent tail of this specified source pathway; no optimization established globally minimal repair or growth rescue.

Prepared-draft gene PP_5011 is already represented by `NP_747113_1` and is reused. The other four candidate loci are absent from its current gene set and retain their PP locus IDs, which the existing mapper understands. Their current RefSeq protein identities are recorded separately rather than creating duplicate aliases. No locus-to-protein or existing-model-ID ambiguity was detected in the mapped source route.

## Precursors and existing alternatives

The [32-row pathway map](quinone_pathway_map.json) includes chorismate production from PEP/E4P, the MEP isoprenoid route from G3P/pyruvate, IPP/DMAPP, DMATT/GRTT, OCTDPS, prenyl attachment and ring modification. It explicitly stops at these central-carbon interfaces; it is not a reconstruction of all central metabolism or cofactor synthesis.

`DHQTi_copy1/copy2`, `PSCVT_copy1/copy2` and `MECDPS_copy1/copy2` in the curated source have stoichiometric counterparts under unsuffixed draft identifiers. Their different IDs do not establish missing chemistry. Balanced `MEPCT` is present; the source's additional `MEPCT_1` lacks a proton and has H/charge residuals of +1. That duplicate is not included in the candidate.

The draft also has alternative 4-hydroxybenzoate interfaces: CHRPL, reversible uptake `UHBZ1t_pp`, aldehyde oxidation `VNDH_2`, and reversible `SUCBZT2`. SAM production (`METAT`) and S-adenosylhomocysteine handling (`AHCYSNS`, `AHCi`) are represented. Adjacency establishes modeled reactions and directions, not their capacity to supply flux under a particular medium.

## Unresolved upstream redox representation

The curated MEP reaction `MECDPDH5` uses reduced flavodoxin:

`2mecdp_c + 2 flxr_c + h_c → 2 flxso_c + h2mb4p_c + h2o_c`.

The draft lacks `flxr_c`, `flxso_c` and that exact reaction. It has two alternatives:

- `MECDPDH`: `2mecdp_c + h_c → h2mb4p_c + h2o_c`, no explicit reducing donor; draft GPR maps to PP_0853.
- `MECDPDH2`: `2mecdp_c + nadh_c → h2mb4p_c + h2o_c + nad_c`; GPR maps to PP_0543.

Both appear balanced using the draft's zero-charge metabolite fields. Rechecking their common metabolites with the curated charges gives **−2** for `MECDPDH` and **0** for `MECDPDH2`. The first reaction therefore masks a redox/charge inconsistency; the second's balance does not independently establish its enzyme assignment. Equations and the charge calculation are retained in [upstream_redox_uncertainty.json](upstream_redox_uncertainty.json). No metadata or reactions were changed. These inherited alternatives remain separate from the frozen five-step transfer.

## Limits of the copied representation

The ten draft reactions touching `q8_c`/`q8h2_c` conserve their combined pool; none synthesizes net quinone. The copied tail introduces a source-model synthesis route without an artificial input. Physical tests are needed to determine the resulting model's capabilities.

The source oxygenase equations retain half-O2 abstractions without explicit reducing-donor demand. Atom balance does not validate those enzyme-level costs. Q8 chain length, UbiD/UbiX cofactor dependence, the UbiB accessory role, and the substrate specificity of the terminal hydroxylase remain unresolved by copying a model. The separate [biological source audit](../../../docs/evidence/09-putida-quinone-repair-sources.md) distinguishes exact transfer from the sequence-informed UbiI/UbiH/Coq7 hypothesis.

Exact source/draft Boolean rules and direction comparisons are in [quinone_pathway_gpr_details.tsv](quinone_pathway_gpr_details.tsv). The `draft_gprs_browser` column in the original summary TSV lists participating genes, not Boolean structure; use the detailed table or JSON when interpreting AND/OR relationships. Full identities, CDS records and source hashes are preserved in the JSON outputs.
