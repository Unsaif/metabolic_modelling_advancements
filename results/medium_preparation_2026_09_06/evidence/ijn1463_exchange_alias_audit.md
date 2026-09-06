# iJN1463 exchange identifier exception

6 September 2026. Additive source audit prompted by the preparation guard's refusal of `EX_AEP_e`. No optimization, phenotype analysis or model amendment was performed. The original inventory and primary attempt remain unchanged.

An exhaustive check of all **348 `EX_` reactions** in the exact local iJN1463 finds **one** reaction-name/metabolite-name mismatch:

| Field | Stored value |
|---|---|
| Reaction | `EX_AEP_e`, named “AEP exchange” |
| Single metabolite | `2ameph_e`, named “2-Aminoethylphosphonate” |
| Coefficient | −1 |
| Compartment | `e` |
| Formula / charge | `C2H7NO3P` / −1 |
| Original bounds | [0, 999999] |
| Gene rule | Empty |
| Reaction identifiers | BiGG `EX_AEP_e`; MetaNetX `MNXR95499`; SBO exchange `0000627` |
| Metabolite identifiers | BiGG `2ameph`; KEGG `C03557`; MetaNetX `MNXM1692` |

Every exchange has a single metabolite with coefficient −1, including this exception. There is no `EX_2ameph_e` reaction, no `AEP_e` metabolite and no second exchange on `2ameph_e`. Thus the exception concerns the correspondence between existing identifiers; the source does not contain two competing exchanges or an ambiguous mixture of metabolites.

Neither identifier is requested among the curated control's 43 mapped carbon conditions. Neither `AEP` nor `2ameph` is a component of its declared MOPS or RCH2 base-medium mappings. A narrowly declared `EX_AEP_e → 2ameph_e` identity adapter can therefore preserve reaction identity and chemistry while allowing the same exchange to be recognized and reset. It need not add a nutrient, rename a reaction, add a transporter or change a biological assumption. This local source correspondence is not independent evidence of physiological transport.

The earlier inventory checked exchange **shape**, not whether the reaction stem equals the metabolite ID. Its claim that no nonunit or multi-metabolite exchange was found remains correct. The preparation guard detected the additional representational requirement before any curated-model optimization. Any adapter must be explicitly declared, tested and frozen in a fresh attempt; the partial primary results must be retained separately.

The model SHA256 is `d573833328ffae0dfa752a1fa3262ed939ed5862288beab287fca30d0fefb4a1`. The [reproducible audit](ijn1463_exchange_alias_audit.py) and [JSON evidence](ijn1463_exchange_alias_audit.json) preserve every checked pair, source annotations, medium-membership checks and four input hashes. The source is the already downloaded [BiGG iJN1463](https://bigg.ucsd.edu/models/iJN1463); existing model-source notices apply.
