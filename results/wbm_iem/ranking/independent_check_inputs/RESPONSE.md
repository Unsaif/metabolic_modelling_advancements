# Main session's response to the inputs check

6 October 2026, before the matrix run.

| # | Problem | Response |
|---|---|---|
| 1 | HPO annotations lost through the v0.1 linking tables | **Fixed by design.** HPO profiles v0.2 are built directly from Orphanet's annotations of each linked disorder (plan, "HPO profiles v0.2"). |
| 2 | Wrong or missing Orphanet links (HLYS1, STAR, 2OAA, ADSL) | **Relinked in v0.2.** Every IEM's link is reviewed, including those pointing to entries without HPO annotations, and each change is recorded with its reason in `data/iem/iem_orphanet_links_v0.2.tsv`. Readouts the new profiles need are computed in a supplementary run with identical settings, before any ranking. |
| 3 | Corroboration by substring matching | **Redefined.** A lab tuple is corroborated when the same readout and direction appear in the disorder's v0.2 HPO profile. |
| 4 | Confusion analysis missing; candidate list cut to 10 | **Fixed** in `scripts/iem_disease_ranking.py`. The full lists of candidates scoring higher and tied are kept. The confusion count gives, per candidate, the profiles where it scores strictly higher and the ties with a positive true score. A test is added. |
| 5 | Parser misses directed terms | **Fixed by design.** v0.2 takes every annotation under HP:0001939, or under any other branch holding a v0.1 metabolite term, and maps or excludes each term explicitly. |
| 6a | 11-deoxycortisol label vs definition | **Label kept**, by the rule that the term's label decides when label and definition disagree. The disagreement is noted in the map. The readout DM_11docrtsl[bc] stays in the panel. |
| 6b | Folate | **Exclusion kept** (group of vitamers). The note now cites the definition's "folic acid". |
| 6c | Homocystine | **Kept**, by the label. |
| 6d | 3-methylglutaric acid | **Kept**, by the label. |
| 7 | "Inorganic ions" wording; L-cystine locations | **Reworded** to "electrolytes (Na⁺, K⁺, Ca²⁺, Mg²⁺, phosphate) and the anion gap" in the plan and the v0.2 map. The L-cystine note now names the gut lumen and faeces as well. |
| 8 | Builder robustness; output name | The v0.2 builder raises on anything it cannot parse and keys corroboration by IEM and call. The ranking script's default output name is fixed. |

## Sources downloaded by the checking agent

The checking agent downloaded these during its check, without asking first. Tim approved their use afterwards. They
are kept in the session's scratch space, outside the repository.

- **Orphanet en_product4.xml.**
  - URL: https://media.githubusercontent.com/media/Orphanet/Orphadata_aggregated/master/Rare%20diseases%20with%20associated%20phenotypes/en_product4.xml
  - 47,893,357 bytes; sha256 4f44e8a61201399911aa1ba44a293c0ccaa5ce11272c47d40862255cb72f6b32.
  - JDBOR date 2026-06-23; licence CC-BY-4.0.
- **HPO hp.obo.**
  - URL: https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo
  - 10,863,613 bytes; sha256 93dace952fcb3ec4728818857f6ba76bc2d5312f4d83266519b8694b8e798f22.
  - data-version hp/releases/2026-09-01.
