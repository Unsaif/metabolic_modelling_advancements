#!/usr/bin/env python3
"""Collect the independent check's numbers into independent_check_summary.json (claim, paper value, recomputed value).
Run after check_reference_comparison.py, check_coverage_and_inputs.py, check_timeline_and_locks.py and
check_ceiling_and_equivalence.py. No gembench or scripts/ imports."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
L = lambda n: json.load(open(os.path.join(HERE, n)))  # noqa: E731
R, C, T, E, W = (L("check_reference_comparison.json"), L("check_coverage_and_inputs.json"), L("check_timeline_and_locks.json"),
                 L("check_ceiling_and_equivalence.json"), L("web_checks.json"))
cm = R["check2_common"]
cov = C["coverage"]
r3 = lambda x: round(x, 3)  # noqa: E731

claims = [
    {"id": 1, "claim": "iML1515 grows in all 32 E. coli conditions; MCC 0.59",
     "paper": "32/32; 0.59", "recomputed": f"{R['check1_iML1515']['n_conditions_wt_grows']}/{R['check1_iML1515']['n_conditions_mapped']}; "
                                         f"{R['check1_iML1515']['mcc_le']:.4f}", "verdict": "confirmed"},
    {"id": 2, "claim": "common genes 764, 28 conditions; B0 0.53, U' 0.57, M 0.61, iJN1463 0.56",
     "paper": "764/28; 0.53/0.57/0.61/0.56",
     "recomputed": f"{cm['B0']['n_common_genes']} (M: {cm['M']['n_common_genes']})/{cm['B0']['n_conditions_both_grow']}; "
                   f"{r3(cm['B0']['mcc_arm_le'])}/{r3(cm[chr(85)+chr(39)]['mcc_arm_le'])}/{r3(cm['M']['mcc_arm_le'])}/{r3(cm['B0']['mcc_curated_le'])}",
     "verdict": "confirmed (M's common set has 765 genes)"},
    {"id": "2b", "claim": "curated - B0 +0.035 [-0.063, 0.140]; curated - M -0.042 [-0.139, 0.054]",
     "paper": "+0.035 [-0.063, 0.140]; -0.042 [-0.139, 0.054]",
     "recomputed": {"B0": [r3(cm["B0"]["curated_minus_arm"]), [r3(x) for x in cm["B0"]["ci95"]]],
                    "U'": [r3(cm["U'"]["curated_minus_arm"]), [r3(x) for x in cm["U'"]["ci95"]]],
                    "M": [r3(cm["M"]["curated_minus_arm"]), [r3(x) for x in cm["M"]["ci95"]]],
                    "ci95_by_seed": {k: v["ci95_by_seed"] for k, v in cm.items()}},
     "verdict": "confirmed in substance (intervals differ by seed only)"},
    {"id": 3, "claim": "union: iJN1463 0.45 against 0.49 to 0.58", "paper": "0.45; 0.49-0.58",
     "recomputed": {k: r3(v) for k, v in R["check3_union"]["four_way"]["mcc"].items()},
     "pairwise_union_curated_minus_arm": {k: [r3(v["curated_minus_arm"]), [r3(x) for x in v["ci95"]]] for k, v in R["check3_union"]["pairwise"].items()},
     "verdict": "confirmed"},
    {"id": 4, "claim": "iJN1463 grows in 36 of 43, drafts 28; 28 subset of 36",
     "recomputed": {"iJN1463": R["check4_coverage"]["per_model"]["iJN1463"]["grows"], "drafts": R["check4_coverage"]["per_model"]["B0"]["grows"],
                    "subset": R["check4_coverage"]["draft_28_subset_of_curated_36"],
                    "screened_carbon_source_conditions": R["check4_coverage"]["conditions_total_in_card_counts"]["iJN1463"]},
     "verdict": "numbers confirmed; '43 conditions in which P. putida was screened' is wrong (57 screened, 43 mapped to BiGG)"},
    {"id": 5, "claim": "12 organisms: 33%-75%, median 60%, 150 of 267; U' 155 of 267",
     "recomputed": {k: cov[k] for k in ("B0_fraction_min", "B0_fraction_max", "B0_fraction_median", "B0_grows_total", "U_grows_total",
                                         "mapped_total", "conditions_total_screened_in_mapped_media", "B0_share_of_all_screened_in_mapped_media")},
     "verdict": "numbers confirmed; denominator is mapped conditions, not screened (>=321 screened)"},
    {"id": 6, "claim": "same protocol, no gap-fill, identity mapping, iML1515 sha256 and source",
     "recomputed": {k: {x: C["reference_inputs"][k][x] for x in ("params_equal_to_putida_arms", "media_table_equal", "carbon_table_equal",
                                                                    "applied_transforms", "identity_mapping_in_matrices", "sha256_matches",
                                                                    "fitness_rebuilt_from_raw", "software")} for k in ("iJN1463", "iML1515")},
     "putida_arms_software": C["reference_inputs"]["software_putida_arms"],
     "iML1515_github_sha256_equal": W["iML1515_source"]["equals_local_file_and_card"],
     "verdict": "confirmed; environment differs (Python 3.13.16 vs 3.11.15) and BW25113 deletions not applied to iML1515"},
    {"id": 7, "claim": "94% = (273+94)/(274+117)", "recomputed": C["changed_predictions"]["recomputed"],
     "share": C["changed_predictions"]["recomputed_share"], "verdict": "confirmed (recomputed from matrices)"},
    {"id": 8, "claim": "Supplementary Table S1", "recomputed": {k: {x: v[x] for x in ("created_at", "commit_time", "external_minus_commit_s")}
                                                               for k, v in T["manifests"].items()},
     "downloads": {k: [v["first_start"], v["last_finish"]] for k, v in T["downloads"].items()},
     "verdict": "confirmed except 13:07:20 (record 13:07:19.691; other times truncated)"},
    {"id": 9, "claim": "only P. putida and M. tuberculosis have hand-curated models in BiGG", "recomputed": W["bigg_models_list"],
     "verdict": "confirmed for BiGG; abstract/Limits wording broader than BiGG is wrong (iAH991, iSO783, iGD1575)"},
    {"id": 10, "claim": "references 17 and 18", "recomputed": [W["reference_17_crossref"], W["reference_18_crossref"]],
     "verdict": "confirmed against Crossref; ref 18 names the model iJN1462"},
]
context = {
    "own_scope_mcc_putida": {k: r3(v["mcc_le"]) for k, v in R["own_scope_mcc_putida"].items()},
    "error_overlap_common": R["error_overlap_common"],
    "union_decomposition_Uprime_vs_iJN1463": R["union_decomposition_Uprime_vs_iJN1463"],
    "curated_by_condition_set": R["curated_by_condition_set"],
    "replicate_agreement": {k: {x: E[k][x] for x in ("mcc_all_replicate_pairs", "by_set")} for k in E if k.startswith("replicate")},
    "equivalence": E["equivalence"],
    "late_H3_arms": "MR1 12:28:54, Cola 12:29:12, RPal_CGA009 12:29:25 (M_noIonR6), after panel A analysis, before the 12:31:09 plan lock",
}
json.dump({"claims": claims, "context": context}, open(os.path.join(HERE, "independent_check_summary.json"), "w"), indent=1, default=float)
print("written", os.path.join(HERE, "independent_check_summary.json"))
