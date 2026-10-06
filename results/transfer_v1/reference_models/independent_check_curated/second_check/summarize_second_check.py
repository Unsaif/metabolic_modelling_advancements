"""Collect the second check's reproduced numbers next to the reported ones (comparison.json, Paper 1 text).
Usage: python -I summarize_second_check.py   (writes summary_second_check.json)
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from common import ROOT  # noqa: E402


def j(name):
    return json.load(open(os.path.join(HERE, name)))


def main():
    comp = json.load(open(os.path.join(ROOT, "results/transfer_v1/reference_models/comparison.json")))["organisms"]
    e = comp["Btheta/iAH991_exploratory_B12_H2S"]
    b = j("recompute_btheta.json")["all"]
    rows = []

    def add(q, rep, mine):
        ok = abs(rep - mine) < 1e-9 if isinstance(rep, float) else rep == mine
        rows.append({"quantity": q, "reported": rep, "reproduced": mine, "match": ok})

    add("iAH991 exploratory own MCC", e["own_mcc"]["iAH991_exploratory_B12_H2S"], b["own_iAH991"]["mcc"])
    for k, kk in (("iAH991_exploratory_B12_H2S", "iAH991"), ("B0", "B0"), ("U'", "U'"), ("M", "M")):
        add(f"coverage {kk}", e["coverage"][k]["n_conditions_grows"], b["coverage"][kk])
    add("four-way union genes", e["four_way"]["union_genes"], b["four_way"]["union_genes"])
    add("four-way conditions", e["four_way"]["conditions_all_grow"], b["four_way"]["conditions_all_grow"])
    add("four-way observations", e["four_way"]["n_observations"], b["four_way"]["n_observations"])
    for k, kk in (("B0", "B0"), ("U'", "U'"), ("M", "M"), ("iAH991_exploratory_B12_H2S", "iAH991")):
        add(f"four-way MCC {kk}", e["four_way"]["mcc"][k], b["four_way"]["mcc"][kk])
    for n in ("B0", "U'", "M"):
        rc, mc = e["pairwise_curated_minus_arm"][n]["common"], b["pairwise"][n]["common"]
        add(f"common genes {n}", rc["n_common_genes"], mc["n_genes"])
        add(f"common curated - {n}", rc["difference"]["point"], mc["curated_minus_arm"])
        add(f"common curated - {n} ci95", [round(x, 10) for x in rc["difference"]["ci95"]], [round(x, 10) for x in mc["ci95_by_seed"][0]])
        ru, mu = e["pairwise_curated_minus_arm"][n]["union"], b["pairwise"][n]["union"]
        add(f"union curated - {n}", ru["union_delta"], mu["curated_minus_arm"])
        add(f"union curated - {n} ci95", [round(x, 10) for x in ru["union_delta_ci95"]], [round(x, 10) for x in mu["ci95_by_seed"][0]])
    for k in ("n_genes", "n_conditions", "n_observations", "wrong_Uprime", "wrong_curated", "wrong_both"):
        add(f"error overlap {k}", e["error_overlap_common_genes_Uprime_vs_curated"][k], b["error_overlap"][k])
    add("error overlap expected", e["error_overlap_common_genes_Uprime_vs_curated"]["expected_both_if_independent"], b["error_overlap"]["expected_both_if_independent"])
    for k in ("n_genes", "n_conditions", "important_calls", "important_calls_confirmed"):
        add(f"curated-only {k}", e["curated_only_genes"][k], b["curated_only"][k])
    out = {"btheta_exploratory": rows, "all_match": all(r["match"] for r in rows),
           "rebuild_vs_pdf": {k: j("compare_rebuild_to_pdf.json")[k] for k in ("s10a_rows", "s10a_rows_compared", "problem_counts")},
           "poppler_spot_check": j("spot_check_poppler.json")["summary"],
           "rebuild_validation": j("validate_rebuild.json"),
           "iah991_growth": {k: v for k, v in j("check_iah991_growth.json").items() if k != "conditions"},
           "iah991_rescore": j("rescore_iah991.json"),
           "drafts_supplemented": {a: j(f"rescore_drafts_supplemented_{a}.json") for a in ("B0", "UNQ", "M")},
           "call_agreement": j("call_agreement.json")}
    json.dump(out, open(os.path.join(HERE, "summary_second_check.json"), "w"), indent=1, default=float)
    for r in rows:
        print(f"{'OK  ' if r['match'] else 'DIFF'} {r['quantity']:<36} {str(r['reported'])[:42]:<42} {str(r['reproduced'])[:42]}")
    print("all match:", out["all_match"])


if __name__ == "__main__":
    main()
