"""Collect the reported numbers (comparison.json, card.json) next to this check's reproduction and the
corrected-gene-map variant. Usage: python -I summarize.py   (writes summary.json)
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import OUT, ROOT, dump  # noqa: E402


def j(name):
    return json.load(open(os.path.join(OUT, name)))


def main():
    rep = json.load(open(os.path.join(ROOT, "results/transfer_v1/reference_models/comparison.json")))["organisms"]["MR1/iSO783"]
    a, b = j("recompute_comparison.json"), j("recompute_comparison_corrected_gene_map.json")
    rows = []

    def add(name, reported, mine, corrected):
        rows.append({"quantity": name, "reported": reported, "reproduced": mine, "corrected_gene_map": corrected,
                     "match": (abs(reported - mine) < 1e-9) if isinstance(reported, float) else reported == mine})

    add("iSO783 own MCC", rep["own_mcc"]["iSO783"], a["own"]["iSO783"]["mcc"], b["own"]["iSO783"]["mcc"])
    add("iSO783 genes scored", 560, a["own"]["iSO783"]["n_genes"], b["own"]["iSO783"]["n_genes"])
    for k in ("B0", "U'", "M"):
        add(f"{k} own MCC", rep["own_mcc"][k], a["own"][k]["mcc"], b["own"][k]["mcc"])
    for k in ("iSO783", "B0", "U'", "M"):
        add(f"{k} coverage (grows / mapped)", f"{rep['coverage'][k]['n_conditions_grows']}/{rep['coverage'][k]['n_conditions_mapped']}",
            f"{a['coverage'][k]['grows']}/{a['coverage'][k]['mapped']}", f"{b['coverage'][k]['grows']}/{b['coverage'][k]['mapped']}")
    add("four-way union genes", rep["four_way"]["union_genes"], a["four_way"]["union_genes"], b["four_way"]["union_genes"])
    add("four-way conditions", rep["four_way"]["conditions_all_grow"], a["four_way"]["conditions_all_grow"], b["four_way"]["conditions_all_grow"])
    for k in ("B0", "U'", "M", "iSO783"):
        add(f"four-way MCC {k}", rep["four_way"]["mcc"][k], a["four_way"]["mcc"][k], b["four_way"]["mcc"][k])
    for k in ("U'", "M"):
        add(f"share of gap closed {k} (not in plan)", rep["four_way"]["share_of_gap_closed"][k],
            a["four_way"]["share_of_gap_closed"][k], b["four_way"]["share_of_gap_closed"][k])
    for k in ("B0", "U'", "M"):
        rc, ac, bc = rep["pairwise_curated_minus_arm"][k]["common"], a["pairwise_curated_minus_arm"][k]["common"], b["pairwise_curated_minus_arm"][k]["common"]
        add(f"common genes ({k})", rc["n_common_genes"], ac["n_genes"], bc["n_genes"])
        add(f"common conditions ({k})", rc["n_conditions_both_grow"], ac["n_conditions"], bc["n_conditions"])
        add(f"common MCC arm {k}", rc["mcc_arm"], ac["mcc_arm"], bc["mcc_arm"])
        add(f"common MCC iSO783 vs {k}", rc["mcc_curated"], ac["mcc_curated"], bc["mcc_curated"])
        add(f"common curated - {k}", rc["difference"]["point"], ac["curated_minus_arm"], bc["curated_minus_arm"])
        rows[-1]["reported_ci95"] = rc["difference"]["ci95"]
        rows[-1]["reproduced_ci95_seed0"] = ac["ci95_by_seed"][0]
        rows[-1]["corrected_ci95_seed0"] = bc["ci95_by_seed"][0]
        ru, au, bu = rep["pairwise_curated_minus_arm"][k]["union"], a["pairwise_curated_minus_arm"][k]["union"], b["pairwise_curated_minus_arm"][k]["union"]
        add(f"union curated - {k}", ru["union_delta"], au["curated_minus_arm"], bu["curated_minus_arm"])
        rows[-1]["reported_ci95"] = ru["union_delta_ci95"]
        rows[-1]["reproduced_ci95_seed0"] = au["ci95_by_seed"][0]
        rows[-1]["corrected_ci95_seed0"] = bu["ci95_by_seed"][0]
    ro, ao, bo = rep["error_overlap_common_genes_Uprime_vs_curated"], a["error_overlap_common_genes_Uprime_vs_curated"], b["error_overlap_common_genes_Uprime_vs_curated"]
    for k in ("n_genes", "n_conditions", "n_observations", "wrong_Uprime", "wrong_curated", "wrong_both", "expected_both_if_independent"):
        add(f"error overlap {k}", ro[k], ao[k], bo[k])
    rc, ac, bc = rep["curated_only_genes"], a["curated_only_genes"], b["curated_only_genes"]
    for k in ("n_genes", "n_conditions", "important_calls", "important_calls_confirmed"):
        add(f"curated-only {k}", rc[k], ac[k], bc[k])
    out = {"rows": rows, "all_reported_numbers_reproduced": all(r["match"] for r in rows),
           "seed_sensitivity": j("seed_sensitivity.json"),
           "translation": {k: {kk: v[kk] for kk in ("mismatch_counts", "collapsed_two_step_exchanges_verified", "collapsed_flipped",
                                                    "growth_tests_summary", "single_gene_deletions_all_equal", "genes_identical")}
                           for k, v in j("check_translation.json").items()},
           "gate": {k: v["values"] for k, v in j("check_energy_gate.json").items()},
           "gene_mapping": {k: {kk: v[kk] for kk in ("model_genes", "mapped", "mapped_with_fitness_data")} for k, v in j("check_gene_mapping.json").items()},
           "rescore_reproduction": j("rescore_iso783.json")["part_A_reproduction"]}
    dump(out, "summary.json")
    for r in rows:
        print(f"{r['quantity']:<42} reported={r['reported']!s:<24.24} reproduced={r['reproduced']!s:<24.24} corrected={r['corrected_gene_map']!s:<24.24} match={r['match']}")
    print("all reproduced:", out["all_reported_numbers_reproduced"])


if __name__ == "__main__":
    main()
