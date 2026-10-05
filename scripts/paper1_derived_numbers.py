"""Derived numbers quoted in the paper 1 draft that are not stored in an existing result file.

Computed only from existing results (no simulation):
- genes and gene x condition calls changed by U' (B0 -> UNQ), by the menaquinone rule (U -> UQ), by removals and
  joins (UNQ -> M_R1R2), with their agreement with the fitness data (from paired.json);
- t-based 95% intervals of the organism-level mean differences (sensitivity to the percentile bootstrap);
- the genes whose calls the menaquinone rule changed, with their Fitness Browser descriptions (from the run matrices).

Usage: python scripts/paper1_derived_numbers.py  -> results/transfer_v1/paper1_derived_numbers.json
"""
from __future__ import annotations

import csv
import json
import math
import os
import statistics

import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "transfer_v1")
NEW = [("evaluation_panel_A", ["Cola", "Dino", "Dyella79", "MycoTube", "RPal_CGA009"]),
       ("evaluation_panel_B", ["Caulo", "Cup4G11", "Marino", "PV4", "SB2B"])]
GROWTH = 0.001


def paired(role):
    with open(os.path.join(R, role, "paired.json")) as fh:
        return json.load(fh)


def changes(a, b):
    out = {"genes_changed_common": 0, "to_important": 0, "to_important_agree": 0, "to_unimportant": 0,
           "to_unimportant_agree": 0, "per_organism": {}}
    for role, orgs in NEW:
        for p in paired(role):
            if p["org"] in orgs and p["A"] == a and p["B"] == b:
                out["genes_changed_common"] += p["genes_changed"]
                out["to_important"] += p["to_no_growth"]; out["to_important_agree"] += p["to_no_growth_agree"]
                out["to_unimportant"] += p["to_growth"]; out["to_unimportant_agree"] += p["to_growth_agree"]
                out["per_organism"][p["org"]] = {"genes_changed_common": p["genes_changed"], "union_delta": p["union_delta"]}
    return out


def t_interval(a, b):
    res = {}
    allx = []
    for role, orgs in NEW:
        xs = [p["union_delta"] for p in paired(role) if p["org"] in orgs and p["A"] == a and p["B"] == b]
        allx += xs
        res[role] = _t(xs)
    res["pooled"] = _t(allx)
    return res


def _t(xs):
    m, s, n = statistics.mean(xs), statistics.stdev(xs), len(xs)
    h = stats.t.ppf(0.975, n - 1) * s / math.sqrt(n)
    return {"n": n, "mean": m, "t95": [m - h, m + h]}


def menaquinone_genes():
    out = {}
    for role, orgs in NEW:
        for org in orgs:
            pa = os.path.join(R, role, org, "U", "matrices.npz")
            pb = os.path.join(R, role, org, "UQ", "matrices.npz")
            if not (os.path.exists(pa) and os.path.exists(pb)):
                continue
            a, b = np.load(pa), np.load(pb)
            ga = {g: i for i, g in enumerate(a["browser_genes"])}
            gb = {g: i for i, g in enumerate(b["browser_genes"])}
            ca = {c: i for i, c in enumerate(a["conditions"])}
            cb = {c: i for i, c in enumerate(b["conditions"])}
            conds = [c for c in a["conditions"] if c in cb and a["wt_growth"][ca[c]] > GROWTH and b["wt_growth"][cb[c]] > GROWTH]
            desc = {}
            with open(os.path.join(ROOT, "data", "fitness_browser_panel", org, "genes.tsv")) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    desc[r["locusId"]] = r["desc"]
            changed = []
            for g in sorted(set(ga) & set(gb)):
                if any((a["sim_growth"][ga[g], ca[c]] < GROWTH) != (b["sim_growth"][gb[g], cb[c]] < GROWTH) for c in conds):
                    changed.append({"locus": str(g), "description": desc.get(str(g))})
            if changed:
                out[org] = changed
    return out


def main() -> None:
    report = {"note": "Derived from existing transfer_v1 results only; see the module docstring.",
              "U_prime_vs_B0": changes("B0", "UNQ"), "menaquinone_rule_U_vs_UQ": changes("U", "UQ"),
              "removals_joins_UNQ_vs_M_R1R2": changes("UNQ", "M_R1R2"),
              "t_intervals": {"H1_U_prime_vs_B0": t_interval("B0", "UNQ"), "H2_M_vs_U_prime": t_interval("UNQ", "M")},
              "menaquinone_rule_changed_genes": menaquinone_genes()}
    report["menaquinone_rule_n_genes"] = sum(len(v) for v in report["menaquinone_rule_changed_genes"].values())
    out = os.path.join(R, "paper1_derived_numbers.json")
    with open(out, "w") as fh:
        json.dump(report, fh, indent=1)
        fh.write("\n")
    print(json.dumps({k: v for k, v in report.items() if k != "menaquinone_rule_changed_genes"}, indent=1)[:2500])


if __name__ == "__main__":
    main()
