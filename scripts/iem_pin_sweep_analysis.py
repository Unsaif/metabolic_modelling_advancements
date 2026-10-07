"""Development analysis of the pin sweep: calls, own-biomarker matches and the disease ranking at each pin level.

Usage:
  python scripts/iem_pin_sweep_analysis.py [results/wbm_iem/pin_sweep/Harvey_1_03d_pin_sweep_development_v1.json]

Development only (docs/studies/wbm-iem-healthy-reference-dev-plan.md): the IEMs are those of the sweep (the
development set); candidates are those IEMs' knockouts and profiles are those IEMs' profiles. For each alpha:
  - calls: the protocol's rule on the sweep's healthy maximum and the matrix's disease maximum (the disease state does
    not depend on alpha); the material rule as a sensitivity;
  - each knockout's share of panel readouts called increased and decreased;
  - own-biomarker matches for the lab and HPO profiles, with increases and decreases separately and balanced accuracy;
  - the ranking among the development candidates (plain and adjusted scores, as in scripts/iem_disease_ranking.py,
    with 10,000 Monte Carlo draws).
The alpha = 1 calls are compared with the matrix's, as a check of the rebuild.
Writes results/wbm_iem/pin_sweep/<sweep stem>_analysis.json.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("iem_disease_ranking", os.path.join(ROOT, "scripts", "iem_disease_ranking.py"))
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)
TOL = 1e-6


def protocol_call(h, d):
    if h is None or d is None or not (math.isfinite(h) and math.isfinite(d)):
        return None
    diff = d - h
    return 1 if diff > TOL else -1 if diff < -TOL else 0


def calls_at(sweep, matrix_by_iem, iems, alpha, rule="protocol"):
    """{iem: {readout: +1/0/-1}} and NA counts, from the sweep's healthy maxima at alpha (one value for every IEM, or
    {iem: alpha}) and the matrix's disease maxima."""
    calls, n_na = {}, {}
    for iem in iems:
        key = float(alpha[iem] if isinstance(alpha, dict) else alpha)
        disease = {e["reaction"]: e for e in matrix_by_iem[iem]["readouts"]}
        c, na = {}, 0
        for rid, per_iem in sweep["readouts"].items():
            entry = next((x for x in per_iem[iem] if float(x["alpha"]) == key), None)
            dz = disease.get(rid)
            if entry is None or dz is None or entry["status"] != "Optimal" or dz["status_disease"] != "Optimal":
                na += 1
                continue
            h, d = entry["value"], dz["disease"]
            v = protocol_call(h, d) if rule == "protocol" else R.material_call(h, d)
            if v is None:
                na += 1
            else:
                c[rid] = v
        calls[iem], n_na[iem] = c, na
    return calls, n_na


def own_matches(calls, profiles):
    out = {"n": 0, "correct": 0, "inc_n": 0, "inc_correct": 0, "dec_n": 0, "dec_correct": 0}
    for iem, prof in profiles.items():
        if iem not in calls:
            continue
        for rid, s in prof:
            if rid not in calls[iem]:
                continue
            ok = calls[iem][rid] == R.SIGN[s]
            out["n"] += 1; out["correct"] += ok
            k = "inc" if s == "Increased" else "dec"
            out[f"{k}_n"] += 1; out[f"{k}_correct"] += ok
    rec_inc = out["inc_correct"] / out["inc_n"] if out["inc_n"] else None
    rec_dec = out["dec_correct"] / out["dec_n"] if out["dec_n"] else None
    out["accuracy"] = out["correct"] / out["n"] if out["n"] else None
    out["balanced_accuracy"] = (rec_inc + rec_dec) / 2 if rec_inc is not None and rec_dec is not None else None
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep", nargs="?", default=os.path.join(ROOT, "results", "wbm_iem", "pin_sweep",
                                                              "Harvey_1_03d_pin_sweep_development_v1.json"))
    ap.add_argument("--profiles", default=os.path.join(ROOT, "data", "iem", "iem_ranking_profiles_v0.2.json"))
    ap.add_argument("--draws", type=int, default=10_000)
    args = ap.parse_args()
    with open(args.sweep) as fh:
        sweep = json.load(fh)
    prov = sweep["provenance"]
    with open(os.path.join(ROOT, prov["matrix"])) as fh:
        matrix = json.load(fh)
    if R.sha256(os.path.join(ROOT, prov["matrix"])) != prov["matrix_sha256"]:
        raise SystemExit("the matrix differs from the one the sweep was built from")
    by_iem = {r["iem"]: r for r in matrix}
    iems = sorted(prov["iems"], key=lambda i: by_iem[i]["call_index"])
    with open(args.profiles) as fh:
        P = json.load(fh)
    panel = list(sweep["readouts"])
    out = {"sweep": os.path.relpath(os.path.abspath(args.sweep), ROOT), "iems": iems, "alphas": prov["alphas"],
           "candidates": len(iems), "check_alpha_1": {}, "by_alpha": {}}

    # Check: alpha = 1 against the matrix.
    calls1, _ = calls_at(sweep, by_iem, iems, 1.0)
    diffs, differing = [], []
    for iem in iems:
        m = {e["reaction"]: e for e in by_iem[iem]["readouts"]}
        for rid, per_iem in sweep["readouts"].items():
            e = next(x for x in per_iem[iem] if float(x["alpha"]) == 1.0)
            if e["status"] == "Optimal" and m[rid]["status_healthy"] == "Optimal":
                diffs.append(abs(e["value"] - m[rid]["healthy"]))
                mc = {"Increased": 1, "Decreased": -1, "Unchanged": 0}.get(m[rid]["predicted"])
                if rid in calls1[iem] and calls1[iem][rid] != mc:
                    differing.append([iem, rid, m[rid]["healthy"], e["value"], m[rid]["disease"]])
    out["check_alpha_1"] = {"n": len(diffs), "max_abs_diff": max(diffs) if diffs else None,
                            "n_abs_diff_above_1e-6": sum(d > 1e-6 for d in diffs), "differing_calls": differing}

    order = iems
    for alpha in prov["alphas"]:
        res = {}
        for rule in ("protocol", "material"):
            calls, n_na = calls_at(sweep, by_iem, iems, alpha, rule)
            shares = {iem: {"increased": sum(v == 1 for v in calls[iem].values()) / max(1, len(calls[iem])),
                            "decreased": sum(v == -1 for v in calls[iem].values()) / max(1, len(calls[iem]))}
                      for iem in iems}
            block = {"median_share_increased": float(np.median([s["increased"] for s in shares.values()])),
                     "median_share_decreased": float(np.median([s["decreased"] for s in shares.values()])),
                     "shares": shares, "na": n_na}
            for pset in ("lab", "hpo"):
                profiles = {i: p for i, p in P["profiles"][pset].items() if i in iems}
                block[f"own_{pset}"] = own_matches(calls, profiles)
                for adjusted in (False, True):
                    a = R.analyse(order, calls, n_na, profiles, adjusted=adjusted, panel=panel, draws=args.draws)
                    s = a["summary"]
                    block[f"ranking_{pset}_{'adjusted' if adjusted else 'plain'}"] = {
                        k: s[k] for k in ("n_candidates", "n_profiles", "mrr", "null_mrr", "p_mrr", "expected_top1",
                                          "expected_top5", "median_expected_rank")}
                    block[f"ranks_{pset}_{'adjusted' if adjusted else 'plain'}"] = {
                        x["iem"]: x["expected_rank"] for x in a["profiles"]}
            res[rule] = block
        out["by_alpha"][str(alpha)] = res
        pr = res["protocol"]
        print(f"alpha {alpha:>6}: increased {pr['median_share_increased']:.2f} decreased {pr['median_share_decreased']:.3f} | "
              f"own lab {pr['own_lab']['correct']}/{pr['own_lab']['n']} (inc {pr['own_lab']['inc_correct']}/{pr['own_lab']['inc_n']}, "
              f"dec {pr['own_lab']['dec_correct']}/{pr['own_lab']['dec_n']}) | MRR lab plain {pr['ranking_lab_plain']['mrr']:.3f} "
              f"adjusted {pr['ranking_lab_adjusted']['mrr']:.3f} | hpo plain {pr['ranking_hpo_plain']['mrr']:.3f} "
              f"adjusted {pr['ranking_hpo_adjusted']['mrr']:.3f} (chance {pr['ranking_lab_plain']['null_mrr']:.3f})", flush=True)
    # A capped pin (development view): each IEM at the largest swept alpha whose pin is at most P (alpha = 1 when
    # v_max <= P; the smallest swept alpha when even that pin is above P).
    swept = sorted((a for a in prov["alphas"] if a > 0), reverse=True)
    out["cap"] = {}
    for cap in (10, 30, 100, 300, 1000):
        alpha_by = {}
        for iem in iems:
            vmax = by_iem[iem]["vmax_healthy"]
            fit = [a for a in swept if a * vmax <= cap]
            alpha_by[iem] = fit[0] if fit else swept[-1]
        res = {"alpha_by_iem": alpha_by}
        for rule in ("protocol", "material"):
            calls, n_na = calls_at(sweep, by_iem, iems, alpha_by, rule)
            block = {"median_share_increased": float(np.median([sum(v == 1 for v in calls[i].values()) / max(1, len(calls[i]))
                                                                for i in iems]))}
            for pset in ("lab", "hpo"):
                profiles = {i: p for i, p in P["profiles"][pset].items() if i in iems}
                block[f"own_{pset}"] = own_matches(calls, profiles)
                for adjusted in (False, True):
                    a = R.analyse(order, calls, n_na, profiles, adjusted=adjusted, panel=panel, draws=args.draws)
                    block[f"ranking_{pset}_{'adjusted' if adjusted else 'plain'}"] = {
                        k: a["summary"][k] for k in ("mrr", "null_mrr", "expected_top1", "expected_top5", "median_expected_rank")}
                    block[f"ranks_{pset}_{'adjusted' if adjusted else 'plain'}"] = {x["iem"]: x["expected_rank"] for x in a["profiles"]}
            res[rule] = block
        out["cap"][str(cap)] = res
        pr = res["protocol"]
        print(f"cap {cap:>5}: increased {pr['median_share_increased']:.2f} | own lab {pr['own_lab']['correct']}/{pr['own_lab']['n']} "
              f"(inc {pr['own_lab']['inc_correct']}/{pr['own_lab']['inc_n']}, dec {pr['own_lab']['dec_correct']}/{pr['own_lab']['dec_n']}) | "
              f"MRR lab plain {pr['ranking_lab_plain']['mrr']:.3f} adjusted {pr['ranking_lab_adjusted']['mrr']:.3f} | hpo plain "
              f"{pr['ranking_hpo_plain']['mrr']:.3f} adjusted {pr['ranking_hpo_adjusted']['mrr']:.3f}", flush=True)
    c = out["check_alpha_1"]
    print(f"check alpha 1 vs matrix: {c['n']} values, max |diff| {c['max_abs_diff']}, calls differing {len(c['differing_calls'])}")
    stem = os.path.splitext(os.path.basename(args.sweep))[0]
    path = os.path.join(os.path.dirname(os.path.abspath(args.sweep)), f"{stem}_analysis.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"-> {os.path.relpath(path, ROOT)}")


if __name__ == "__main__":
    main()
