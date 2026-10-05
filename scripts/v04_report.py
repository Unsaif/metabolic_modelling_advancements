"""IEM v0.4: run every analysis declared in docs/studies/wbm-iem-v0.4-plan.md and apply its reading rules.

Inputs (results/wbm_iem/): Harvetta_1_03d_iem_results_v0.4.json (run A), Harvey_1_03d_iem_results_v0.4_min.json
(run B), Harvetta_1_03d_iem_results_v0.4_min.json (run C), Harvey_1_03d_iem_results_v0.3.json (Harvey maxima), and
the MATLAB reference outputs in matlab_reference/harvetta_step2 (run M) and step2 (Harvey, v0.3).

Analyses: A1 accuracy and coverage; A2 Python vs MATLAB; A3 error anatomy; A4 effect sizes; A5 flux-range rules.
Readings: Q1 (a)-(e) on Harvetta, Q2 R1 and R2 on both models. Also reports whether every input is complete
(57 IEMs) and carries one run fingerprint.

Usage: python scripts/v04_report.py [--out results/wbm_iem/v0.4_report.json]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
W = os.path.join(ROOT, "results", "wbm_iem")
sys.path.insert(0, ROOT)


def _load_module(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, "scripts", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


anat = _load_module("iem_error_anatomy")
eff = _load_module("iem_effect_size_sensitivity")
rng = _load_module("iem_range_calls")
cmp_ = _load_module("compare_matlab_reference")


def load(name):
    with open(os.path.join(W, name)) as fh:
        return json.load(fh)


def run_status(results):
    fps = {r.get("run_fingerprint") for r in results}
    return {"n_iems": len(results), "complete": len(results) == 57, "fingerprints": sorted(f[:12] for f in fps if f),
            "one_fingerprint": len(fps) == 1, "statuses": sorted({r.get("status") for r in results})}


def effect_sizes(results):
    rows = list(eff.scored(results))
    out = {"n_scored": len(rows)}
    for tau in eff.TAUS:
        out[f"tau_{tau}"] = sum(eff.call(b["healthy"], b["disease"], tau) == b["expected"] for _, b in rows)
    small = 0
    for _, b in rows:
        if b["correct"]:
            scale = max(abs(b["healthy"]), abs(b["disease"]))
            rel = abs(b["disease"] - b["healthy"]) / scale if scale > 0 else 0.0
            small += rel < 0.05
    out["correct_below_5pct"] = small
    return out


def matlab_comparison(matlab_dir, model, python_results):
    sol = cmp_.parse_iemsol(os.path.join(matlab_dir, f"matlab_iem_results_{model}.mat"))
    res = cmp_.compare_results(sol, python_results)
    both_finite = [e for e in res["rows"] if e["matlab_finite"] and e.get("python_call") not in (None, "absent", "NA")]
    res["n_both_have_values"] = len(both_finite)
    res["n_same_call_both_have_values"] = sum(e["python_call"] == e["matlab_call"] for e in both_finite)
    nonfinite_directional = [e for e in res["rows"] if not e["matlab_finite"] and e["expected"] != "Unchanged"]
    res["n_matlab_nonfinite_directional"] = len(nonfinite_directional)
    res.pop("rows")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(W, "v0.4_report.json"))
    args = ap.parse_args()
    harvey_max = load("Harvey_1_03d_iem_results_v0.3.json")
    harvetta_max = load("Harvetta_1_03d_iem_results_v0.4.json")
    harvey_min = load("Harvey_1_03d_iem_results_v0.4_min.json")
    harvetta_min = load("Harvetta_1_03d_iem_results_v0.4_min.json")
    report = {"inputs": {"A_harvetta_max": run_status(harvetta_max), "B_harvey_min": run_status(harvey_min),
                         "C_harvetta_min": run_status(harvetta_min), "harvey_max_v0.3": run_status(harvey_max)}}
    # A1 + A3
    report["A1_A3"] = {"Harvey_v0.3": {k: v for k, v in anat.anatomy(harvey_max).items() if k.startswith("n_")},
                       "Harvetta_v0.4": {k: v for k, v in anat.anatomy(harvetta_max).items() if k.startswith("n_")}}
    hv = anat.anatomy(harvetta_max)
    report["A3_harvetta_lists"] = {"opposite": hv["opposite"], "no_change": hv["no_change"]}
    # A2
    report["A2"] = {"Harvetta": matlab_comparison(os.path.join(W, "matlab_reference", "harvetta_step2"), "Harvetta_1_03d", harvetta_max),
                    "Harvey": matlab_comparison(os.path.join(W, "matlab_reference", "step2"), "Harvey_1_03d", harvey_max)}
    # A4
    report["A4"] = {"Harvey_v0.3": effect_sizes(harvey_max), "Harvetta_v0.4": effect_sizes(harvetta_max)}
    # A5
    a5 = {}
    for label, mx, mn in (("Harvey", harvey_max, harvey_min), ("Harvetta", harvetta_max, harvetta_min)):
        r = rng.analyse(rng.pair_records(mx, mn))
        a5[label] = {k: r[k] for k in ("n_scored_four_optima", "n_directional_excluded_missing_or_nonoptimal", "correct",
                                       "net_gain_vs_max", "reading", "transitions_vs_max", "n_conflicting_R1",
                                       "equal_positive_maxima", "minima_differ")}
        a5[label]["changed_calls"] = r["changed_calls"]
    report["A5"] = a5
    # Readings
    acc = report["A1_A3"]["Harvetta_v0.4"]
    m2 = report["A2"]["Harvetta"]
    python_acc = acc["n_correct"] / acc["n_scored"] if acc["n_scored"] else None
    q1 = {
        "a_setup_identical": "replicated (checked before the freeze; see harvetta_step1_setup_bounds_comparison.json)",
        "b_same_calls_where_matlab_finite": ("replicated" if m2["n_same_call_both_have_values"] == m2["n_both_have_values"]
                                             else f"not replicated: {m2['n_both_have_values'] - m2['n_same_call_both_have_values']} differ"),
        "c_solver_failures_lower_runiem_figure": ("replicated" if m2["n_matlab_nonfinite_directional"] >= 1
                                                  and m2["matlab_accuracy_runiem_definition"] < python_acc else "not replicated"),
        "d_majority_no_change_errors_capped": ("replicated" if acc["n_no_change_equal_positive_maxima"] > acc["n_no_change"] / 2
                                               else "not replicated"),
        "e_correct_calls_below_5pct": report["A4"]["Harvetta_v0.4"]["correct_below_5pct"],
        "python_accuracy_harvetta": python_acc, "matlab_runiem_accuracy_harvetta": m2["matlab_accuracy_runiem_definition"]}
    r1 = [a5[m]["net_gain_vs_max"]["R1"] for m in ("Harvey", "Harvetta")]
    r1_reading = "supported" if all(g > 0 for g in r1) else "harmful" if any(g < 0 for g in r1) else "inconclusive"
    r2_useful = all(a5[m]["transitions_vs_max"]["R2"]["wrong_to_right"] > a5[m]["transitions_vs_max"]["R2"]["wrong_to_other_wrong"]
                    for m in ("Harvey", "Harvetta"))
    report["readings"] = {"Q1": q1, "Q2": {"R1": r1_reading, "R1_net_gain": dict(zip(("Harvey", "Harvetta"), r1)),
                                            "R2": "useful" if r2_useful else "not useful"}}
    with open(args.out, "w") as fh:
        json.dump(report, fh, indent=1, default=float)
        fh.write("\n")
    print(json.dumps({"inputs": report["inputs"], "A1_A3": report["A1_A3"], "A4": report["A4"],
                      "A2": {k: {kk: vv for kk, vv in v.items() if kk != "differences"} for k, v in report["A2"].items()},
                      "A5": {k: {kk: vv for kk, vv in v.items() if kk != "changed_calls"} for k, v in a5.items()},
                      "readings": report["readings"]}, indent=1, default=float))


if __name__ == "__main__":
    main()
