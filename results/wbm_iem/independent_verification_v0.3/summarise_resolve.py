"""Summarise resolve_nan_healthy.jsonl: my healthy optimum per MATLAB-NaN biomarker, an interval for the healthy maximum
[value of my feasible vector, weak-duality upper bound], and the direction implied by that interval and MATLAB's own
(finite) disease value with the 1e-6 rule. Writes resolve_summary.json next to this script."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    recs = [json.loads(l) for l in open(os.path.join(HERE, "resolve_nan_healthy.jsonl"))]
    mc = json.load(open(os.path.join(HERE, "matlab_comparison.json")))
    mat = {(f["iem"], f["reaction"]): f for f in mc["matlab_nonfinite"]}
    diffs = {(d["iem"], d["reaction"]): d for d in mc["python"]["v0.3"]["differences"]}
    rows, vm = [], []
    for r in recs:
        v = r["vmax"]
        vm.append({"iem": r["iem"], "status": v["status"], "mine": v["value"], "python_v0.3": v["python_v0.3"],
                   "abs_diff": abs(v["value"] - v["python_v0.3"]), "iem_rxns_equal_python": r["iem_rxns_equal_python_record"],
                   "max_row_violation": v.get("max_row_violation"), "max_bound_violation": v.get("max_bound_violation")})
        for b in r["biomarkers"]:
            key = (r["iem"], b["biomarker"])
            md = float(mat[key]["disease"])
            lo, hi = b["value"], b.get("lagrangian_upper_bound")
            # implied MATLAB-style call from healthy interval and MATLAB disease value
            if hi is not None and md - hi > 1e-6:
                implied = "Increased"
            elif md - lo < -1e-6:
                implied = "Decreased"
            else:
                implied = "undetermined"
            rows.append({"iem": r["iem"], "biomarker": b["biomarker"], "status": b["status"], "mine": lo,
                         "python_v0.3": b["python_v0.3_healthy"], "abs_diff": b["abs_diff_vs_python"],
                         "rel_diff": (b["abs_diff_vs_python"] / max(abs(lo), abs(b["python_v0.3_healthy"]))
                                      if max(abs(lo), abs(b["python_v0.3_healthy"])) > 0 else 0.0),
                         "upper_bound": hi, "max_row_violation": b["max_row_violation"],
                         "max_bound_violation": b["max_bound_violation"], "matlab_disease": md,
                         "python_call": diffs[key]["py_call"] if key in diffs else "(same as MATLAB)",
                         "implied_call": implied, "seconds": b["seconds"]})
    out = {"vmax": vm, "healthy": rows,
           "n_healthy": len(rows), "n_optimal": sum(r["status"] == "Optimal" for r in rows),
           "worst_row_violation": max(r["max_row_violation"] for r in rows),
           "worst_bound_violation": max(r["max_bound_violation"] for r in rows),
           "max_abs_diff_vs_v0.3": max(r["abs_diff"] for r in rows),
           "max_rel_diff_vs_v0.3_nonzero": max(r["rel_diff"] for r in rows if abs(r["python_v0.3"]) > 1e-3),
           "n_implied_equals_python_call": sum(r["implied_call"] == r["python_call"] for r in rows if r["python_call"] != "(same as MATLAB)"),
           "n_differing": sum(r["python_call"] != "(same as MATLAB)" for r in rows)}
    with open(os.path.join(HERE, "resolve_summary.json"), "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    for v in vm:
        print(f"vmax {v['iem']:6s} {v['status']} mine {v['mine']:.9f} v0.3 {v['python_v0.3']:.9f} diff {v['abs_diff']:.1e} "
              f"rxns equal {v['iem_rxns_equal_python']} row {v['max_row_violation']:.1e} bound {v['max_bound_violation']:.1e}")
    for r in rows:
        print(f"{r['iem']:6s} {r['biomarker']:15s} {r['status']} mine {r['mine']:.9g} v0.3 {r['python_v0.3']:.9g} diff {r['abs_diff']:.1e} "
              f"UB {r['upper_bound']:.4g} row {r['max_row_violation']:.1e} bound {r['max_bound_violation']:.1e} "
              f"MATLAB disease {r['matlab_disease']:.6g} py {r['python_call']} implied {r['implied_call']}")
    print({k: v for k, v in out.items() if k not in ("vmax", "healthy")})


if __name__ == "__main__":
    main()
