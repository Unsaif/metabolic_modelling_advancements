"""Check results/wbm_iem/matlab_reference/certified_python_optima_v0.3.json (claim 7). No gembench imports.

For each biomarker whose MATLAB healthy value is not finite (list taken from my own matlab_comparison.json), check:
  healthy_solve.status == 'Optimal'; certificate max_row_violation <= 1e-6 and max_bound_violation <= 1e-6;
  certificate objective_recomputed equal to the healthy value reported in the v0.3 results file within 1e-6
  (absolute or relative), and to the certify run's own healthy value.
Also: whether the certify re-run reproduces the v0.3 values and calls for every biomarker of these IEMs, and the worst
violations over every solve in the file. Writes certified_check.json next to this script.
"""
import json
import math
import os
from collections import OrderedDict

ROOT = "/home/claude/mma"
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(ROOT, "results/wbm_iem/matlab_reference/certified_python_optima_v0.3.json")
NEEDED = ["EF", "ASNSD", "GMT", "PHOX1", "SUCLA", "BTD"]
TOL = 1e-6


def close(a, b):
    if a is None or b is None:
        return False
    d = abs(a - b)
    return d <= TOL or d <= TOL * max(abs(a), abs(b))


def main():
    cert = json.load(open(CERT))
    mc = json.load(open(os.path.join(HERE, "matlab_comparison.json")))
    v03 = {r["iem"]: r for r in json.load(open(os.path.join(ROOT, "results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json")))}
    nan_rows = [(f["iem"], f["reaction"]) for f in mc["matlab_nonfinite"] if f["healthy"] in ("NaN", "NA")]
    diff_rows = [(d["iem"], d["reaction"]) for d in mc["python"]["v0.3"]["differences"]]
    rep = OrderedDict(iems_present=list(cert["iems"]), missing=[i for i in NEEDED if i not in cert["iems"]],
                      n_matlab_healthy_nonfinite=len(nan_rows), n_differing_calls=len(diff_rows))
    rows = []
    for iem, rxn in nan_rows:
        if iem not in cert["iems"]:
            rows.append({"iem": iem, "biomarker": rxn, "available": False})
            continue
        b = next((x for x in cert["iems"][iem]["biomarkers"] if x["biomarker"] == rxn), None)
        if b is None:
            rows.append({"iem": iem, "biomarker": rxn, "available": False})
            continue
        hs = b["healthy_solve"] or {}
        c = hs.get("certificate") or {}
        reported = next(x["healthy"] for x in v03[iem]["biomarkers"] if x["reaction"] == rxn)
        rec = c.get("objective_recomputed")
        ok = (hs.get("status") == "Optimal" and c.get("max_row_violation", 1) <= TOL and c.get("max_bound_violation", 1) <= TOL
              and close(rec, reported))
        rows.append({"iem": iem, "biomarker": rxn, "available": True, "in_differing_calls": (iem, rxn) in diff_rows,
                     "status": hs.get("status"), "objective_reactions": hs.get("objective_reactions"),
                     "max_row_violation": c.get("max_row_violation"), "max_bound_violation": c.get("max_bound_violation"),
                     "objective_recomputed": rec, "v0.3_reported_healthy": reported,
                     "abs_diff_vs_v0.3": abs(rec - reported) if rec is not None and reported is not None else None,
                     "certify_run_healthy": b["healthy"], "certify_run_matches_v0.3": b["healthy"] == reported,
                     "passes": ok})
    rep["nan_biomarkers"] = rows
    avail = [r for r in rows if r["available"]]
    rep["n_available"] = len(avail)
    rep["n_pass"] = sum(r["passes"] for r in avail)
    rep["n_differing_available_pass"] = sum(r["passes"] for r in avail if r["in_differing_calls"])
    rep["worst_row_violation_nan_set"] = max((r["max_row_violation"] for r in avail), default=None)
    rep["worst_bound_violation_nan_set"] = max((r["max_bound_violation"] for r in avail), default=None)
    rep["worst_objective_abs_diff_nan_set"] = max((r["abs_diff_vs_v0.3"] for r in avail), default=None)
    # every solve in the file
    allc = []
    repro = []
    for iem, e in cert["iems"].items():
        for s in e.get("setup_solves", []):
            allc.append((iem, "setup", s))
        py = {x["reaction"]: x for x in v03[iem]["biomarkers"]}
        for b in e["biomarkers"]:
            for side in ("healthy_solve", "disease_solve"):
                if b.get(side):
                    allc.append((iem, f"{b['biomarker']} {side}", b[side]))
            p = py[b["biomarker"]]
            repro.append({"iem": iem, "biomarker": b["biomarker"], "healthy_equal": b["healthy"] == p["healthy"],
                          "disease_equal": b["disease"] == p["disease"], "call_equal": b["call"] == p["predicted"],
                          "healthy_abs_diff": None if b["healthy"] is None or p["healthy"] is None else abs(b["healthy"] - p["healthy"]),
                          "disease_abs_diff": None if b["disease"] is None or p["disease"] is None else abs(b["disease"] - p["disease"])})
        rep.setdefault("iem_status", {})[iem] = {"status": e["status"], "n_solves": e["n_solves"], "seconds": e["seconds"],
                                                  "vmax": e["setup_solves"][0]["objective"] if e.get("setup_solves") else None,
                                                  "vmax_v0.3": v03[iem]["vmax_healthy"]}
    certs = [s["certificate"] for _, _, s in allc if s.get("certificate")]
    rep["all_solves"] = {"n": len(allc), "n_with_certificate": len(certs),
                         "statuses": sorted({s["status"] for _, _, s in allc}),
                         "worst_row_violation": max(c["max_row_violation"] for c in certs) if certs else None,
                         "worst_bound_violation": max(c["max_bound_violation"] for c in certs) if certs else None,
                         "worst_objective_abs_diff": max(c["objective_abs_diff"] or 0 for c in certs) if certs else None}
    rep["certify_rerun_vs_v0.3"] = {"n": len(repro), "n_healthy_equal": sum(r["healthy_equal"] for r in repro),
                                    "n_disease_equal": sum(r["disease_equal"] for r in repro),
                                    "n_call_equal": sum(r["call_equal"] for r in repro),
                                    "max_healthy_abs_diff": max((r["healthy_abs_diff"] for r in repro if r["healthy_abs_diff"] is not None), default=None),
                                    "max_disease_abs_diff": max((r["disease_abs_diff"] for r in repro if r["disease_abs_diff"] is not None), default=None),
                                    "not_equal": [r for r in repro if not (r["healthy_equal"] and r["disease_equal"] and r["call_equal"])]}
    with open(os.path.join(HERE, "certified_check.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
        fh.write("\n")
    print(json.dumps({k: v for k, v in rep.items() if k not in ("nan_biomarkers",)}, indent=1)[:5000])
    for r in rows:
        if r["available"]:
            print(f"  {r['iem']:6s} {r['biomarker']:16s} {r['status']} row {r['max_row_violation']:.2e} bound {r['max_bound_violation']:.2e} "
                  f"obj {r['objective_recomputed']!r} v0.3 {r['v0.3_reported_healthy']!r} diff {r['abs_diff_vs_v0.3']:.2e} "
                  f"differing {r['in_differing_calls']} pass {r['passes']}")
        else:
            print(f"  {r['iem']:6s} {r['biomarker']:16s} not available yet")


if __name__ == "__main__":
    main()
