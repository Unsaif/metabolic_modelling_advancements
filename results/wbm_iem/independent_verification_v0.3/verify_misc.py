"""Record-level cross-checks against the protocol and the MATLAB header rows. No gembench imports.

1. For every IEM record in each run: the IEM reaction set recomputed from the protocol patterns (substring match on the
   model's reaction ids, minus exclude patterns and 'Micro_') equals the record's iem_reactions.
2. MATLAB header rows vs the v0.3 records: 'WB obj - Disease' (feasibility) and 'IEM Rxns All obj - Disease'.
3. Per-run IEM status counts and the sum of n_solves.
Writes misc_checks.json next to this script.
"""
import json
import os
from collections import Counter, OrderedDict

import numpy as np
import scipy.io as sio

ROOT = "/home/claude/mma"
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = OrderedDict([
    ("v0.2", "results/wbm_iem/Harvey_1_03d_iem_results_v0.2.json"),
    ("v0.2b", "results/wbm_iem/Harvey_1_03d_iem_results_v0.2b.json"),
    ("v0.3", "results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json"),
])


def cellstr(a):
    return [str(np.asarray(x).ravel()[0]) if np.asarray(x).size else "" for x in np.asarray(a, dtype=object).ravel()]


def main():
    g = sio.loadmat(os.path.join(ROOT, "results/wbm_iem/matlab_reference/step2/matlab_global_bounds_Harvey_1_03d.mat"))
    rxns = cellstr(g["rxns"])        # same order as the model file (checked in verify_bounds.py)
    protocol = json.load(open(os.path.join(ROOT, "data/iem/iem_protocol_v0.2.json")))
    rep = OrderedDict()
    expected = {}
    for p in protocol:
        s = [r for r in rxns if any(q in r for q in p["include_patterns"])]
        s = [r for r in s if not any(q in r for q in p["exclude_patterns"]) and "Micro_" not in r]
        expected[p["iem"]] = sorted(s)
    for label, path in RUNS.items():
        recs = json.load(open(os.path.join(ROOT, path)))
        bad = [r["iem"] for r in recs if sorted(r["iem_reactions"]) != expected[r["iem"]]]
        rep[label] = {"n_records": len(recs), "iem_reaction_set_mismatch": bad,
                      "n_iem_reactions_total": sum(len(r["iem_reactions"]) for r in recs),
                      "status": dict(Counter(r["status"] for r in recs)),
                      "n_solves_total": sum(r["n_solves"] for r in recs),
                      "wb_objective_disease_feasible": dict(Counter(r["wb_objective_disease_feasible"] for r in recs)),
                      "vmax_disease_nonzero": [(r["iem"], r["vmax_disease"]) for r in recs
                                               if r["vmax_disease"] is None or abs(r["vmax_disease"]) > 1e-6]}
    d = sio.loadmat(os.path.join(ROOT, "results/wbm_iem/matlab_reference/step2/matlab_iem_results_Harvey_1_03d.mat"))
    hdr = {}
    for k in d:
        if k.startswith("IEMSol_"):
            c = np.asarray(d[k], dtype=object)
            txt = lambda x: str(np.asarray(x).ravel()[0]) if np.asarray(x).size else ""
            hdr[k[7:]] = {"vmax_disease": txt(c[1, 1]), "wb_disease": txt(c[3, 1]), "wb_disease_status": txt(c[3, 2])}
    rep["matlab_header"] = {"wb_disease_values": dict(Counter(h["wb_disease"] for h in hdr.values())),
                            "wb_disease_status": dict(Counter(h["wb_disease_status"] for h in hdr.values())),
                            "vmax_disease_values": dict(Counter(h["vmax_disease"] for h in hdr.values()))}
    with open(os.path.join(HERE, "misc_checks.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
        fh.write("\n")
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
