"""Supplementary check of the claim 'each such [capped] value sits at a cap set by the constraints' (Paper 2).

For every capped no-change error (equal positive maxima) on Harvey v0.3 and Harvetta v0.4, look for reactions of the
same metabolite (reaction ID contains the metabolite ID, case-insensitive, as a whole token) whose lower or upper bound
in the LP (MATLAB global bounds, shown bit-identical to the Python LP bounds by verify_v04_bounds.py) has the same
magnitude as the capped maximum (relative tolerance 1e-6). A match is evidence of a single binding cap; no match means
the cap, if any, is not a single bound of that metabolite (it may be a combination of constraints).
Imports nothing from gembench or scripts/. Writes caps_v04.json here.
"""
import json
import os
import re
import sys
from collections import OrderedDict

import numpy as np

from vcommon import ROOT, dump, load_bounds

GLOB = {"Harvey": "results/wbm_iem/matlab_reference/step2/matlab_global_bounds_Harvey_1_03d.mat",
        "Harvetta": "results/wbm_iem/matlab_reference/harvetta_step2/matlab_global_bounds_Harvetta_1_03d.mat"}
KEY = {"Harvey": "harvey_max_v0.3", "Harvetta": "harvetta_max_A"}


def met_of(rxn):
    m = re.match(r"^(?:EX|DM)_(.+?)\[(u|bc|csf|c)\]$", rxn)
    return m.group(1) if m else None


def main():
    scores = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "scores_v04.json")))
    out = OrderedDict()
    for model, path in GLOB.items():
        rx, lb, ub = load_bounds(os.path.join(ROOT, path))
        capped = [x for x in scores[KEY[model]]["no_change_errors"] if x[3] == "capped"]
        res = []
        for iem, rxn, exp, sub, h, d in capped:
            met = met_of(rxn)
            pat = re.compile(r"(^|[_(\[])" + re.escape(met) + r"($|[_(\[\]])", re.IGNORECASE)
            val = max(h, d)
            same_met = []
            any_rxn = 0
            for i, r in enumerate(rx):
                for bname, b in (("lb", lb[i]), ("ub", ub[i])):
                    if b != 0 and abs(abs(b) - val) <= 1e-6 * val:
                        any_rxn += 1
                        if pat.search(r):
                            same_met.append((r, bname, float(b)))
            res.append(dict(iem=iem, reaction=rxn, value=val, n_same_metabolite_bound_matches=len(same_met),
                            same_metabolite_matches=same_met[:6], n_any_reaction_bound_matches=any_rxn))
        out[model] = dict(n_capped=len(res), n_with_same_metabolite_bound=sum(r["n_same_metabolite_bound_matches"] > 0 for r in res),
                          n_with_any_bound=sum(r["n_any_reaction_bound_matches"] > 0 for r in res), rows=res)
        print(f"== {model}: capped {len(res)}; with a bound of the same metabolite equal to the value: "
              f"{out[model]['n_with_same_metabolite_bound']}; with any bound equal: {out[model]['n_with_any_bound']}")
        for r in res:
            print(f"   {r['iem']:6s} {r['reaction']:18s} {r['value']:12.6g}  same-met {r['n_same_metabolite_bound_matches']:2d} "
                  f"any {r['n_any_reaction_bound_matches']:4d}  {r['same_metabolite_matches'][:3]}")
    dump("caps_v04.json", out)


if __name__ == "__main__":
    sys.exit(main())
