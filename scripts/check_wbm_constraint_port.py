"""Check the constraint port against the bounds stored in a whole-body model file.

Harvey 1.03d ships with the physiological and diet constraints already applied (its SetupInfo records
the diet). Re-applying the port therefore shows (a) whether the port reproduces the stored bounds and
(b) exactly what re-application by the current Toolbox code changes. Two physiological variants are run:
  current  - physiologicalConstraintsHMDBbased at the pinned commit (GFR 90 ml/min, CSF export from
             CSFBloodFlowRate 0.52 ml/min);
  legacy   - GFR as 20 percent of renal plasma flow and CSF export upper bounds from CSFFlowRate
             (0.35 ml/min), as in physiologicalConstraintsHMDBbased_old.m.
The bile-duct step of runIEM_HH.m is compared with the v0.2 implementation as well.

Usage: python scripts/check_wbm_constraint_port.py Harvey_1_03d
Writes results/wbm_iem/<model>_constraint_port_check_v0.3.json
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from gembench import wbm as W  # noqa: E402
from gembench import wbm_constraints as C  # noqa: E402
from gembench import wbm_iem as I  # noqa: E402


def net_changes(rxns, lb0, ub0, lb1, ub1, log=None):
    last = {}
    for section, r, bound, _, _ in (log.changes if log else []):
        last[(r, bound)] = section
    out = []
    for name, a, b in (("lb", lb0, lb1), ("ub", ub0, ub1)):
        for i in np.where(a != b)[0]:
            out.append({"rxn": str(rxns[i]), "bound": name, "stored": float(a[i]), "port": float(b[i]),
                        "section": last.get((str(rxns[i]), name), "unknown")})
    return out


def group(changes):
    by = collections.defaultdict(list)
    for c in changes:
        by[c["section"]].append(c)
    out = {}
    for section, items in sorted(by.items()):
        ratios = collections.Counter(round(c["port"] / c["stored"], 4) if c["stored"] != 0 else None for c in items)
        out[section] = {"n": len(items), "bounds": dict(collections.Counter(c["bound"] for c in items)),
                        "ratio_port_to_stored": {str(k): v for k, v in ratios.most_common(5)},
                        "examples": items[:5] if len(items) > 12 else items}
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model")
    ap.add_argument("--model-file")
    args = ap.parse_args()
    model_file = args.model_file or os.path.join(ROOT, "external", "COBRA.models", "mat", f"{args.model}.mat")
    m = W.load_wbm(model_file)
    inputs = C.load_inputs()
    sex = m.meta.get("sex") or "male"
    report = {"model": args.model, "model_sha256": hashlib.sha256(open(model_file, "rb").read()).hexdigest(),
              "constraint_inputs_sha256": hashlib.sha256(open(C.DEFAULT_INPUTS, "rb").read()).hexdigest(),
              "toolbox_commit": inputs["provenance"]["toolbox_commit"], "sex": sex}

    params = C.default_parameters(inputs, sex)
    lb1, ub1, log1 = C.physiological_constraints(m.rxns, m.mets, m.S, m.lb, m.ub, params, inputs)
    report["physiology_current"] = group(net_changes(m.rxns, m.lb, m.ub, lb1, ub1, log1))
    legacy = C.default_parameters(inputs, sex)
    legacy.gfr = C.legacy_gfr(legacy)
    legacy.csf_export_ub_flow_rate = legacy.csf_flow_rate
    lb2, ub2, log2 = C.physiological_constraints(m.rxns, m.mets, m.S, m.lb, m.ub, legacy, inputs)
    report["physiology_legacy"] = {"gfr_ml_min": legacy.gfr, "csf_export_ub_flow_ml_min": legacy.csf_export_ub_flow_rate,
                                   "differences": net_changes(m.rxns, m.lb, m.ub, lb2, ub2, log2)}
    lb3, ub3, log3 = C.diet_constraints(m.rxns, m.lb, m.ub, inputs["diet_eu_average"], inputs["agora_essential"])
    report["diet_current"] = group(net_changes(m.rxns, m.lb, m.ub, lb3, ub3, log3))
    report["diet_warnings"] = log3.warnings
    lb4, ub4, rep4 = C.runiem_model_setup(m.rxns, m.mets, m.S, m.lb, m.ub, sex, inputs)
    report["full_setup_current"] = {"n_lb_changed": int((lb4 != m.lb).sum()), "n_ub_changed": int((ub4 != m.ub).sum()),
                                    "warnings": rep4["warnings"]}
    hw = I.HighsWBM(m)
    base_ub = hw.ub.copy()
    I.apply_runiem_global_constraints(hw, bile_duct="toolbox")
    toolbox_ub = hw.ub.copy()
    hw.set_bounds(range(len(base_ub)), lb=m.lb, ub=base_ub)
    I.apply_runiem_global_constraints(hw, bile_duct="v0.2_all")
    v02_ub = hw.ub.copy()
    bile = [i for i, r in enumerate(m.rxns) if r.startswith("BileDuct_EX_") and r.endswith("[bd]_[luSI]")]
    differ = [i for i in bile if toolbox_ub[i] != v02_ub[i]]
    report["bile_duct"] = {"n_bile_duct_exits": len(bile), "n_toolbox_list": len(I.RUNIEM_BILE_DUCT_UB100),
                           "n_differ_v02_vs_toolbox": len(differ),
                           "stored_ub_of_differing": dict(collections.Counter(str(base_ub[i]) for i in differ)),
                           "toolbox_list_stored_ub": dict(collections.Counter(str(base_ub[m.rxn_index(r)]) for r in I.RUNIEM_BILE_DUCT_UB100))}
    out = os.path.join(ROOT, "results", "wbm_iem", f"{args.model}_constraint_port_check_v0.3.json")
    with open(out, "w") as fh:
        json.dump(report, fh, indent=1)
        fh.write("\n")
    print(json.dumps({k: (v if k not in ("physiology_current", "diet_current") else {s: g["n"] for s, g in v.items()})
                      for k, v in report.items() if k != "physiology_legacy"}, indent=1))
    print("legacy differences:", [(d["rxn"], d["bound"], d["stored"], d["port"]) for d in report["physiology_legacy"]["differences"]])


if __name__ == "__main__":
    main()
