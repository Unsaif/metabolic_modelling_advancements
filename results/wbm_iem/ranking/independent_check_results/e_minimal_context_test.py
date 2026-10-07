"""Task E (exploratory): own-biomarker calls of a few IEMs in the 'minimal' context (no demand sink open except the
one maximised), HiGHS interior point for every solve, compared with the protocol-context calls in the matrix.
Uses the project's protocol code read-only; writes only to this scratch directory."""
import importlib.util
import json
import os
import sys
import time

sys.dont_write_bytecode = True
ROOT = "/home/claude/mma"
sys.path.insert(0, ROOT)
spec = importlib.util.spec_from_file_location("rx", os.path.join(ROOT, "scripts", "run_wbm_iem_cross.py"))
X = importlib.util.module_from_spec(spec); spec.loader.exec_module(X)
from gembench import wbm as W, wbm_iem as I, wbm_constraints as WC  # noqa: E402

iems = sys.argv[1:] or ["STAR", "LTC4S", "CYP21D"]
protocol = json.load(open(os.path.join(ROOT, "data", "iem", "iem_protocol_v0.2.json")))
extra = [l.strip() for l in open(os.path.join(ROOT, "data", "iem", "iem_ranking_extra_readouts_v0.2.txt")) if l.strip()]
panel = list(dict.fromkeys([rid for p in protocol for rid, _ in p["biomarkers"]] + extra))
M = {r["iem"]: r for r in json.load(open(os.path.join(ROOT, "results", "wbm_iem", "Harvey_1_03d_iem_cross_ranking_v1.json")))}

t0 = time.time()
m = W.load_wbm(os.path.join(ROOT, "external", "COBRA.models", "mat", "Harvey_1_03d.mat"))
inputs = WC.load_inputs(WC.DEFAULT_INPUTS)
lb, ub, _ = WC.runiem_model_setup(m.rxns, m.mets, m.S, m.lb, m.ub, sex=m.meta.get("sex") or "male", inputs=inputs)
m.lb, m.ub = lb, ub
hw = X.make_backend("highs", m, 2)
I.apply_runiem_global_constraints(hw, bile_duct="toolbox")
sink_mets = list(dict.fromkeys([rid[3:] for rid in panel if rid.startswith("DM_")] +
                               [met for p in protocol for met in (p.get("demand_metabolites") or [])]))
for met in sink_mets:
    if met in hw.met_pos:
        hw.add_demand(met, lb=0.0, ub=0.0)
print(f"set-up {time.time() - t0:.0f}s", flush=True)
out = {}
for name in iems:
    p = next(x for x in protocol if x["iem"] == name)
    own = [rid for rid, _ in p["biomarkers"]]
    t1 = time.time()
    rec = X.run_one(hw, p, own, "minimal", "ipm", log=False)
    prot = {e["reaction"]: e for e in M[name]["readouts"]}
    rows = []
    for e in rec["readouts"]:
        q = prot[e["reaction"]]
        rows.append({"reaction": e["reaction"], "expected": q["expected"], "protocol_call": q["predicted"],
                     "minimal_call": e["predicted"], "protocol": [q["healthy"], q["disease"]],
                     "minimal": [e["healthy"], e["disease"]]})
    out[name] = {"status": rec["status"], "vmax_minimal": rec["vmax_healthy"], "vmax_protocol": M[name]["vmax_healthy"],
                 "rows": rows, "time_s": round(time.time() - t1)}
    same = sum(r["protocol_call"] == r["minimal_call"] for r in rows)
    print(f"{name}: status {rec['status']} vmax minimal {rec['vmax_healthy']:.6g} protocol {M[name]['vmax_healthy']:.6g}; "
          f"own calls same {same}/{len(rows)}; {out[name]['time_s']}s", flush=True)
    for r in rows:
        print(f"   {r['reaction']:22s} expected {r['expected']:9s} protocol {r['protocol_call']:9s} minimal {r['minimal_call']:9s} "
              f"prot h/d {r['protocol'][0]:.6g}/{r['protocol'][1]:.6g} min h/d {r['minimal'][0]:.6g}/{r['minimal'][1]:.6g}", flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "e_minimal_context_test.json"), "w"), indent=1)
print(f"done {time.time() - t0:.0f}s")
