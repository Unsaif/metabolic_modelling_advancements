"""Task E: does every IEM reaction's healthy-state bound interval contain 0? (Then, without the pin, the disease
state's feasible set lies inside the healthy state's and no increase is possible.) Uses the project's loader and
Toolbox set-up read-only; no LP is solved."""
import json, sys
sys.path.insert(0, "/home/claude/mma")
import numpy as np
from gembench import wbm as W, wbm_iem as I, wbm_constraints as WC
m = W.load_wbm("/home/claude/mma/external/COBRA.models/mat/Harvey_1_03d.mat")
inputs = WC.load_inputs(WC.DEFAULT_INPUTS)
lb, ub, rep = WC.runiem_model_setup(m.rxns, m.mets, m.S, m.lb, m.ub, sex=m.meta.get("sex") or "male", inputs=inputs)
M = json.load(open("/home/claude/mma/results/wbm_iem/Harvey_1_03d_iem_cross_ranking_v1.json"))
pos = {r: i for i, r in enumerate(m.rxns)}
bad = []
n = 0
for r in M:
    for rx in r["iem_reactions"]:
        i = pos[rx]; n += 1
        if not (lb[i] <= 0 <= ub[i]):
            bad.append((r["iem"], rx, lb[i], ub[i]))
print("IEM reactions checked", n, "with 0 outside [lb, ub] after the Toolbox set-up:", bad)
# global constraints only set lb=0, lb=ub=0 or ub=100, and bound tweaks only set lb=0 / ub=0 (see protocol), so 0 stays inside.
print("any model reaction with lb > 0:", int((lb > 0).sum()), "ub < 0:", int((ub < 0).sum()))
iem_set = {rx for r in M for rx in r["iem_reactions"]}
print("IEM reactions among lb>0 or ub<0:", [m.rxns[i] for i in np.where((lb > 0) | (ub < 0))[0] if m.rxns[i] in iem_set])
# coupling constraints: are they all homogeneous (rhs 0), so v_IEM = 0 keeps them satisfiable?
for attr in ("C", "d", "dsense", "ctrs"):
    print(attr, hasattr(m, attr))
