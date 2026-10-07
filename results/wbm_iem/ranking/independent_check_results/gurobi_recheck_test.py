"""Does a 'recheck' (set Method=2 then optimize() again, nothing else changed) re-solve in Gurobi?"""
import time
import numpy as np
import gurobipy as gp
from gurobipy import GRB

rng = np.random.default_rng(1)
m, n = 300, 600
A = rng.random((m, n)) * (rng.random((m, n)) < 0.05)
b = rng.random(m) * 10 + 1
c = rng.random(n)
env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
mod = gp.Model(env=env)
x = mod.addMVar(n, lb=0, ub=10)
mod.addMConstr(A, x, "<", b)
mod.setObjective(gp.quicksum([]) , GRB.MAXIMIZE)
mod.setAttr("Obj", x.tolist(), c.tolist())
mod.ModelSense = GRB.MAXIMIZE
mod.Params.Method = 2
mod.optimize()
print("first ipm: status", mod.Status, "obj", mod.ObjVal, "iters", mod.IterCount, "bar", mod.BarIterCount, "rt", mod.Runtime)
# warm solve after a bound change, as in the runner (concurrent)
x.tolist()[0].UB = 5.0
mod.Params.Method = 3
mod.optimize()
print("warm concurrent: status", mod.Status, "obj", mod.ObjVal, "iters", mod.IterCount, "bar", mod.BarIterCount, "rt", mod.Runtime)
v_warm = mod.ObjVal
# 'recheck' exactly as run_wbm_iem_cross.run_readout_major does: Method=2 and optimize(), no reset
mod.Params.Method = 2
t = time.time(); mod.optimize(); dt = time.time() - t
print("recheck as coded: status", mod.Status, "obj", mod.ObjVal, "iters", mod.IterCount, "bar", mod.BarIterCount, "rt", mod.Runtime, "wall", round(dt, 4), "same value", mod.ObjVal == v_warm)
# a real from-scratch recheck
mod.reset()
t = time.time(); mod.optimize(); dt = time.time() - t
print("recheck with reset(): status", mod.Status, "obj", mod.ObjVal, "iters", mod.IterCount, "bar", mod.BarIterCount, "rt", mod.Runtime, "wall", round(dt, 4))
