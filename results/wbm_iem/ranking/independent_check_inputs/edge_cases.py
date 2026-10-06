import importlib.util, numpy as np
spec = importlib.util.spec_from_file_location("R", "/home/claude/mma/scripts/iem_disease_ranking.py")
R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)
def ro(r, pred): return {"reaction": r, "predicted": pred, "healthy": 0.0, "disease": 0.0, "status_healthy": "Optimal", "status_disease": "Optimal"}
# float tie: A has 7/10 up, 5/10 down (0.7-0.5 = 0.19999999999999996); B has 2/10 up (0.2). Both match the profile readout P.
A = [ro("P", "Increased")] + [ro(f"a{i}", "Increased") for i in range(6)] + [ro(f"b{i}", "Decreased") for i in range(5)]
B = [ro("P", "Increased")] + [ro(f"a{i}", "Increased") for i in range(1)] + [ro(f"b{i}", "Unchanged") for i in range(10)]
A = A[:10]; B = B[:10]
mat = [{"iem": "A", "call_index": 1, "readouts": A}, {"iem": "B", "call_index": 2, "readouts": B}, {"iem": "C", "call_index": 3, "readouts": [ro("P", "Unchanged")]}]
order, calls, n_na = R.calls_from_matrix(mat, "protocol")
panel = sorted({x["reaction"] for r in mat for x in r["readouts"]})
names, S = R.score_matrix(order, calls, {"A": [["P", "Increased"]]}, adjusted=True, panel=panel)
pa = [sum(1 for v in calls[c].values() if v == 1) / len(calls[c]) - sum(1 for v in calls[c].values() if v == -1) / len(calls[c]) for c in ("A", "B")]
print("p_up-p_down A,B:", pa, "adjusted scores:", S[0].tolist(), "exact diff A-B:", S[0][0] - S[0][1])
print("rank_stats A:", R.rank_stats(S[0], 0))
# material-rule boundaries
for h, d in [(0.0, 1e-3), (0.0, 1.0011e-3), (1.0, 1.0526316), (1.0, 1.05), (1.0, 1.06), (0.5, 0.0), (0.0005, 0.0), (-2.0, -1.0), (None, 1.0)]:
    print(f"material_call(h={h}, d={d}) = {R.material_call(h, d)}")
# all-NA profile
mat2 = [{"iem": f"D{i}", "call_index": i, "readouts": [{"reaction": "X", "predicted": "NA", "healthy": None, "disease": None,
         "status_healthy": "absent", "status_disease": "absent"}]} for i in range(57)]
order, calls, n_na = R.calls_from_matrix(mat2, "material")
out = R.analyse(order, calls, n_na, {"D5": [["X", "Increased"]]}, draws=1000)
p = out["profiles"][0]; H57 = sum(1 / i for i in range(1, 58))
print("all-NA profile:", p["higher"], p["tied"], p["expected_rank"], round(p["expected_reciprocal_rank"], 6), round(H57 / 57, 6), round(p["p_top1"], 6), "n_listed", len(p["candidates_scoring_at_least_as_high"]), "of", p["n_candidates_scoring_at_least_as_high"])
