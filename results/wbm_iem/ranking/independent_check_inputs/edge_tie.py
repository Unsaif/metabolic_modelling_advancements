import importlib.util
spec = importlib.util.spec_from_file_location("R", "/home/claude/mma/scripts/iem_disease_ranking.py")
R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)
def ro(r, pred): return {"reaction": r, "predicted": pred, "healthy": 0.0, "disease": 0.0, "status_healthy": "Optimal", "status_disease": "Optimal"}
panel = ["P"] + [f"x{i}" for i in range(11)]
A = [ro("P", "Increased")] + [ro(f"x{i}", "Increased") for i in range(6)] + [ro(f"x{i}", "Decreased") for i in range(6, 11)]   # 7 up, 5 down of 12
B = [ro("P", "Increased")] + [ro("x0", "Increased")] + [ro(f"x{i}", "Unchanged") for i in range(1, 11)]                        # 2 up, 0 down of 12
mat = [{"iem": "A", "call_index": 1, "readouts": A}, {"iem": "B", "call_index": 2, "readouts": B}]
order, calls, n_na = R.calls_from_matrix(mat, "protocol")
names, S = R.score_matrix(order, calls, {"A": [["P", "Increased"]], "B": [["P", "Increased"]]}, adjusted=True, panel=panel)
print("adjusted scores", S[0].tolist(), "float diff", S[0][0] - S[0][1])
print("A:", R.rank_stats(S[0], 0)); print("B:", R.rank_stats(S[1], 1))
