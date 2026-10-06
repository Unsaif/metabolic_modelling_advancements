import importlib.util, math
import numpy as np
spec = importlib.util.spec_from_file_location("R", "/home/claude/mma/scripts/iem_disease_ranking.py")
R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)
spec2 = importlib.util.spec_from_file_location("M", "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_ranking/my_ranking.py")
M = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(M)
m, prof, panel = M.synth(n=9, nr=8, seed=101, n_profiles=7)
order, calls, n_na = R.calls_from_matrix(m, "protocol")
names, S = R.score_matrix(order, calls, prof)
cands, mine = M.my_scores(m, prof, "protocol", False, panel)
Rf = np.array([[R.rank_stats(S[a], b)["p_top5"] for b in range(9)] for a in range(len(names))])
Ri = np.array([[int(M.my_rank(mine[nm], b)[5] * 2520) for b in range(9)] for nm in names])
obs_i = sum(int(M.my_rank(mine[nm], cands.index(nm))[5] * 2520) for nm in names)
res = R.analyse(order, calls, n_na, prof, draws=1_000_000)["summary"]
rng = np.random.default_rng(20261006)
idx = rng.integers(0, 9, size=(1_000_000, len(names)))
rows = np.arange(len(names))
f = Rf[rows, idx].sum(axis=1); i = Ri[rows, idx].sum(axis=1)
pf = (1 + np.sum(f >= res["expected_top5"] - 1e-12)) / (1 + 1e6); pi = (1 + np.sum(i >= obs_i)) / (1 + 1e6)
print("script p", res["p_top5"], "replayed float p", pf, "replayed exact-integer p", pi, "misclassified draws:", int(np.sum((f >= res['expected_top5'] - 1e-12) != (i >= obs_i))))
# distribution of z over many independent seeds for this case (exact p known)
pe = 0.8191501972937729
zs = []
for sd in range(100, 140):
    idx = np.random.default_rng(sd).integers(0, 9, size=(200_000, len(names)))
    p = np.mean(Ri[rows, idx].sum(axis=1) >= obs_i)
    zs.append((p - pe) / math.sqrt(pe * (1 - pe) / 200_000))
print("40 independent seeds: mean z %.2f, sd z %.2f" % (np.mean(zs), np.std(zs)))
