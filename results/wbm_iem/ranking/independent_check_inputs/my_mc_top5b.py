import importlib.util, math, collections
import numpy as np
from fractions import Fraction
spec = importlib.util.spec_from_file_location("R", "/home/claude/mma/scripts/iem_disease_ranking.py")
R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)
spec2 = importlib.util.spec_from_file_location("M", "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_ranking/my_ranking.py")
M = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(M)
m, prof, panel = M.synth(n=9, nr=8, seed=101, n_profiles=7)
order, calls, n_na = R.calls_from_matrix(m, "protocol")
names, S = R.score_matrix(order, calls, prof)
cands, mine = M.my_scores(m, prof, "protocol", False, panel)
# element-wise comparison of the null value tables
Rs = np.array([[R.rank_stats(S[a], b)["p_top5"] for b in range(len(order))] for a in range(len(names))])
Rm = np.array([[float(M.my_rank(mine[nm], b)[5]) for b in range(len(order))] for nm in names])
print("null table max |diff|:", np.abs(Rs - Rm).max(), "; names order equal:", names == list(mine.keys()))
obs = sum(M.my_rank(mine[nm], cands.index(nm))[5] for nm in names)
# exact distribution
dist = {Fraction(0): Fraction(1)}
for nm in names:
    vals = [M.my_rank(mine[nm], b)[5] for b in range(len(cands))]
    new = collections.defaultdict(Fraction)
    for s, p in dist.items():
        for v in vals: new[s + v] += p / len(vals)
    dist = new
pe = float(sum(p for s, p in dist.items() if s >= obs)); patom = float(dist.get(obs, 0))
print("obs", obs, float(obs), "exact P(>=obs)", pe, "P(==obs)", patom)
# independent MC, numpy Generator PCG64 with other seeds and 4e6 draws
for sd in (11, 12):
    rng = np.random.default_rng(sd); N = 4_000_000
    idx = rng.integers(0, len(order), size=(N, len(names)))
    null = Rm[np.arange(len(names)), idx].sum(axis=1)
    p = (1 + np.sum(null >= float(obs) - 1e-12)) / (1 + N)
    print(f"my MC seed {sd}: p {p:.5f} z={(p-pe)/math.sqrt(pe*(1-pe)/N):+.2f}; P(==obs) MC {np.mean(np.abs(null-float(obs))<1e-12):.5f}")
res = R.analyse(order, calls, n_na, prof, draws=4_000_000, seed=13)["summary"]
print(f"script seed 13, 4e6 draws: p_top5 {res['p_top5']:.5f} z={(res['p_top5']-pe)/math.sqrt(pe*(1-pe)/4e6):+.2f}")
