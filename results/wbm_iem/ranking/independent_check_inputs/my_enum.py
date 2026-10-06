import importlib.util, math, itertools
import numpy as np
spec = importlib.util.spec_from_file_location("R", "/home/claude/mma/scripts/iem_disease_ranking.py")
R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)
spec2 = importlib.util.spec_from_file_location("M", "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_ranking/my_ranking.py")
M = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(M)
for case in (101, 103):
    m, prof, panel = M.synth(n=9, nr=8, seed=case, n_profiles=7)
    cands, mine = M.my_scores(m, prof, "protocol", False, panel)
    names = list(mine.keys()); n = len(cands)
    # integer-scaled values: p_top5 = (5-a)/t clipped -> multiply by LCM of t's (<=9): 2520
    V = np.array([[int(M.my_rank(mine[nm], b)[5] * 2520) for b in range(n)] for nm in names], dtype=np.int64)
    obs = sum(int(M.my_rank(mine[nm], cands.index(nm))[5] * 2520) for nm in names)
    tot = np.zeros(1, dtype=np.int64)
    for row in V:                                # all n^k combinations, exact integers
        tot = (tot[:, None] + row[None, :]).ravel()
    pe = np.mean(tot >= obs)
    print(f"case {case}: {tot.size} combinations; exact P(top5 sum >= obs) = {pe:.6f}")
