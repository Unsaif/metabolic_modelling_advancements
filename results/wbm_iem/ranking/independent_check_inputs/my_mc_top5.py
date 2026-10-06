import importlib.util, math, collections
spec = importlib.util.spec_from_file_location("R", "/home/claude/mma/scripts/iem_disease_ranking.py")
R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)
spec2 = importlib.util.spec_from_file_location("M", "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_ranking/my_ranking.py")
M = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(M)
from fractions import Fraction
def exact_tail(vals_list, obs):
    dist = {Fraction(0): Fraction(1)}
    for vals in vals_list:
        new = collections.defaultdict(Fraction)
        for s, p in dist.items():
            for v in vals: new[s + v] += p / len(vals)
        dist = new
    return float(sum(p for s, p in dist.items() if s >= obs))
for seed in (101, 103):
    m, prof, panel = M.synth(n=9, nr=8, seed=seed, n_profiles=7)
    order, calls, n_na = R.calls_from_matrix(m, "protocol")
    cands, mine = M.my_scores(m, prof, "protocol", False, panel)
    t5 = [[M.my_rank(row, j)[5] for j in range(len(cands))] for row in mine.values()]
    obs = sum(M.my_rank(row, cands.index(tr))[5] for tr, row in mine.items())
    pe = exact_tail(t5, obs)
    for s2, dr in ((20261006, 1_000_000), (1, 1_000_000), (2, 1_000_000)):
        res = R.analyse(order, calls, n_na, prof, draws=dr, seed=s2)["summary"]
        print(f"case {seed} seed {s2} draws {dr}: p_top5 script {res['p_top5']:.5f} exact {pe:.5f} z={(res['p_top5']-pe)/math.sqrt(pe*(1-pe)/dr):+.2f}")
