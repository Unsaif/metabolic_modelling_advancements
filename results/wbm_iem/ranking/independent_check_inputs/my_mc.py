import importlib.util, math, random, collections, sys
import numpy as np
sys.path.insert(0, "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_ranking")
spec = importlib.util.spec_from_file_location("R", "/home/claude/mma/scripts/iem_disease_ranking.py")
R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)
spec2 = importlib.util.spec_from_file_location("M", "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_ranking/my_ranking.py")
# reuse only synth() and my_scores/my_rank from my own module (importing runs its self-test quickly)
M = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(M)

def exact_tail(per_profile_values, obs, stat="mean"):
    """P(statistic >= obs) when each profile draws one of its n values uniformly (exact convolution)."""
    dist = {0.0: 1.0}
    for vals in per_profile_values:
        new = collections.defaultdict(float)
        for s, p in dist.items():
            for v in vals:
                new[round(s + v, 10)] += p / len(vals)
        dist = new
    k = len(per_profile_values)
    tgt = obs * k if stat == "mean" else obs
    return sum(p for s, p in dist.items() if s >= tgt - 1e-9)

for seed in (101, 102, 103, 104):
    m, prof, panel = M.synth(n=9, nr=8, seed=seed, n_profiles=7)
    order, calls, n_na = R.calls_from_matrix(m, "protocol")
    res = R.analyse(order, calls, n_na, prof, draws=100_000)
    s = res["summary"]
    cands, mine = M.my_scores(m, prof, "protocol", False, panel)
    rr_vals = [[float(M.my_rank(row, j)[3]) for j in range(len(cands))] for row in mine.values()]
    t1_vals = [[float(M.my_rank(row, j)[4]) for j in range(len(cands))] for row in mine.values()]
    t5_vals = [[float(M.my_rank(row, j)[5]) for j in range(len(cands))] for row in mine.values()]
    pe = exact_tail(rr_vals, s["mrr"]); p1 = exact_tail(t1_vals, s["expected_top1"], "sum"); p5 = exact_tail(t5_vals, s["expected_top5"], "sum")
    se = lambda p: math.sqrt(p * (1 - p) / 1e5)
    print(f"seed {seed}: MRR {s['mrr']:.4f} p_script {s['p_mrr']:.5f} p_exact {pe:.5f} (z={(s['p_mrr']-pe)/max(se(pe),1e-9):+.2f}); "
          f"top1 p {s['p_top1']:.5f} vs {p1:.5f}; top5 p {s['p_top5']:.5f} vs {p5:.5f}; null mean MRR {np.mean([np.mean(v) for v in rr_vals]):.6f} vs H_n/n {s['null_mrr']:.6f}")

# full-size: independent MC with Python's random module
m, prof, panel = M.synth(seed=7)
order, calls, n_na = R.calls_from_matrix(m, "protocol")
res = R.analyse(order, calls, n_na, prof, draws=100_000)
s = res["summary"]
cands, mine = M.my_scores(m, prof, "protocol", False, panel)
rr_vals = [[float(M.my_rank(row, j)[3]) for j in range(len(cands))] for row in mine.values()]
rnd = random.Random(424242); N = 200_000; ge = 0; tot = 0.0
for _ in range(N):
    v = sum(rnd.choice(vals) for vals in rr_vals) / len(rr_vals); tot += v
    ge += v >= s["mrr"] - 1e-12
pm = (1 + ge) / (1 + N)
print(f"57x57 synthetic: MRR {s['mrr']:.4f}, script p {s['p_mrr']:.5f}, my MC p {pm:.5f} (SE {math.sqrt(pm*(1-pm)/N):.5f}); my null mean {tot/N:.5f} vs {s['null_mrr']:.5f}")
# determinism with the fixed seed
res2 = R.analyse(order, calls, n_na, prof, draws=100_000)
print("same seed reproduces p:", res2["summary"]["p_mrr"] == s["p_mrr"], res2["summary"]["p_top5"] == s["p_top5"])
