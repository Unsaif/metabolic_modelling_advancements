"""Development comparison of candidate rules on the full pin sweep (development IEMs only; 19 candidates)."""
import json, importlib.util, math, os
import numpy as np
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/"
spec = importlib.util.spec_from_file_location('r', R + 'scripts/iem_disease_ranking.py'); RK = importlib.util.module_from_spec(spec); spec.loader.exec_module(RK)
sw = json.load(open(R + 'results/wbm_iem/pin_sweep/Harvey_1_03d_pin_sweep_development_v1.json'))
m = {r['iem']: r for r in json.load(open(R + 'results/wbm_iem/Harvey_1_03d_iem_cross_ranking_v1.json'))}
P = json.load(open(R + 'data/iem/iem_ranking_profiles_v0.2.json'))['profiles']
dev = sorted(sw['provenance']['iems'], key=lambda i: m[i]['call_index'])
TOL = 1e-6
def c(h, d):
    x = d - h
    return 1 if x > TOL else -1 if x < -TOL else 0
H = {}   # (iem, rid, alpha) -> healthy value
D = {}
for rid, per in sw['readouts'].items():
    for iem in dev:
        e = next(x for x in m[iem]['readouts'] if x['reaction'] == rid)
        if e['status_disease'] != 'Optimal':
            continue
        D[(iem, rid)] = e['disease']
        for x in per[iem]:
            if x['status'] == 'Optimal':
                H[(iem, rid, float(x['alpha']))] = x['value']
panel = list(sw['readouts'])
def rule_calls(rule):
    calls = {i: {} for i in dev}
    for (iem, rid), d in D.items():
        hs = {a: H.get((iem, rid, a)) for a in (1.0, 0.5, 0.1, 0.01, 0.001, 0.0)}
        if any(v is None for v in hs.values()):
            continue
        if rule == 'a1':
            v = c(hs[1.0], d)
        elif rule == 'a0.5':
            v = c(hs[0.5], d)
        elif rule == 'coupled':        # (a) increases only where coupled at a small pin; decreases from the same pin
            v = c(hs[0.001], d)
        elif rule == 'capped_indet':   # (b) alpha 1 calls; an increase whose healthy maximum reaches the disease maximum at alpha 0.5 is indeterminate (0)
            v = c(hs[1.0], d)
            if v == 1 and abs(hs[0.5] - d) <= TOL:
                v = 0
        calls[iem][rid] = v
    return calls
na = {i: 0 for i in dev}
res = {}
for rule in ('a1', 'a0.5', 'coupled', 'capped_indet'):
    calls = rule_calls(rule)
    res[rule] = {}
    for pset in ('lab', 'hpo'):
        prof = {i: p for i, p in P[pset].items() if i in dev}
        for adj in (False, True):
            names, S = RK.score_matrix(dev, calls, prof, adj, panel)
            rr = {n: RK.rank_stats(S[a], dev.index(n))['expected_reciprocal_rank'] for a, n in enumerate(names)}
            res[rule][(pset, adj)] = rr
    own = [calls[i].get(r) == RK.SIGN[s] for i in dev for r, s in P['lab'][i] if r in calls[i]]
    inc_share = np.median([sum(v == 1 for v in calls[i].values()) / len(calls[i]) for i in dev])
    print(f"{rule:13s} own {sum(own)}/{len(own)}  inc share {inc_share:.2f}  " + '  '.join(
        f"{p}{'-adj' if a else ''} {np.mean(list(res[rule][(p, a)].values())):.3f}" for p in ('lab', 'hpo') for a in (False, True)))
# paired comparison alpha 0.5 vs 1 (and rules vs 1), lab plain and adjusted
from scipy.stats import wilcoxon
rng = np.random.default_rng(20261008)
for rule in ('a0.5', 'coupled', 'capped_indet'):
    for pset in ('lab', 'hpo'):
        for adj in (False, True):
            a, b = res['a1'][(pset, adj)], res[rule][(pset, adj)]
            ks = sorted(a)
            diff = np.array([b[k] - a[k] for k in ks])
            boot = [np.mean(rng.choice(diff, len(diff))) for _ in range(10000)]
            better, worse = int((diff > 1e-9).sum()), int((diff < -1e-9).sum())
            try:
                p = wilcoxon(diff[np.abs(diff) > 1e-9]).pvalue
            except ValueError:
                p = float('nan')
            print(f"  {rule:13s} vs alpha 1, {pset}{'-adj' if adj else '    '}: mean diff {diff.mean():+.3f} [95% boot {np.percentile(boot,2.5):+.3f}, {np.percentile(boot,97.5):+.3f}]  better {better} worse {worse}  Wilcoxon p {p:.3f}")
print('\nPer-profile expected reciprocal rank, lab plain (alpha 1 -> 0.5):')
for k in sorted(res['a1'][('lab', False)]):
    print(f"  {k:7s} vmax {m[k]['vmax_healthy']:9.0f}: {res['a1'][('lab', False)][k]:.2f} -> {res['a0.5'][('lab', False)][k]:.2f}")
