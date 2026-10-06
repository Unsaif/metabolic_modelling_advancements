"""Independent implementation of the ranking rules (plan section 'Scoring and ranking' / 'Metrics and test'),
compared with scripts/iem_disease_ranking.py on synthetic matrices in the runner's output schema."""
import importlib.util, math, random, statistics
from fractions import Fraction
import numpy as np
spec = importlib.util.spec_from_file_location("R", "/home/claude/mma/scripts/iem_disease_ranking.py")
R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)

# ---------------- my implementation -----------------
def my_call(x, rule):
    if rule == "protocol":
        return {"Increased": 1, "Decreased": -1, "Unchanged": 0}.get(x["predicted"])  # NA/absent -> None
    if x["status_healthy"] != "Optimal" or x["status_disease"] != "Optimal":
        return None
    h, d = x["healthy"], x["disease"]
    if h is None or d is None: return None
    delta = d - h
    m = max(abs(h), abs(d))
    if abs(delta) <= 1e-3 or m == 0: return 0
    if abs(delta) / m <= 0.05: return 0
    return 1 if delta > 0 else -1

def my_scores(matrix, profiles, rule, adjusted, panel):
    cands = [r["iem"] for r in sorted(matrix, key=lambda r: r["call_index"])]
    calls = {r["iem"]: {x["reaction"]: my_call(x, rule) for x in r["readouts"]} for r in matrix}
    out = {}
    for true, prof in profiles.items():
        if not prof: continue
        nu = sum(s == "Increased" for _, s in prof); nd = len(prof) - nu
        row = []
        for d in cands:
            c = calls[d]
            S = Fraction(sum((1 if s == "Increased" else -1) * (c.get(r) or 0) for r, s in prof))
            if adjusted:
                av = [c[r] for r in panel if c.get(r) is not None]
                if av:
                    S -= (nu - nd) * Fraction(sum(v == 1 for v in av) - sum(v == -1 for v in av), len(av))
            row.append(S)
        out[true] = row
    return cands, out

def my_rank(row, j):
    s = row[j]
    a = sum(x > s for x in row); t = sum(x == s for x in row)
    err = Fraction(sum(Fraction(1, i) for i in range(a + 1, a + t + 1)), t)
    clip = lambda v: min(Fraction(1), max(Fraction(0), v))
    return a, t, Fraction(2 * a + t + 1, 2), err, clip(Fraction(1 - a, t)), clip(Fraction(5 - a, t))

def compare(matrix, profiles, panel, label, draws=4000):
    worst = 0.0
    for rule in ("protocol", "material"):
        order, calls, n_na = R.calls_from_matrix(matrix, rule)
        for adjusted in (False, True):
            res = R.analyse(order, calls, n_na, profiles, adjusted, panel, draws=draws)
            cands, mine = my_scores(matrix, profiles, rule, adjusted, panel)
            assert order == cands
            per = {p["iem"]: p for p in res["profiles"]}
            assert set(per) == set(mine), (set(per) ^ set(mine))
            rr, t1, t5, ranks = [], [], [], []
            for true, row in mine.items():
                a, t, er, err, p1, p5 = my_rank(row, cands.index(true))
                p = per[true]
                assert (p["higher"], p["tied"]) == (a, t), (label, rule, adjusted, true, (p["higher"], p["tied"]), (a, t))
                for k, v in (("expected_rank", er), ("expected_reciprocal_rank", err), ("p_top1", p1), ("p_top5", p5)):
                    worst = max(worst, abs(p[k] - float(v)))
                assert abs(p["score"] - float(row[cands.index(true)])) < 1e-9
                ge = sorted(d for d, x in zip(cands, row) if x >= row[cands.index(true)] and d != true)
                assert p["n_candidates_scoring_at_least_as_high"] == len(ge)
                assert sorted(x[0] for x in p["candidates_scoring_at_least_as_high"]) == ge[:10] or len(ge) > 10
                rr.append(err); t1.append(p1); t5.append(p5); ranks.append(er)
            s = res["summary"]; n = len(cands)
            assert abs(s["mrr"] - float(sum(rr) / len(rr))) < 1e-12
            assert abs(s["expected_top1"] - float(sum(t1))) < 1e-9 and abs(s["expected_top5"] - float(sum(t5))) < 1e-9
            assert abs(s["median_expected_rank"] - float(statistics.median(ranks))) < 1e-12
            assert abs(s["null_mrr"] - sum(1 / i for i in range(1, n + 1)) / n) < 1e-15
            assert s["n_true_disease_among_top_scores"] == sum(my_rank(r, cands.index(tr))[0] == 0 for tr, r in mine.items())
    return worst

def synth(n=57, nr=177, seed=1, p_na=0.03, n_profiles=57, all_na_readout=True):
    rng = random.Random(seed)
    readouts = [f"R{i}" for i in range(nr)]
    matrix = []
    for k in range(n):
        rows = []
        # candidate-specific promiscuity
        pu, pdn = rng.uniform(0, 0.4), rng.uniform(0, 0.2)
        for r in readouts:
            if r == "R0" and all_na_readout:      # absent from model for everyone (like EX_25aics[u])
                rows.append({"reaction": r, "predicted": "NA", "healthy": None, "disease": None,
                             "status_healthy": "absent", "status_disease": "absent"}); continue
            if rng.random() < p_na:
                rows.append({"reaction": r, "predicted": "NA", "healthy": None, "disease": None,
                             "status_healthy": "Optimal", "status_disease": "Infeasible"}); continue
            h = rng.choice([0.0, rng.uniform(0, 2), rng.uniform(0, 2000)])
            u = rng.random()
            if u < pu: d = h + rng.choice([1.2e-6, 5e-4, 2e-3, h * 0.03 + 2e-3, h * 0.5 + 1])
            elif u < pu + pdn: d = max(0.0, h - rng.choice([1.2e-6, 5e-4, 2e-3, h * 0.03 + 2e-3, h * 0.5]))
            else: d = h
            diff = d - h
            pred = "Increased" if diff > 1e-6 else "Decreased" if diff < -1e-6 else "Unchanged"
            rows.append({"reaction": r, "predicted": pred, "healthy": h, "disease": d,
                         "status_healthy": "Optimal", "status_disease": "Optimal"})
        matrix.append({"iem": f"D{k}", "call_index": k + 1, "readouts": rows})
    profiles = {}
    for k in range(n_profiles):
        size = rng.choice([1, 1, 2, 3, 4, 5, 6, 8, 12])
        rs = rng.sample(readouts, min(size, nr))
        profiles[f"D{k}"] = [[r, rng.choice(["Increased", "Increased", "Increased", "Decreased"])] for r in rs]
    profiles["D3"] = [["R0", "Increased"]]                          # profile NA for every candidate
    return matrix, profiles, readouts

worst = 0
for seed in range(1, 31):
    m, prof, panel = synth(seed=seed)
    worst = max(worst, compare(m, prof, panel, f"seed{seed}", draws=500))
# small-n case with many ties
for seed in range(31, 61):
    m, prof, panel = synth(n=8, nr=6, seed=seed, n_profiles=8)
    worst = max(worst, compare(m, prof, panel, f"small{seed}", draws=500))
print("all per-profile ranks, scores, summaries agree; max |float - exact| =", worst)
