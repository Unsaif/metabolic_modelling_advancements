"""Task C: independent recomputation of every ranking metric (no project code imported).

Calls are rebuilt from the stored healthy/disease values. Adjusted scores use exact fractions so ties are exact.
"""
import json
import math
import statistics
import sys
from fractions import Fraction
from collections import Counter, defaultdict

import numpy as np

ROOT = "/home/claude/mma"
RES = f"{ROOT}/results/wbm_iem"
load = lambda p: json.load(open(p))
M = load(f"{RES}/Harvey_1_03d_iem_cross_ranking_v1.json")
P = load(f"{ROOT}/data/iem/iem_ranking_profiles_v0.2.json")
R = load(f"{RES}/ranking/Harvey_1_03d_ranking_v1.json")

cands = [c["iem"] for c in sorted(P["candidates"], key=lambda c: c["call_index"])]
assert cands == [r["iem"] for r in M], "candidate order"
N = len(cands)
TOL = 1e-6


def call_protocol(h, d, ok):
    if not ok:
        return None
    diff = d - h
    return 1 if diff > TOL else -1 if diff < -TOL else 0


def call_material(h, d, ok, rel_den="max"):
    if not ok:
        return None
    diff = d - h
    if abs(diff) <= 1e-3:
        return 0
    den = max(abs(h), abs(d)) if rel_den == "max" else abs(h)
    rel = diff / den if den > 0 else math.copysign(math.inf, diff)
    if rel > 0.05:
        return 1
    if rel < -0.05:
        return -1
    return 0


def build_calls(rule, **kw):
    C = {}
    for r in M:
        row = {}
        for e in r["readouts"]:
            ok = e["status_healthy"] == "Optimal" and e["status_disease"] == "Optimal"
            h, d = e["healthy"], e["disease"]
            row[e["reaction"]] = (call_protocol(h, d, ok) if rule == "protocol" else call_material(h, d, ok, **kw))
        C[r["iem"]] = row
    return C


def promiscuity(C):
    pu, pd = {}, {}
    for c in cands:
        avail = [v for v in C[c].values() if v is not None]
        pu[c] = Fraction(sum(v == 1 for v in avail), len(avail))
        pd[c] = Fraction(sum(v == -1 for v in avail), len(avail))
    return pu, pd


def H(n):
    return sum(Fraction(1, i) for i in range(1, n + 1))


def rank_stats(a, t):
    er = a + Fraction(t + 1, 2)
    err = sum(Fraction(1, i) for i in range(a + 1, a + t + 1)) / t
    p1 = min(max(Fraction(1 - a, t), 0), 1)
    p5 = min(max(Fraction(5 - a, t), 0), 1)
    return er, err, p1, p5


def analyse(profiles, C, score_kind):
    pu, pd = promiscuity(C)
    out = {}
    for T, tuples in profiles.items():
        n_up = sum(1 for _, s in tuples if s == "Increased")
        n_dn = sum(1 for _, s in tuples if s == "Decreased")
        scores = {}
        for c in cands:
            S = 0
            for rid, s in tuples:
                v = C[c].get(rid, "missing")
                if v == "missing":
                    raise KeyError(rid)
                v = v or 0
                S += (1 if s == "Increased" else -1) * v
            S = Fraction(S)
            if score_kind == "adjusted":
                S = S - (n_up - n_dn) * (pu[c] - pd[c])
            scores[c] = S
        sT = scores[T]
        a = sum(1 for c in cands if scores[c] > sT)
        t = sum(1 for c in cands if scores[c] == sT)
        er, err, p1, p5 = rank_stats(a, t)
        # every candidate's stats as if it were the true disease (for the Monte Carlo null)
        allc = []
        for c in cands:
            ac = sum(1 for x in cands if scores[x] > scores[c]); tc = sum(1 for x in cands if scores[x] == scores[c])
            allc.append(rank_stats(ac, tc)[1:])
        own = sum(1 for rid, s in tuples if (C[T][rid] or 0) == (1 if s == "Increased" else -1))
        out[T] = dict(size=len(tuples), n_up=n_up, n_dn=n_dn, score=sT, a=a, t=t, er=er, err=err, p1=p1, p5=p5,
                      higher=[c for c in cands if scores[c] > sT], tied=[c for c in cands if scores[c] == sT and c != T],
                      own=own, allc=allc, scores=scores)
    return out


def summarize(res, rng, draws):
    names = list(res)
    n = len(names)
    mrr = sum(res[k]["err"] for k in names) / n
    top1 = sum(res[k]["p1"] for k in names)
    top5 = sum(res[k]["p5"] for k in names)
    med = statistics.median([float(res[k]["er"]) for k in names])
    # Monte Carlo null: replace the true disease by a uniformly random candidate in each profile
    E = np.array([[float(x[0]) for x in res[k]["allc"]] for k in names])   # n x 57
    P1 = np.array([[float(x[1]) for x in res[k]["allc"]] for k in names])
    P5 = np.array([[float(x[2]) for x in res[k]["allc"]] for k in names])
    idx = rng.integers(0, N, size=(draws, n))
    rows = np.arange(n)[None, :]
    null_mrr = E[rows, idx].mean(axis=1)
    null_t1 = P1[rows, idx].sum(axis=1)
    null_t5 = P5[rows, idx].sum(axis=1)
    eps = 1e-12
    k_mrr = int((null_mrr >= float(mrr) - eps).sum())
    k_t1 = int((null_t1 >= float(top1) - eps).sum())
    k_t5 = int((null_t5 >= float(top5) - eps).sum())
    return dict(n=n, mrr=float(mrr), top1=float(top1), top5=float(top5), median_rank=med,
                mean_rank=float(sum(res[k]["er"] for k in names) / n),
                n_no_higher=sum(1 for k in names if res[k]["a"] == 0),
                n_unique_first=sum(1 for k in names if res[k]["a"] == 0 and res[k]["t"] == 1),
                null_mrr=float(H(N) / N), null_top1=n / N, null_top5=5 * n / N,
                mc_k=(k_mrr, k_t1, k_t5), mc_p=tuple((1 + k) / (1 + draws) for k in (k_mrr, k_t1, k_t5)),
                mc_null_max=(float(null_mrr.max()), float(null_t1.max()), float(null_t5.max())),
                mc_null_mean=(float(null_mrr.mean()), float(null_t1.mean()), float(null_t5.mean())))


def strata(res):
    bins = {"1": lambda s: s == 1, "2-3": lambda s: 2 <= s <= 3, "4-6": lambda s: 4 <= s <= 6, ">=7": lambda s: s >= 7}
    out = {}
    for b, f in bins.items():
        ks = [k for k in res if f(res[k]["size"])]
        if not ks:
            out[b] = dict(n=0, mrr=float("nan"), top1=0.0, top5=0.0, median=float("nan"))
            continue
        out[b] = dict(n=len(ks), mrr=float(sum(res[k]["err"] for k in ks) / len(ks)), top1=float(sum(res[k]["p1"] for k in ks)),
                      top5=float(sum(res[k]["p5"] for k in ks)), median=statistics.median([float(res[k]["er"]) for k in ks]))
    return out


if __name__ == "__main__":
    draws = 100_000
    rng = np.random.default_rng(987654321)   # my own seed, not the project's
    calls = {"protocol": build_calls("protocol"), "material": build_calls("material")}
    # sanity: protocol calls equal the stored predicted calls
    code = {"Increased": 1, "Decreased": -1, "Unchanged": 0, "NA": None}
    mism = sum(1 for r in M for e in r["readouts"] if code[e["predicted"]] != calls["protocol"][r["iem"]][e["reaction"]])
    print("protocol calls vs stored predicted: mismatches", mism)
    # profile readouts all in panel?
    panel = {e["reaction"] for e in M[0]["readouts"]}
    for s, prof in P["profiles"].items():
        miss = sorted({rid for t in prof.values() for rid, _ in t if rid not in panel})
        dup = [k for k, t in prof.items() if len({rid for rid, _ in t}) != len(t)]
        print(f"profile set {s}: {len(prof)} profiles, {sum(len(t) for t in prof.values())} tuples; readouts not in panel {miss}; "
              f"profiles with repeated readout {dup}; non-candidates {[k for k in prof if k not in cands]}")
    allres = {}
    table = []
    for rule in ("protocol", "material"):
        for kind in ("plain", "adjusted"):
            for pset in ("lab", "hpo", "hpo_frequent", "lab_hpo_corroborated"):
                key = f"{pset}|{rule}|{kind}"
                res = analyse(P["profiles"][pset], calls[rule], kind)
                allres[key] = res
                sm = summarize(res, rng, draws)
                st = strata(res)
                # compare with the project's output
                ps = R["analyses"][key]["summary"]
                pp = {x["iem"]: x for x in R["analyses"][key]["profiles"]}
                dif = []
                for f_me, f_them in (("mrr", "mrr"), ("top1", "expected_top1"), ("top5", "expected_top5"),
                                     ("median_rank", "median_expected_rank"), ("mean_rank", "mean_expected_rank"),
                                     ("n_no_higher", "n_true_disease_among_top_scores"), ("null_mrr", "null_mrr"),
                                     ("null_top1", "null_top1"), ("null_top5", "null_top5")):
                    if abs(sm[f_me] - ps[f_them]) > 1e-9:
                        dif.append((f_me, sm[f_me], ps[f_them]))
                prof_dif = 0
                for T, x in res.items():
                    y = pp[T]
                    if (x["a"] != y["higher"] or x["t"] != y["tied"] or abs(float(x["er"]) - y["expected_rank"]) > 1e-9
                            or abs(float(x["err"]) - y["expected_reciprocal_rank"]) > 1e-9 or abs(float(x["p1"]) - y["p_top1"]) > 1e-9
                            or abs(float(x["p5"]) - y["p_top5"]) > 1e-9 or abs(float(x["score"]) - y["score"]) > 1e-9
                            or sorted(x["higher"]) != sorted(y["candidates_scoring_higher"])
                            or sorted(x["tied"]) != sorted(y["candidates_tied"]) or x["own"] != y["own_matches"]):
                        prof_dif += 1
                        if key == "lab|protocol|plain" and prof_dif <= 3:
                            print("   DIFF", T, dict(a=x["a"], t=x["t"], er=float(x["er"]), own=x["own"], score=float(x["score"]),
                                  nh=len(x["higher"]), nt=len(x["tied"])),
                                  {k: y[k] for k in ("higher", "tied", "expected_rank", "own_matches", "score")},
                                  len(y["candidates_scoring_higher"]), len(y["candidates_tied"]))
                st_them = ps.get("strata_by_profile_size", {})
                st_dif = []
                for b, v in st.items():
                    w = st_them.get(b if b != ">=7" else next((k for k in st_them if "7" in k), b))
                    if v["n"] == 0:
                        if w not in (None, {}) and w.get("n_profiles", 0) != 0:
                            st_dif.append((b, v, w))
                        continue
                    if w is None or w["n_profiles"] != v["n"] or abs(w["mrr"] - v["mrr"]) > 1e-9 or \
                            abs(w["expected_top1"] - v["top1"]) > 1e-9 or abs(w["expected_top5"] - v["top5"]) > 1e-9 or \
                            abs(w["median_expected_rank"] - v["median"]) > 1e-9:
                        st_dif.append((b, v, w))
                print(f"{key:40s} n={sm['n']:2d} MRR={sm['mrr']:.4f} top1={sm['top1']:.2f} top5={sm['top5']:.2f} "
                      f"med={sm['median_rank']:.1f} mean={sm['mean_rank']:.2f} nohigher={sm['n_no_higher']} uniq1={sm['n_unique_first']} "
                      f"MC k={sm['mc_k']} p={['%.2g' % p for p in sm['mc_p']]} nullmax={tuple(round(x, 3) for x in sm['mc_null_max'])} "
                      f"| project p=({ps['p_mrr']:.2g},{ps['p_top1']:.2g},{ps['p_top5']:.2g}) summary diffs {dif} profile diffs {prof_dif} strata diffs {len(st_dif)}")
                table.append((key, sm, st))
    # strata for the primary
    print("\nstrata, primary:")
    for k, v in table[0][2].items():
        print(f"  {k:4s} n={v['n']:2d} MRR={v['mrr']:.3f} top1={v['top1']:.2f} top5={v['top5']:.2f} median={v['median']:.1f}")
    print("strata keys in project output:", list(R["analyses"]["lab|protocol|plain"]["summary"]["strata_by_profile_size"].keys()))
    print("exact null MRR H57/57 =", float(H(57) / 57))
    import pickle
    with open("c_ranking_results.pkl", "wb") as fh:
        pickle.dump({"allres": allres, "table": table, "calls": calls}, fh)
