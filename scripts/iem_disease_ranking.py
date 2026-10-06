"""Rank the simulated IEMs against each disease's biomarker profile (IEM disease-ranking study).

Usage:
  python scripts/iem_disease_ranking.py results/wbm_iem/Harvey_1_03d_iem_cross_ranking_v1.json
         [--profiles data/iem/iem_ranking_profiles_v0.2.json] [--out results/wbm_iem/ranking/Harvey_1_03d_ranking_v1.json]

The analysis follows docs/studies/wbm-iem-ranking-plan.md and is fixed with it:
  calls      c(d, r) in {+1, 0, -1}: the protocol's call for candidate IEM d and readout r (NA counts as 0);
             the "material" rule instead needs |disease - healthy| > 1e-3 and a relative change above 5%.
  score      S(d | P) = sum over the profile's (r, s) of s * c(d, r), s = +1 for Increased, -1 for Decreased;
             the "adjusted" score subtracts (n_up(P) - n_down(P)) * (p_up(d) - p_down(d)), where p_up and p_down
             are the fractions of d's available panel readouts called Increased and Decreased.
  rank       a = candidates scoring higher than the true disease, t = candidates tied with it (itself included);
             expected rank under random tie-breaking a + (t + 1) / 2, expected reciprocal rank
             (1/t) * sum_{i=a+1..a+t} 1/i, P(rank <= k) = clip((k - a) / t, 0, 1).
  test       one-sided Monte Carlo: each profile's true disease replaced by a uniformly drawn candidate;
             100,000 draws, seed 20261006; p = (1 + #(null >= observed)) / (1 + draws). The null mean of the MRR is
             exactly H_n / n for n candidates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from typing import Dict, List, Tuple

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROFILE_SETS = ["lab", "hpo", "hpo_frequent", "lab_hpo_corroborated"]
SIGN = {"Increased": 1, "Decreased": -1}
SEED = 20261006
DRAWS = 100_000


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def material_call(h, d, rel_tau=0.05, abs_floor=1e-3):
    if h is None or d is None or not (math.isfinite(h) and math.isfinite(d)):
        return None
    diff = d - h
    scale = max(abs(h), abs(d))
    rel = diff / scale if scale > 0 else 0.0
    if diff > abs_floor and rel > rel_tau:
        return 1
    if diff < -abs_floor and rel < -rel_tau:
        return -1
    return 0


def calls_from_matrix(results: List[dict], rule: str) -> Tuple[List[str], Dict[str, Dict[str, int]], Dict[str, int]]:
    """Candidate order (by call_index), calls {iem: {readout: +1/0/-1}} without NA, and NA counts per candidate."""
    results = sorted(results, key=lambda r: r["call_index"])
    order = [r["iem"] for r in results]
    calls, n_na = {}, {}
    for r in results:
        c, na = {}, 0
        for x in r.get("readouts", []):
            if rule == "protocol":
                v = {"Increased": 1, "Decreased": -1, "Unchanged": 0}.get(x["predicted"])
            elif rule == "material":
                ok = x["status_healthy"] == x["status_disease"] == "Optimal"
                v = material_call(x["healthy"], x["disease"]) if ok else None
            else:
                raise ValueError(rule)
            if v is None:
                na += 1
            else:
                c[x["reaction"]] = v
        calls[r["iem"]], n_na[r["iem"]] = c, na
    return order, calls, n_na


def score_matrix(order, calls, profiles: Dict[str, list], adjusted=False, panel=None):
    """Scores [profile x candidate] for the profiles whose true disease is a candidate."""
    names = [i for i in order if profiles.get(i)]
    S = np.zeros((len(names), len(order)))
    for a, iem in enumerate(names):
        prof = profiles[iem]
        n_up = sum(1 for _, s in prof if s == "Increased")
        n_down = sum(1 for _, s in prof if s == "Decreased")
        for b, cand in enumerate(order):
            c = calls[cand]
            S[a, b] = sum(SIGN[s] * c.get(r, 0) for r, s in prof)
            if adjusted:
                avail = [c[r] for r in (panel or c) if r in c]
                if avail:
                    p_up = sum(1 for v in avail if v == 1) / len(avail)
                    p_down = sum(1 for v in avail if v == -1) / len(avail)
                    S[a, b] -= (n_up - n_down) * (p_up - p_down)
    return names, S


def rank_stats(scores: np.ndarray, j: int, eps=1e-9):
    s = scores[j]
    a = int(np.sum(scores > s + eps))
    t = int(np.sum(np.abs(scores - s) <= eps))
    rr = sum(1.0 / i for i in range(a + 1, a + t + 1)) / t
    top = {k: min(1.0, max(0.0, (k - a) / t)) for k in (1, 5)}
    return {"higher": a, "tied": t, "expected_rank": a + (t + 1) / 2, "expected_reciprocal_rank": rr,
            "p_top1": top[1], "p_top5": top[5]}


def analyse(order, calls, n_na, profiles, adjusted=False, panel=None, draws=DRAWS, seed=SEED):
    names, S = score_matrix(order, calls, profiles, adjusted, panel)
    n = len(order)
    per, R = [], np.zeros((len(names), n, 3))
    for a, iem in enumerate(names):
        for b in range(n):
            st = rank_stats(S[a], b)
            R[a, b] = (st["expected_reciprocal_rank"], st["p_top1"], st["p_top5"])
        j = order.index(iem)
        st = rank_stats(S[a], j)
        higher = sorted(((order[b], float(S[a, b])) for b in range(n) if b != j and S[a, b] > S[a, j] + 1e-9),
                        key=lambda x: -x[1])
        tied = [order[b] for b in range(n) if b != j and abs(S[a, b] - S[a, j]) <= 1e-9]
        prof = profiles[iem]
        own = calls[iem]
        per.append({"iem": iem, "profile_size": len(prof),
                    "n_increased": sum(1 for _, s in prof if s == "Increased"),
                    "n_decreased": sum(1 for _, s in prof if s == "Decreased"),
                    "own_matches": sum(1 for r, s in prof if own.get(r, 0) == SIGN[s]),
                    "own_contradictions": sum(1 for r, s in prof if own.get(r, 0) == -SIGN[s]),
                    "score": float(S[a, j]), "best_score": float(S[a].max()),
                    **st, "candidates_scoring_higher": higher, "candidates_tied": tied})
    obs = {"mrr": float(np.mean([p["expected_reciprocal_rank"] for p in per])) if per else None,
           "top1": float(np.sum([p["p_top1"] for p in per])), "top5": float(np.sum([p["p_top5"] for p in per]))}
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(draws, len(names)))
    rows = np.arange(len(names))
    null_rr = R[rows, idx, 0].mean(axis=1) if len(names) else np.zeros(draws)
    null_top1 = R[rows, idx, 1].sum(axis=1) if len(names) else np.zeros(draws)
    null_top5 = R[rows, idx, 2].sum(axis=1) if len(names) else np.zeros(draws)

    def p(null, o):
        return float((1 + np.sum(null >= o - 1e-12)) / (1 + draws))
    h_n = sum(1.0 / i for i in range(1, n + 1))
    # Confusions (descriptive): how often each candidate scores strictly above a profile's true disease, and how
    # often it ties with a true disease whose score is positive (ties at zero or below carry no information).
    confusers = {c: {"strictly_higher": 0, "tied_positive": 0} for c in order}
    for x in per:
        for c, _ in x["candidates_scoring_higher"]:
            confusers[c]["strictly_higher"] += 1
        if x["score"] > 0:
            for c in x["candidates_tied"]:
                confusers[c]["tied_positive"] += 1
    confusers = sorted(([c, v["strictly_higher"], v["tied_positive"]] for c, v in confusers.items()
                        if v["strictly_higher"] or v["tied_positive"]), key=lambda r: (-r[1], -r[2], r[0]))
    strata = {}
    for lo, hi, label in ((1, 1, "1"), (2, 3, "2-3"), (4, 6, "4-6"), (7, 10 ** 6, ">=7")):
        sel = [x for x in per if lo <= x["profile_size"] <= hi]
        if sel:
            strata[label] = {"n_profiles": len(sel),
                             "mrr": float(np.mean([x["expected_reciprocal_rank"] for x in sel])),
                             "expected_top1": float(sum(x["p_top1"] for x in sel)),
                             "expected_top5": float(sum(x["p_top5"] for x in sel)),
                             "median_expected_rank": float(np.median([x["expected_rank"] for x in sel]))}
    summary = {"n_candidates": n, "n_profiles": len(per),
               "mrr": obs["mrr"], "null_mrr": h_n / n, "p_mrr": p(null_rr, obs["mrr"]) if per else None,
               "expected_top1": obs["top1"], "null_top1": len(per) / n, "p_top1": p(null_top1, obs["top1"]) if per else None,
               "expected_top5": obs["top5"], "null_top5": 5 * len(per) / n, "p_top5": p(null_top5, obs["top5"]) if per else None,
               "mean_expected_rank": float(np.mean([x["expected_rank"] for x in per])) if per else None,
               "median_expected_rank": float(np.median([x["expected_rank"] for x in per])) if per else None,
               "null_mean_rank": (n + 1) / 2,
               "n_true_disease_among_top_scores": sum(1 for x in per if x["higher"] == 0),
               "strata_by_profile_size": strata,
               "confusers_[candidate, profiles_strictly_higher, profiles_tied_with_positive_true_score]": confusers,
               "n_na_readouts_by_candidate": {k: v for k, v in n_na.items() if v}}
    return {"summary": summary, "profiles": per}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("matrix")
    ap.add_argument("--profiles", default=os.path.join(ROOT, "data", "iem", "iem_ranking_profiles_v0.2.json"))
    ap.add_argument("--out")
    ap.add_argument("--draws", type=int, default=DRAWS)
    args = ap.parse_args()
    with open(args.matrix) as fh:
        results = json.load(fh)
    with open(args.profiles) as fh:
        prof = json.load(fh)
    panel = prof["panel_protocol"] + prof["extra_readouts"]
    expected = [c["iem"] for c in sorted(prof["candidates"], key=lambda c: c["call_index"])]
    got = [r["iem"] for r in sorted(results, key=lambda r: r["call_index"])]
    if got != expected:
        raise SystemExit(f"matrix has {len(got)} IEMs; the plan needs all {len(expected)} candidates")
    readouts = {x["reaction"] for r in results for x in r["readouts"]}
    if set(panel) - readouts:
        raise SystemExit(f"matrix lacks panel readouts: {sorted(set(panel) - readouts)}")
    analyses = {}
    for rule in ("protocol", "material"):
        order, calls, n_na = calls_from_matrix(results, rule)
        for adjusted in (False, True):
            for ps in PROFILE_SETS:
                key = f"{ps}|{rule}|{'adjusted' if adjusted else 'plain'}"
                analyses[key] = analyse(order, calls, n_na, prof["profiles"][ps], adjusted, panel, draws=args.draws)
    out = {"plan": "docs/studies/wbm-iem-ranking-plan.md", "primary": "lab|protocol|plain",
           "matrix": os.path.relpath(os.path.abspath(args.matrix), ROOT), "matrix_sha256": sha256(args.matrix),
           "profiles": os.path.relpath(os.path.abspath(args.profiles), ROOT), "profiles_sha256": sha256(args.profiles),
           "seed": SEED, "draws": args.draws, "analyses": analyses}
    path = args.out or os.path.join(ROOT, "results", "wbm_iem", "ranking",
                                    os.path.basename(args.matrix).replace("_iem_cross_ranking", "_ranking")
                                    .replace("_iem_cross", "_ranking"))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    for key, a in analyses.items():
        s = a["summary"]
        if s["n_profiles"]:
            print(f"{key:42s} n={s['n_profiles']:2d} MRR={s['mrr']:.3f} (null {s['null_mrr']:.3f}, p={s['p_mrr']:.2g}) "
                  f"top1={s['expected_top1']:.1f} (null {s['null_top1']:.1f}) top5={s['expected_top5']:.1f} "
                  f"(null {s['null_top5']:.1f}) median rank {s['median_expected_rank']:.1f}")
    print("written", os.path.relpath(path, ROOT))


if __name__ == "__main__":
    main()
