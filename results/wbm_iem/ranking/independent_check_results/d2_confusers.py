import json, pickle
from collections import Counter
R = json.load(open("/home/claude/mma/results/wbm_iem/ranking/Harvey_1_03d_ranking_v1.json"))
D = pickle.load(open("c_ranking_results.pkl", "rb"))
K = "confusers_[candidate, profiles_strictly_higher, profiles_tied_with_positive_true_score]"
bad = 0
for key, res in D["allres"].items():
    higher = Counter(); tiepos = Counter()
    for T, x in res.items():
        for c in x["higher"]: higher[c] += 1
        for c in x["tied"]:
            if x["score"] > 0: tiepos[c] += 1
    theirs = {c: (h, t) for c, h, t in R["analyses"][key]["summary"][K]}
    mine = {c: (higher[c], tiepos[c]) for c in theirs}
    allc = set(higher) | set(tiepos)
    missing = [c for c in allc if c not in theirs and (higher[c] or tiepos[c])]
    d = {c: (mine[c], theirs[c]) for c in theirs if mine[c] != theirs[c]}
    if d or missing:
        bad += 1; print(key, "differs", list(d.items())[:5], "missing", missing)
print("analyses with confuser-count differences:", bad)
# adjusted, lab: top list with promiscuity shares
calls = D["calls"]["protocol"]
cands = list(calls)
share = {c: sum(v == 1 for v in calls[c].values() if v is not None) / sum(v is not None for v in calls[c].values()) for c in cands}
med = sorted(share.values())[len(share)//2]
print("median increase share", round(med, 3))
theirs = R["analyses"]["lab|protocol|adjusted"]["summary"][K]
for c, h, t in theirs[:15]:
    print(f"  {c:7s} higher {h:2d} tiedpos {t}  share Increased {share[c]:.2f}")
