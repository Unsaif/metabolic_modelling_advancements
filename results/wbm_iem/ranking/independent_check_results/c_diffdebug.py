import json, pickle
from fractions import Fraction
R = json.load(open("/home/claude/mma/results/wbm_iem/ranking/Harvey_1_03d_ranking_v1.json"))
D = pickle.load(open("c_ranking_results.pkl", "rb"))
from collections import Counter
fields = Counter()
for key, res in D["allres"].items():
    pp = {x["iem"]: x for x in R["analyses"][key]["profiles"]}
    for T, x in res.items():
        y = pp[T]
        checks = {"a": x["a"] == y["higher"], "t": x["t"] == y["tied"], "er": abs(float(x["er"]) - y["expected_rank"]) < 1e-9,
                  "err": abs(float(x["err"]) - y["expected_reciprocal_rank"]) < 1e-9, "p1": abs(float(x["p1"]) - y["p_top1"]) < 1e-9,
                  "p5": abs(float(x["p5"]) - y["p_top5"]) < 1e-9, "score": abs(float(x["score"]) - y["score"]) < 1e-9,
                  "higher_set": sorted(x["higher"]) == sorted(c for c, _ in y["candidates_scoring_higher"]),
                  "tied_set": sorted(x["tied"]) == sorted(y["candidates_tied"]), "own": x["own"] == y["own_matches"],
                  "size": x["size"] == y["profile_size"], "nup": x["n_up"] == y["n_increased"], "ndn": x["n_dn"] == y["n_decreased"]}
        for f, ok in checks.items():
            if not ok:
                fields[(key, f)] += 1
                if fields[(key, f)] == 1:
                    print(key, T, f, {k: x[k] for k in ("a", "t", "own", "size")}, float(x["score"]),
                          {k: y[k] for k in ("higher", "tied", "own_matches", "score", "profile_size")})
print(fields)
