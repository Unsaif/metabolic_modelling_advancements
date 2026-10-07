import pickle, statistics
from collections import defaultdict
D = pickle.load(open("c_ranking_results.pkl", "rb"))
calls = D["calls"]["protocol"]; cands = list(calls)
res = D["allres"]["lab|protocol|adjusted"]; plain = D["allres"]["lab|protocol|plain"]
def shares(c):
    av = [v for v in calls[c].values() if v is not None]
    return sum(v == 1 for v in av) / len(av), sum(v == -1 for v in av) / len(av), sum(v != 0 for v in av)
net = {c: shares(c)[0] - shares(c)[1] for c in cands}
nchg = {c: shares(c)[2] for c in cands}
print("median net (p_up - p_down)", round(statistics.median(net.values()), 3), "median n changed", statistics.median(nchg.values()))
where = defaultdict(list)
for T, x in res.items():
    for c in x["higher"]:
        where[c].append(T)
for c in ["AADC", "BTD", "SUCLA", "TETB", "AGAT", "GMT", "LTC4S", "PC", "ARG", "CPS1", "OTC", "CIT1", "FED"]:
    ts = where[c]
    kinds = sum(1 for T in ts if res[T]["n_dn"] == 0)
    print(f"{c:6s} net {net[c]:.2f} changed {nchg[c]:3d}/186 outranks in {len(ts):2d} profiles ({kinds} all-increase); "
          f"plain score of {c} in those: {sorted(float(plain[T]['scores'][c]) for T in ts)[:6]}... true plain scores {sorted(float(plain[T]['score']) for T in ts)[:6]}")
