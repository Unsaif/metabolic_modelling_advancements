"""Task D: every number and claim in the results write-up, recomputed independently."""
import json
import pickle
import statistics
from collections import Counter, defaultdict

ROOT = "/home/claude/mma"
RES = f"{ROOT}/results/wbm_iem"
load = lambda p: json.load(open(p))
M = load(f"{RES}/Harvey_1_03d_iem_cross_ranking_v1.json")
P = load(f"{ROOT}/data/iem/iem_ranking_profiles_v0.2.json")
D = pickle.load(open("c_ranking_results.pkl", "rb"))
allres, calls = D["allres"], D["calls"]
cands = [r["iem"] for r in M]

print("== precise values for rounding")
for key, sm, st in D["table"]:
    print(f"  {key:40s} MRR {sm['mrr']:.5f} top1 {sm['top1']:.4f} top5 {sm['top5']:.4f} med {sm['median_rank']}")
    if key == "lab|protocol|plain":
        for b, v in st.items():
            print(f"     stratum {b}: n {v['n']} MRR {v['mrr']:.4f} top1 {v['top1']:.4f} top5 {v['top5']:.4f} median {v['median']}")

prim = allres["lab|protocol|plain"]
adj = allres["lab|protocol|adjusted"]
print("adjusted/plain MRR ratio", D["table"][4][1]["mrr"] / D["table"][0][1]["mrr"])

print("\n== per-profile table (lab): IEM size(up/down) own plain adjusted")
wr = {}
for line in open(f"{ROOT}/docs/studies/wbm-iem-ranking-results.md"):
    if line.startswith("| ") and "(" in line and "/" in line and line.count("|") == 6 and not line.startswith("| IEM"):
        parts = [x.strip() for x in line.strip().strip("|").split("|")]
        wr[parts[0]] = parts
bad = []
for T in cands:
    x, y = prim[T], adj[T]
    mine = [T, f"{x['size']} ({x['n_up']}/{x['n_dn']})", str(x["own"]), f"{float(x['er']):.1f}", f"{float(y['er']):.1f}"]
    if wr.get(T) != mine:
        bad.append((wr.get(T), mine))
print("rows in write-up:", len(wr), "rows differing:", bad)

print("\n== decreases and unique first (plain, lab)")
with_dec = [T for T in cands if prim[T]["n_dn"] > 0]
all_inc = [T for T in cands if prim[T]["n_dn"] == 0]
uniq = [T for T in cands if prim[T]["a"] == 0 and prim[T]["t"] == 1]
print("profiles with a decrease", len(with_dec), "all-increase", len(all_inc))
print("uniquely first (plain):", uniq, "with decrease:", sum(T in with_dec for T in uniq), "all-increase:", sum(T in all_inc for T in uniq))
print("no candidate higher (plain):", sum(prim[T]["a"] == 0 for T in cands))
ties_allinc = [prim[T]["t"] - 1 for T in all_inc]
print("all-increase profiles: median number of others tied", statistics.median(ties_allinc),
      "| among those with none higher:", statistics.median([prim[T]["t"] - 1 for T in all_inc if prim[T]["a"] == 0]),
      "n none higher", sum(prim[T]["a"] == 0 for T in all_inc))
uniq_adj = [T for T in cands if adj[T]["a"] == 0 and adj[T]["t"] == 1]
print("uniquely first (adjusted):", len(uniq_adj), "all-increase among them:", sum(T in all_inc for T in uniq_adj),
      sorted(T for T in uniq_adj if T in all_inc))
print("with-decrease uniquely first under adjusted:", sorted(T for T in uniq_adj if T in with_dec))
print("ranked below most (plain expected rank > 29):", sorted(((float(prim[T]['er']), T, prim[T]['own'], prim[T]['size']) for T in cands if prim[T]['er'] > 25), reverse=True))

print("\n== confusions")
for key in ("lab|protocol|plain", "lab|protocol|adjusted"):
    res = allres[key]
    higher = Counter(); tiepos = Counter(); tieany = Counter()
    for T, x in res.items():
        for c in x["higher"]:
            higher[c] += 1
        for c in x["tied"]:
            tieany[c] += 1
            if x["score"] > 0:
                tiepos[c] += 1
    hv = [higher[c] for c in cands]; tv = [tiepos[c] for c in cands]
    print(key, "outranks: max", max(hv), "median", statistics.median(hv), "| ties with positive true: median", statistics.median(tv),
          "range", min(tv), max(tv))
    print("   top outrankers:", higher.most_common(14))
    print("   ties (positive) for top outrankers:", {c: tiepos[c] for c, _ in higher.most_common(14)})
    print("   ties (any) for top outrankers:", {c: tieany[c] for c, _ in higher.most_common(14)})

print("\n== generic increases (protocol rule)")
for rule in ("protocol", "material"):
    C = calls[rule]
    up_share = []; dn_share = []
    for c in cands:
        avail = [v for v in C[c].values() if v is not None]
        up_share.append(sum(v == 1 for v in avail) / len(avail)); dn_share.append(sum(v == -1 for v in avail) / len(avail))
    print(rule, "n available", len(avail), "median share Increased", round(100 * statistics.median(up_share), 2),
          "range", round(100 * min(up_share), 2), round(100 * max(up_share), 2),
          "median share Decreased", round(100 * statistics.median(dn_share), 2))
    rx = [e["reaction"] for e in M[0]["readouts"] if e["status_healthy"] != "absent"]
    n_up_by_r = [sum(C[c][r] == 1 for c in cands) for r in rx]
    print("   readouts", len(rx), "median knockouts calling a readout Increased", statistics.median(n_up_by_r),
          "readouts increased by >= 40 knockouts", sum(n >= 40 for n in n_up_by_r))
    low = sorted((round(100 * s, 1), c) for s, c in zip(up_share, cands))[:5]
    print("   lowest increase shares", low)

print("\n== matching known directions across the lab tuples (protocol rule)")
C = calls["protocol"]
own_ok = own_n = oth_ok = oth_n = 0
oth_ok_inc = oth_n_inc = oth_ok_dec = oth_n_dec = 0
for T, tuples in P["profiles"]["lab"].items():
    for rid, s in tuples:
        sv = 1 if s == "Increased" else -1
        if C[T][rid] is None:
            continue
        own_n += 1; own_ok += C[T][rid] == sv
        for c in cands:
            if c == T or C[c][rid] is None:
                continue
            oth_n += 1; oth_ok += C[c][rid] == sv
            if sv == 1:
                oth_n_inc += 1; oth_ok_inc += C[c][rid] == sv
            else:
                oth_n_dec += 1; oth_ok_dec += C[c][rid] == sv
print(f"true disease {own_ok}/{own_n} = {100*own_ok/own_n:.2f}%; other knockouts {oth_ok}/{oth_n} = {100*oth_ok/oth_n:.2f}%")
print(f"   others on known increases {100*oth_ok_inc/oth_n_inc:.1f}% ; on known decreases {100*oth_ok_dec/oth_n_dec:.1f}%")

print("\n== HPO limits")
hpo = P["profiles"]["hpo"]
sizes = Counter(len(t) for t in hpo.values())
print("hpo profile sizes", sorted(sizes.items()), "profiles with 1-2 tuples", sum(v for k, v in sizes.items() if k <= 2))
ann = P["hpo_annotations"]
used_terms = defaultdict(set)
for a in ann:
    if a.get("used"):
        used_terms[a["iem"]].add(a["hpo_id"])
nt = Counter(len(v) for k, v in used_terms.items() if k in hpo)
print("distinct used terms per hpo profile", sorted(nt.items()), "with 1-2 terms", sum(v for k, v in nt.items() if k <= 2),
      "hpo IEMs without used terms", [k for k in hpo if k not in used_terms])
amm = [k for k, t in hpo.items() if any("nh4" in rid for rid, _ in t)]
amm_terms = sorted({(a["hpo_id"], a["hpo_term"], a["readout"]) for a in ann if a.get("used") and "mmon" in a["hpo_term"]})
print("hpo profiles with an ammonia readout", len(amm), amm)
print("hyperammonaemia-type terms used", amm_terms)
amm_iems = sorted({a["iem"] for a in ann if a.get("used") and "mmon" in a["hpo_term"]})
print("IEMs with a used ammonia term", len(amm_iems), amm_iems)
