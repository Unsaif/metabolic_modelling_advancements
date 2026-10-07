"""Task B: own protocol biomarkers in the matrix vs v0.3 (Harvey, HiGHS) and the baseline numbers (task D)."""
import json
from collections import Counter

ROOT = "/home/claude/mma"
RES = f"{ROOT}/results/wbm_iem"
load = lambda p: json.load(open(p))
merged = load(f"{RES}/Harvey_1_03d_iem_cross_ranking_v1.json")
protocol = load(f"{ROOT}/data/iem/iem_protocol_v0.2.json")
v03 = load(f"{RES}/Harvey_1_03d_iem_results_v0.3.json")
v04 = load(f"{RES}/Harvetta_1_03d_iem_results_v0.4.json")

ref = {(r["iem"], r["call_index"]): r for r in v03}
n = same = 0
diffs = []
calls_differ = []
max_abs = (0, None); max_rel = (0, None)
for r in merged:
    by = {e["reaction"]: e for e in r["readouts"]}
    rr = ref[(r["iem"], r["call_index"])]
    for b in rr["biomarkers"]:
        e = by[b["reaction"]]
        assert e["own"], (r["iem"], b["reaction"])
        n += 1
        if e["predicted"] == b["predicted"]:
            same += 1
        else:
            calls_differ.append((r["iem"], b["reaction"], b["predicted"], e["predicted"]))
        for s in ("healthy", "disease"):
            a, c = e[s], b[s]
            if a is None or c is None:
                continue
            # compare with v0.3 values zeroed by the same rule
            cz = 0.0 if abs(c) <= 1e-6 else c
            d = abs(a - cz)
            rel = d / max(1.0, abs(cz))
            if d > max_abs[0]:
                max_abs = (d, (r["iem"], b["reaction"], s, a, c))
            if rel > max_rel[0]:
                max_rel = (rel, (r["iem"], b["reaction"], s, a, c))
            diffs.append(d)
    # vmax
print("own biomarker calls compared", n, "same", same, "differ", calls_differ)
print("max abs diff", max_abs)
print("max rel diff (|d|/max(1,|ref|))", max_rel)
# relative using |ref| only (for values > 1e-3)
mr = 0; arg = None
for r in merged:
    by = {e["reaction"]: e for e in r["readouts"]}
    rr = ref[(r["iem"], r["call_index"])]
    for b in rr["biomarkers"]:
        e = by[b["reaction"]]
        for s in ("healthy", "disease"):
            a, c = e[s], b[s]
            if a is None or c is None or abs(c) < 1e-3:
                continue
            x = abs(a - c) / abs(c)
            if x > mr:
                mr, arg = x, (r["iem"], b["reaction"], s, a, c)
print("max rel diff (|d|/|ref|, ref >= 1e-3)", mr, arg)
vm = max(abs(r["vmax_healthy"] - ref[(r["iem"], r["call_index"])]["vmax_healthy"]) for r in merged)
print("max |vmax - v0.3 vmax|", vm)
pins = [(r["iem"], r["healthy_pin"], r["vmax_healthy"]) for r in merged]

# Protocol accuracy and the 'always increased' baseline, for v0.3 and v0.4 (and the new matrix)
def acc(results, label):
    tot = cor = inc = inc_c = dec = dec_c = na = 0
    for r in results:
        for b in r["biomarkers"]:
            if b["predicted"] == "NA":
                na += 1
                continue
            tot += 1
            ok = b["predicted"] == b["expected"]
            cor += ok
            if b["expected"] == "Increased":
                inc += 1; inc_c += ok
            else:
                dec += 1; dec_c += ok
    print(f"{label}: scored {tot} (NA {na}), correct {cor} ({100*cor/tot:.1f}%), always-increased {inc} ({100*inc/tot:.1f}%),"
          f" increases {inc_c}/{inc}, decreases {dec_c}/{dec}")
acc(v03, "Harvey v0.3")
acc(v04, "Harvetta v0.4")
# same for the matrix own calls
mat = [{"biomarkers": [{"predicted": e["predicted"], "expected": e["expected"]} for e in r["readouts"] if e["own"]]} for r in merged]
acc(mat, "matrix own")
# the recorded 'correct' field for v0.3/v0.4
for res, lab in ((v03, "v0.3"), (v04, "v0.4")):
    c = Counter(b["correct"] for r in res for b in r["biomarkers"])
    print(lab, "recorded correct field", c)
