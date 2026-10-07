"""Task A: independent integrity check of the joined cross-disease matrix (no project code imported)."""
import hashlib
import json
import math
import statistics
import sys
from collections import Counter

ROOT = "/home/claude/mma"
RES = f"{ROOT}/results/wbm_iem"


def load(p):
    with open(p) as fh:
        return json.load(fh)


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


merged = load(f"{RES}/Harvey_1_03d_iem_cross_ranking_v1.json")
shards = {k: load(f"{RES}/Harvey_1_03d_iem_cross_ranking_v1_shard{k}of3.json") for k in (1, 2, 3)}
protocol = load(f"{ROOT}/data/iem/iem_protocol_v0.2.json")
extra = [l.strip() for l in open(f"{ROOT}/data/iem/iem_ranking_extra_readouts_v0.2.txt") if l.strip()]

print("merged sha256", sha(f"{RES}/Harvey_1_03d_iem_cross_ranking_v1.json"))
for k in (1, 2, 3):
    print(f"shard{k} sha256", sha(f"{RES}/Harvey_1_03d_iem_cross_ranking_v1_shard{k}of3.json"))
print("merged_from recorded:", [(x['file'].split('/')[-1], x['sha256'][:12]) for x in merged[0]['provenance']['merged_from']])

# panel built independently: protocol biomarker reactions in protocol order (deduplicated), then the extra list
panel = []
for p in protocol:
    for rid, _ in p["biomarkers"]:
        if rid not in panel:
            panel.append(rid)
n_prot = len(panel)
for rid in extra:
    if rid not in panel:
        panel.append(rid)
print("panel size", len(panel), "protocol readouts", n_prot, "extra", len(panel) - n_prot,
      "extra already in protocol:", [r for r in extra if r in panel[:n_prot]])
print("panel sha256 matches provenance:",
      hashlib.sha256(json.dumps(panel).encode()).hexdigest() == merged[0]["provenance"]["panel_sha256"])

# shard k holds panel[k-1::3]
problems = []
for k, recs in shards.items():
    want = panel[k - 1::3]
    for r in recs:
        got = [e["reaction"] for e in r["readouts"]]
        if got != want:
            problems.append(f"shard {k} {r['iem']}: readout list differs from panel[{k-1}::3]")
print("shard readout lists:", "OK" if not problems else problems[:5])

# interleave and compare with the joined file
iems_m = [(r["iem"], r["call_index"]) for r in merged]
for k, recs in shards.items():
    if [(r["iem"], r["call_index"]) for r in recs] != iems_m:
        problems.append(f"shard {k}: IEM order differs")
n_diff_entries = 0
for i, mr in enumerate(merged):
    by = {}
    for k, recs in shards.items():
        for e in recs[i]["readouts"]:
            if e["reaction"] in by:
                problems.append(f"duplicate {e['reaction']}")
            by[e["reaction"]] = e
    inter = [by[r] for r in panel]
    if inter != mr["readouts"]:
        n_diff_entries += sum(1 for a, b in zip(inter, mr["readouts"]) if a != b) + abs(len(inter) - len(mr["readouts"]))
    # set-up fields across shards
    for f in ("iem", "call_index", "iem_reactions", "healthy_pin", "wb_objective_disease_feasible", "status", "notes",
              "vmax_disease", "context"):
        vals = [json.dumps(shards[k][i].get(f)) for k in (1, 2, 3)] + [json.dumps(mr.get(f))]
        if len(set(vals)) != 1:
            problems.append(f"{mr['iem']} field {f} differs: {vals}")
    vm = [shards[k][i]["vmax_healthy"] for k in (1, 2, 3)] + [mr["vmax_healthy"]]
    if max(vm) - min(vm) > 1e-9 * max(1, abs(vm[0])):
        problems.append(f"{mr['iem']} vmax differs {vm}")
print("interleaved entries differing from joined file:", n_diff_entries)
print("set-up field problems:", problems[:10] if problems else "none")

# provenance identical apart from shard fields
SK = {"shard", "readouts_sha256", "n_readouts"}
provs = [json.dumps({a: b for a, b in recs[0]["provenance"].items() if a not in SK}, sort_keys=True) for recs in shards.values()]
print("shard provenance identical apart from shard fields:", len(set(provs)) == 1)
fps = Counter(r["run_fingerprint"] for recs in shards.values() for r in recs)
print("fingerprints per shard:", {k: len({r['run_fingerprint'] for r in recs}) for k, recs in shards.items()})

# completeness, statuses, absent readout
print("IEMs", len(merged), "statuses", Counter(r["status"] for r in merged))
print("readouts per IEM", Counter(len(r["readouts"]) for r in merged),
      "same order for all:", len({tuple(e['reaction'] for e in r['readouts']) for r in merged}) == 1)
st = Counter((e["status_healthy"], e["status_disease"]) for r in merged for e in r["readouts"])
print("status pairs", st)
absent = Counter(e["reaction"] for r in merged for e in r["readouts"] if "absent" in (e["status_healthy"], e["status_disease"]))
print("absent readouts", absent)
pred = Counter(e["predicted"] for r in merged for e in r["readouts"])
print("predicted calls", pred)
fb = sum(e["fallback_healthy"] or e["fallback_disease"] for r in merged for e in r["readouts"])
print("fallbacks", fb)
notes = [(r["iem"], r["notes"]) for r in merged if r["notes"]]
print("notes", notes)

# call rule from stored values
TOL = 1e-6
bad_rule, bad_zero = [], []
for r in merged:
    for e in r["readouts"]:
        h, d = e["healthy"], e["disease"]
        if e["status_healthy"] != "Optimal" or e["status_disease"] != "Optimal":
            want = "NA"
        else:
            for v in (h, d):
                if v is not None and v != 0 and abs(v) <= TOL:
                    bad_zero.append((r["iem"], e["reaction"], v))
            diff = d - h
            want = "Increased" if diff > TOL else "Decreased" if diff < -TOL else "Unchanged"
        if want != e["predicted"]:
            bad_rule.append((r["iem"], e["reaction"], h, d, e["predicted"], want))
print("calls not following the rule:", len(bad_rule), bad_rule[:5])
print("stored values with 0 < |f| <= 1e-6 (should have been zeroed):", len(bad_zero), bad_zero[:5])
# near-threshold differences
near = [(r["iem"], e["reaction"], e["disease"] - e["healthy"]) for r in merged for e in r["readouts"]
        if e["predicted"] != "NA" and 1e-6 < abs(e["disease"] - e["healthy"]) <= 1e-5]
print("calls decided by a difference in (1e-6, 1e-5]:", len(near))
neg = [(r["iem"], e["reaction"], e["healthy"], e["disease"]) for r in merged for e in r["readouts"]
       if e["predicted"] != "NA" and (e["healthy"] < 0 or e["disease"] < 0)]
print("negative maxima:", len(neg), neg[:5])

# own flags and expected directions against the protocol
own_bad = 0
for r, p in zip(merged, protocol):
    assert r["iem"] == p["iem"] and r["call_index"] == p["call_index"]
    exp = {rid: ("Increased" if lab.lower().startswith("increased") else "Decreased" if lab.lower().startswith("decreased") else lab)
           for rid, lab in p["biomarkers"]}
    for e in r["readouts"]:
        if e["own"] != (e["reaction"] in exp) or e["expected"] != exp.get(e["reaction"]):
            own_bad += 1
print("own/expected flags inconsistent with protocol:", own_bad)

# rechecks
rows = []
for r in merged:
    for e in r["readouts"]:
        for s in ("healthy", "disease"):
            info = e[f"solve_{s}"] or {}
            if "recheck" in info:
                rows.append((r["iem"], e["reaction"], s, e[s], info["recheck"], info.get("method")))
print("rechecks", len(rows), "statuses", Counter(x[4]["status"] for x in rows), "methods of checked solves", Counter(x[5] for x in rows))
print("recorded abs_diff max", max(x[4]["abs_diff"] for x in rows))
zero_rec = [x for x in rows if x[4]["value"] is not None and (0 if abs(x[4]["value"]) <= TOL else x[4]["value"]) != x[3]]
print("recheck value (zeroed by rule) != stored value:", len(zero_rec))
print("max |recheck value - stored value| (stored is zeroed):", max(abs(x[4]["value"] - x[3]) for x in rows))
ts = [x[4]["time_s"] for x in rows]
print("recheck time_s: n", len(ts), "max", max(ts), "n == 0.0", sum(t == 0 for t in ts), "mean", statistics.mean(ts))

# expected recheck selection by hash (independent)
def chosen(iem, state, rid):
    return int(hashlib.md5(f"{iem}|{state}|{rid}".encode()).hexdigest(), 16) % 100 == 0
exp_rc = set()
for r in merged:
    for e in r["readouts"]:
        for s in ("healthy", "disease"):
            info = e[f"solve_{s}"] or {}
            if info.get("method") not in (None, "ipm") and chosen(r["iem"], s, e["reaction"]):
                exp_rc.add((r["iem"], e["reaction"], s))
got_rc = {(x[0], x[1], x[2]) for x in rows}
print("rechecks expected by hash", len(exp_rc), "== recorded:", exp_rc == got_rc)

# timing
times = [e[f"time_{s}_s"] for r in merged for e in r["readouts"] for s in ("healthy", "disease")
         if e[f"status_{s}"] == "Optimal"]
print("readout solves", len(times), "median time", statistics.median(times), "mean", statistics.mean(times),
      "sum h", sum(times) / 3600)
ipm_t = [e["time_disease_s"] for r in merged for e in r["readouts"] if (e["solve_disease"] or {}).get("method") == "ipm"]
print("first-of-readout ipm solves", len(ipm_t), "median", statistics.median(ipm_t))
warm_t = [e[f"time_{s}_s"] for r in merged for e in r["readouts"] for s in ("healthy", "disease")
          if (e[f"solve_{s}"] or {}).get("method") == "concurrent"]
print("concurrent solves", len(warm_t), "median", statistics.median(warm_t))
setup = [(r["iem"], r["n_setup_solves"]) for r in merged]
print("setup solves per IEM", Counter(n for _, n in setup), "total", sum(n for _, n in setup), "x3 shards", 3 * sum(n for _, n in setup))
