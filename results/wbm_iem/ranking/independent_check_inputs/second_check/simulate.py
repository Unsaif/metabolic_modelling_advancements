"""Simulate corrections: SUCLA -> ORPHA:17 (with a map row for HP:0002912), and optional map changes."""
import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_v02")
import csv
import json
from collections import defaultdict
from common import get_all, ancestors, REPO

d = get_all()
orpha, terms = d["orpha"], d["terms"]
ROOTS = {"HP:0001939", "HP:0003117", "HP:0040085"}
FREQ_OK = {"Obligate (100%)", "Very frequent (99-80%)", "Frequent (79-30%)"}
protocol = json.load(open(f"{REPO}/data/iem/iem_protocol_v0.2.json"))
links = {r["iem"]: [c for c in r["v02_orpha_codes"].split(";") if c]
         for r in csv.DictReader(open(f"{REPO}/data/iem/iem_orphanet_links_v0.2.tsv"), delimiter="\t")}
tmap = {r["hpo_id"]: dict(r) for r in csv.DictReader(open(f"{REPO}/data/iem/hpo_term_readout_map_v0.2.tsv"), delimiter="\t")}
J = json.load(open(f"{REPO}/data/iem/iem_ranking_profiles_v0.2.json"))
panel = set(J["panel_protocol"])
lab = {k: [tuple(t) for t in v] for k, v in J["profiles"]["lab"].items()}

scen = sys.argv[1:]
if "sucla17" in scen:
    links["SUCLA"] = ["17"]
if "sucla17+1933" in scen:
    links["SUCLA"] = ["17", "1933"]
if "sucla17" in scen or "sucla17+1933" in scen:
    tmap["HP:0002912"] = {"decision": "include", "readout": "DM_HC00900[bc]", "direction": "Increased"}
    for h in ("HP:0012087", "HP:0011923", "HP:0011924", "HP:0008347", "HP:0001397"):
        tmap[h] = {"decision": "exclude"}
if "drop2oaa" in scen:
    links["2OAA"] = []
if "star168558" in scen:
    links["STAR"] = ["90790", "168558"]
    for h in ("HP:0012598", "HP:0003107", "HP:0001941", "HP:0030349", "HP:0008232", "HP:0011969"):
        tmap[h] = {"decision": "exclude"}
if "nobili" in scen:
    tmap["HP:0003265"] = {"decision": "exclude"}
if "gsh" in scen:
    tmap["HP:0034738"] = {"decision": "include", "readout": "DM_gthrd[bc]", "direction": "Decreased"}

tuples = defaultdict(list)
for p in protocol:
    seen = set()
    for c in links[p["iem"]]:
        for hid, term, fr in orpha[c]["annotations"]:
            if not ((ancestors(hid, terms) | {hid}) & ROOTS):
                continue
            if (hid, fr) in seen:
                continue
            seen.add((hid, fr))
            m = tmap.get(hid)
            if m is None:
                print("UNMAPPED", p["iem"], hid, term)
                continue
            if m["decision"] != "include" or fr == "Excluded (0%)":
                continue
            tuples[p["iem"]].append((m["readout"], m["direction"], fr))

def build(freq):
    prof = {}
    for iem, tl in tuples.items():
        dd = defaultdict(set)
        for ro, di, fr in tl:
            if freq and fr not in FREQ_OK:
                continue
            dd[ro].add(di)
        kept = sorted((ro, next(iter(s))) for ro, s in dd.items() if len(s) == 1)
        if kept:
            prof[iem] = kept
    return prof
hpo, hpof = build(False), build(True)
corr = {k: [t for t in v if t in set(hpo.get(k, []))] for k, v in lab.items()}
corr = {k: v for k, v in corr.items() if v}
sz = lambda pr: (len(pr), sum(len(v) for v in pr.values()))
extra = sorted({ro for pr in (hpo, hpof) for v in pr.values() for ro, _ in v} - panel)
print("scenario", scen, "hpo", sz(hpo), "hpo_frequent", sz(hpof), "corr", sz(corr), "extra", len(extra),
      "new vs v0.2 extras", sorted(set(extra) - set(J["extra_readouts"])), "dropped", sorted(set(J["extra_readouts"]) - set(extra)))
for k in ("SUCLA", "PC", "OXOP"):
    print("  ", k, "hpo", hpo.get(k), "\n      freq", hpof.get(k), "\n      corr", corr.get(k))
