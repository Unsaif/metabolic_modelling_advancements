import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_v02")
import csv
import json
from collections import defaultdict, Counter
from common import get_all, ancestors, REPO

d = get_all()
orpha, terms, H = d["orpha"], d["terms"], d["harvey"]
ROOTS = {"HP:0001939", "HP:0003117", "HP:0040085"}
FREQ_OK = {"Obligate (100%)", "Very frequent (99-80%)", "Frequent (79-30%)"}

protocol = json.load(open(f"{REPO}/data/iem/iem_protocol_v0.2.json"))
links = {r["iem"]: r for r in csv.DictReader(open(f"{REPO}/data/iem/iem_orphanet_links_v0.2.tsv"), delimiter="\t")}
tmap = {r["hpo_id"]: r for r in csv.DictReader(open(f"{REPO}/data/iem/hpo_term_readout_map_v0.2.tsv"), delimiter="\t")}
J = json.load(open(f"{REPO}/data/iem/iem_ranking_profiles_v0.2.json"))

# --- lab profiles from the protocol, independently
def dirn(label):
    l = label.lower()
    if "increase" in l:
        return "Increased"
    if "decrease" in l:
        return "Decreased"
    return None

lab = {}
for p in protocol:
    tl = []
    for r, lbl in p["biomarkers"]:
        dd = dirn(lbl)
        assert dd, (p["iem"], r, lbl)
        tl.append((r, dd))
    lab[p["iem"]] = sorted(set(tl))
print("lab", len(lab), sum(len(v) for v in lab.values()), Counter(dd for v in lab.values() for _, dd in v))

# --- annotations
anc_cache = {}
def under(h):
    if h not in anc_cache:
        anc_cache[h] = (h in ROOTS) or bool(ancestors(h, terms) & ROOTS)
    return anc_cache[h]

ann_considered = []
missing_terms, obsolete = [], []
for p in protocol:
    codes = [c for c in links[p["iem"]]["v02_orpha_codes"].split(";") if c]
    for c in codes:
        for hid, term, fr in orpha[c]["annotations"]:
            if hid not in terms:
                missing_terms.append((p["iem"], c, hid, term))
                continue
            if terms[hid]["obsolete"]:
                obsolete.append((p["iem"], c, hid, term))
            if under(hid):
                ann_considered.append((p["iem"], c, hid, term, fr))
print("missing hpo ids", missing_terms)
print("obsolete", obsolete)
print("annotations considered", len(ann_considered), "distinct terms", len({a[2] for a in ann_considered}))
dups = Counter((a[0], a[2]) for a in ann_considered)
print("dup (iem,term)", [k for k, v in dups.items() if v > 1])

# --- map coverage
considered_terms = {a[2] for a in ann_considered}
print("terms in map", len(tmap), "included", sum(1 for r in tmap.values() if r["decision"] == "include"))
print("considered terms not in map:", sorted(considered_terms - set(tmap)))
print("map terms not considered:", sorted(set(tmap) - considered_terms))
# label in map equals obo name?
for h, r in tmap.items():
    if h in terms and terms[h]["name"] != r["hpo_term"]:
        print("map label differs from obo:", h, r["hpo_term"], "|", terms[h]["name"])

tuples = defaultdict(list)
used = 0
per_ann = []
for iem, c, hid, term, fr in ann_considered:
    r = tmap[hid]
    if r["decision"] != "include":
        continue
    if fr == "Excluded (0%)":
        continue
    used += 1
    tuples[iem].append((r["readout"], r["direction"], fr, hid))
print("used annotations", used)

def build(freq_filter):
    prof, confl = {}, {}
    for iem, tl in tuples.items():
        dd = defaultdict(set)
        for ro, di, fr, hid in tl:
            if freq_filter and fr not in FREQ_OK:
                continue
            dd[ro].add(di)
        kept = sorted((ro, next(iter(s))) for ro, s in dd.items() if len(s) == 1)
        dropped = sorted(ro for ro, s in dd.items() if len(s) > 1)
        if kept:
            prof[iem] = kept
        if dropped:
            confl[iem] = dropped
    return prof, confl

hpo, c1 = build(False)
hpof, c2 = build(True)
corr = {}
for iem, tl in lab.items():
    s = set(hpo.get(iem, []))
    cc = [t for t in tl if t in s]
    if cc:
        corr[iem] = cc
def size(pr):
    return len(pr), sum(len(v) for v in pr.values())
print("hpo", size(hpo), "conflicts", c1)
print("hpo_frequent", size(hpof), "conflicts", c2)
print("corroborated", size(corr))

# compare with JSON
def norm(pr):
    return {k: sorted(tuple(t) for t in v) for k, v in pr.items() if v}
for name, mine in (("hpo", hpo), ("hpo_frequent", hpof), ("lab_hpo_corroborated", corr), ("lab", lab)):
    theirs = norm(J["profiles"][name])
    mine_n = norm(mine)
    if theirs == mine_n:
        print(f"MATCH {name}")
    else:
        for k in sorted(set(theirs) | set(mine_n)):
            if theirs.get(k) != mine_n.get(k):
                print(f"DIFF {name} {k}: theirs={theirs.get(k)} mine={mine_n.get(k)}")

# extra readouts
panel = set(J["panel_protocol"])
needed = {ro for pr in (hpo, hpof) for v in pr.values() for ro, _ in v}
extra = sorted(needed - panel)
print("extra", len(extra), extra == J["extra_readouts"])
ext_file = [l.strip() for l in open(f"{REPO}/data/iem/iem_ranking_extra_readouts_v0.2.txt") if l.strip()]
print("extra file equal", ext_file == extra)
# readouts used by profiles but exist?
rx = set(H["rxns"]); mets = set(H["mets"])
def exists(r):
    return r in rx or (r.startswith("DM_") and r[3:] in mets)
print("needed readouts missing in Harvey:", [r for r in sorted(needed) if not exists(r)])
print("panel readouts missing in Harvey:", [r for r in sorted(panel) if not exists(r)])
# panel check: protocol biomarkers
panel_mine = []
for p in protocol:
    for r, _ in p["biomarkers"]:
        if r not in panel_mine:
            panel_mine.append(r)
print("panel", len(panel_mine), set(panel_mine) == panel)
# annotations in JSON
ja = J["hpo_annotations"]
print("JSON annotations", len(ja), "used", sum(a["used"] for a in ja))
mine_set = Counter((a[0], a[1], a[2], a[4]) for a in ann_considered)
their_set = Counter((a["iem"], a["orpha_code"], a["hpo_id"], a["frequency"]) for a in ja)
print("annotation sets equal", mine_set == their_set)
json.dump({"hpo": hpo, "hpo_frequent": hpof, "corr": corr, "tuples": {k: v for k, v in tuples.items()},
           "considered": ann_considered}, open("/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_v02/mine.json", "w"), indent=1)
