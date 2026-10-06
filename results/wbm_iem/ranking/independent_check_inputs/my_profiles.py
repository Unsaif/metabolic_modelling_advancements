"""Independent re-derivation of the ranking profiles from the source tables (not reusing the builder's code)."""
import csv, json, re, collections
D = "/home/claude/mma/data/iem/"
prot = json.load(open(D + "iem_protocol_v0.2.json"))
lab = list(csv.DictReader(open(D + "iem_biomarkers_v0.1_linked.tsv"), delimiter="\t"))
only = list(csv.DictReader(open(D + "hpo_metabolite_tuples_not_in_lab_table.tsv"), delimiter="\t"))
hmap = {r["hpo_id"]: r for r in csv.DictReader(open(D + "hpo_term_readout_map_v0.1.tsv"), delimiter="\t")}
names = [p["iem"] for p in prot]
call_index = {p["iem"]: str(p["call_index"]) for p in prot}

def lab_dir(label):
    w = label.strip().split()[0]
    assert w in ("Increased", "Decreased"), label
    return w

def dedupe(tuples):
    by = collections.defaultdict(set)
    for r, s in tuples: by[r].add(s)
    conflicts = sorted(r for r, s in by.items() if len(s) > 1)
    return sorted((r, next(iter(s))) for r, s in by.items() if len(s) == 1), conflicts

mine = {k: {} for k in ("lab", "lab_hpo_corroborated", "hpo", "hpo_frequent")}
conf = collections.defaultdict(dict)
# lab and corroborated: match lab-table row by (iem, call_index, reaction, direction)
labrow = {}
for r in lab:
    labrow.setdefault((r["iem_abbr"], r["call_index"], r["biomarker_reaction"], r["expected_direction"]), []).append(r)
for p in prot:
    tl = [(rx, lab_dir(l)) for rx, l in p["biomarkers"]]
    mine["lab"][p["iem"]], conf["lab"][p["iem"]] = dedupe(tl)
    corr = []
    for rx, s in tl:
        rows = labrow.get((p["iem"], str(p["call_index"]), rx, s), [])
        assert len(rows) == 1, (p["iem"], rx, len(rows))
        if rows[0]["hpo_check"] == "corroborated":
            corr.append((rx, s))
    t, c = dedupe(corr)
    if t: mine["lab_hpo_corroborated"][p["iem"]] = t
# HPO annotations: every HP id appearing anywhere in hpo_evidence (any separator), plus HPO-only rows
ann = collections.defaultdict(dict)   # iem -> hid -> set(freq)
for r in lab:
    if r["iem_abbr"] not in names: continue
    for hid, lbl, fr in re.findall(r"(HP:\d+)\s+([^\[\];|]*?)\s*\[([^\]]*)\]", r["hpo_evidence"]):
        ann[r["iem_abbr"]].setdefault(hid, set()).add(fr)
    n_ids = len(re.findall(r"HP:\d+", r["hpo_evidence"]))
    assert n_ids == len(re.findall(r"(HP:\d+)\s+([^\[\];|]*?)\s*\[([^\]]*)\]", r["hpo_evidence"])), r["hpo_evidence"]
for r in only:
    if r["iem_abbr"] in names:
        ann[r["iem_abbr"]].setdefault(r["hpo_id"], set()).add(r["frequency"])
assert all(len(f) == 1 for a in ann.values() for f in a.values()), "frequency conflict"
FREQ = {"Obligate (100%)", "Very frequent (99-80%)", "Frequent (79-30%)"}
terms_used = set()
for iem in names:
    tup_all, tup_freq = [], []
    for hid, fr in ann.get(iem, {}).items():
        fr = next(iter(fr)); terms_used.add(hid)
        m = hmap[hid]
        if m["decision"] != "include" or fr == "Excluded (0%)": continue
        tup_all.append((m["readout"], m["direction"]))
        if fr in FREQ: tup_freq.append((m["readout"], m["direction"]))
    for key, tl in (("hpo", tup_all), ("hpo_frequent", tup_freq)):
        t, c = dedupe(tl)
        if t: mine[key][iem] = t
        if c: conf[key][iem] = c
print("distinct HPO terms annotated to protocol IEMs:", len(terms_used), "; map rows:", len(hmap),
      "; map rows unused by protocol IEMs:", sorted(set(hmap) - terms_used))
print("frequencies seen:", collections.Counter(next(iter(f)) for a in ann.values() for f in a.values()))
committed = json.load(open(D + "iem_ranking_profiles_v0.1.json"))
for k in mine:
    cp = {i: sorted(map(tuple, v)) for i, v in committed["profiles"][k].items() if v}
    mp = mine[k]
    print(f"{k:22s} mine {len(mp)} profiles / {sum(map(len, mp.values()))} tuples; committed {len(cp)} / {sum(map(len, cp.values()))}; "
          f"identical={cp == mp}; conflicts={dict((i, c) for i, c in conf[k].items() if c)}")
    for i in sorted(set(cp) | set(mp)):
        if cp.get(i) != mp.get(i): print("   DIFF", i, cp.get(i), mp.get(i))
# extra readouts
panel = list(dict.fromkeys(r for p in prot for r, _ in p["biomarkers"]))
extra = sorted({r for k in ("hpo", "hpo_frequent") for v in mine[k].values() for r, _ in v} - set(panel))
print("extra readouts mine:", len(extra), extra == committed["extra_readouts"],
      extra == [l.strip() for l in open(D + "iem_ranking_extra_readouts_v0.1.txt") if l.strip()])
print("panel", len(panel), panel == committed["panel_protocol"])
print("candidates match protocol order:", [c["iem"] for c in committed["candidates"]] == names)
# frequency-excluded annotations
print("Excluded (0%) annotations among protocol IEMs:", [(i, h) for i, a in ann.items() for h, f in a.items() if "Excluded (0%)" in f])
json.dump({k: v for k, v in mine.items()}, open("/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_ranking/my_profiles.json", "w"), indent=1)
