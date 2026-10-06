"""Corroboration by readout identity: lab tuple (r, s) is corroborated iff an HPO annotation of the IEM's linked
Orphanet disorder (product4, not 'Excluded (0%)') maps (via the HPO map, plus the two lost DPYR terms) to (r, s)."""
import csv, json, re, collections, importlib.util
import xml.etree.ElementTree as ET
SCR = "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad"
D = "/home/claude/mma/data/iem/"
prot = json.load(open(D + "iem_protocol_v0.2.json"))
lab = list(csv.DictReader(open(D + "iem_biomarkers_v0.1_linked.tsv"), delimiter="\t"))
hmap = {r["hpo_id"]: r for r in csv.DictReader(open(D + "hpo_term_readout_map_v0.1.tsv"), delimiter="\t")}
extra_map = {"HP:6000118": ("EX_56dura[u]", "Increased"), "HP:6000331": ("EX_thym[u]", "Increased")}
tree = ET.parse(f"{SCR}/dl_orphadata/en_product4.xml")
ann = collections.defaultdict(list)
for d in tree.getroot().iter("Disorder"):
    o = d.findtext("OrphaCode")
    for a in d.iter("HPODisorderAssociation"):
        ann[o].append((a.findtext("HPO/HPOId"), a.findtext("HPO/HPOTerm"), a.findtext("HPOFrequency/Name")))
orpha = {}
for r in lab:
    orpha.setdefault((r["iem_abbr"], r["call_index"]), r["orpha_code"])
committed = json.load(open(D + "iem_ranking_profiles_v0.1.json"))["profiles"]["lab_hpo_corroborated"]
committed = {i: set(map(tuple, v)) for i, v in committed.items() if v}
mine = {}
for p in prot:
    o = orpha[(p["iem"], str(p["call_index"]))]
    mapped = {}
    for h, t, f in ann.get(o, []):
        if f == "Excluded (0%)": continue
        if h in hmap and hmap[h]["decision"] == "include":
            mapped[(hmap[h]["readout"], hmap[h]["direction"])] = f"{h} {t}"
        elif h in extra_map:
            mapped[extra_map[h]] = f"{h} {t}"
    s = set()
    for rx, l in p["biomarkers"]:
        key = (rx, l.split()[0])
        if key in mapped: s.add(key)
    if s: mine[p["iem"]] = s
    c = committed.get(p["iem"], set())
    if s != c:
        print(f"{p['iem']:7s} only in committed: {sorted(c - s)}   only by readout identity: {sorted(s - c)}"
              + "".join(f"\n          {k} <- {mapped[k]}" for k in sorted(s - c)))
print("readout-identity corroborated:", len(mine), "profiles,", sum(map(len, mine.values())), "tuples; committed:",
      len(committed), "/", sum(map(len, committed.values())))
