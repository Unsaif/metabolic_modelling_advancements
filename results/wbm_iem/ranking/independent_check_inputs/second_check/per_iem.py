import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_v02")
import csv
import json
from common import get_all, ancestors, REPO

d = get_all()
orpha, terms = d["orpha"], d["terms"]
ROOTS = {"HP:0001939", "HP:0003117", "HP:0040085"}
tmap = {r["hpo_id"]: r for r in csv.DictReader(open(f"{REPO}/data/iem/hpo_term_readout_map_v0.2.tsv"), delimiter="\t")}
codes = sys.argv[1:]
for c in codes:
    v = orpha[c]
    print(f"== {c} {v['name']} [{v['type']}] n={len(v['annotations'])}")
    for hid, term, fr in sorted(v["annotations"], key=lambda a: a[1]):
        und = bool((ancestors(hid, terms) | {hid}) & ROOTS)
        m = tmap.get(hid)
        dec = "" if not und else (f"-> {m['readout']} {m['direction']}" if m and m["decision"] == "include" else ("EXCL: " + m["reason"][:50] if m else "NOT IN MAP"))
        flag = "R" if und else " "
        print(f"   {flag} {hid} {term} | {fr} {dec}")
