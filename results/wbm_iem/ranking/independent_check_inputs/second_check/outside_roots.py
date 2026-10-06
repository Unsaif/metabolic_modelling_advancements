import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_v02")
import csv
import json
import re
from collections import defaultdict
from common import get_all, ancestors, REPO

d = get_all()
orpha, terms = d["orpha"], d["terms"]
ROOTS = {"HP:0001939", "HP:0003117", "HP:0040085"}
links = {r["iem"]: r for r in csv.DictReader(open(f"{REPO}/data/iem/iem_orphanet_links_v0.2.tsv"), delimiter="\t")}
out = defaultdict(list)
for iem, r in links.items():
    for c in [x for x in r["v02_orpha_codes"].split(";") if x]:
        for hid, term, fr in orpha[c]["annotations"]:
            anc = ancestors(hid, terms) | {hid}
            if anc & ROOTS:
                continue
            out[hid].append((iem, fr))
# print all distinct non-root terms with top-level branch
TOP = {t for t, v in terms.items() if "HP:0000118" in v["is_a"]}
def tops(h):
    return sorted(terms[t]["name"] for t in (ancestors(h, terms) | {h}) & TOP)
kw = re.compile(r"(uria|emia|aemia|concentration|level|circulating|serum|plasma|CSF|cerebrospinal|urine|urinary|blood|acid|ketosis|acidosis|alkalosis|crystal|stone|lithiasis|deposit|accumul|excretion)", re.I)
allterms = sorted(out, key=lambda h: terms[h]["name"] or "")
print("non-root distinct terms:", len(allterms))
for h in allterms:
    nm = terms[h]["name"]
    if kw.search(nm or "") or (len(sys.argv) > 1 and re.search(sys.argv[1], nm or "", re.I)):
        print(f"{h}\t{nm}\t{tops(h)}\t{out[h]}")
