import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_v02")
import csv
import re
import sys
from common import get_all, REPO

d = get_all()
orpha = d["orpha"]
links = list(csv.DictReader(open(f"{REPO}/data/iem/iem_orphanet_links_v0.2.tsv"), delimiter="\t"))
from collections import Counter
print("n rows", len(links), Counter(r["change"] for r in links))
for r in links:
    v1 = [c for c in r["v01_orpha_codes"].split(";") if c]
    v2 = [c for c in r["v02_orpha_codes"].split(";") if c]
    s1 = "; ".join(f"{c}={orpha[c]['name']!r}[{orpha[c]['type']}] n={len(orpha[c]['annotations'])}" if c in orpha else f"{c}=NOT_IN_P4" for c in v1)
    s2 = "; ".join(f"{c}={orpha[c]['name']!r}[{orpha[c]['type']}] n={len(orpha[c]['annotations'])}" if c in orpha else f"{c}=NOT_IN_P4" for c in v2)
    # check v02 name and count columns
    chk = ""
    if v2:
        nm = "; ".join(orpha[c]["name"] for c in v2 if c in orpha)
        na = sum(len(orpha[c]["annotations"]) for c in v2 if c in orpha)
        if nm != r["v02_orpha_names"]:
            chk += f" NAME_MISMATCH({nm!r} vs {r['v02_orpha_names']!r})"
        if str(na) != r["v02_n_hpo_annotations"]:
            chk += f" COUNT_MISMATCH({na} vs {r['v02_n_hpo_annotations']})"
    print(f"{r['iem']:7s} {r['change']:10s} v01: {s1} || v02: {s2}{chk}")
