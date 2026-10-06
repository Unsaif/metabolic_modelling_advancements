import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_v02")
import re
from common import get_all

d = get_all()
orpha = d["orpha"]
pats = sys.argv[1:]
for p in pats:
    rx = re.compile(p, re.I)
    print(f"== {p}")
    for c, v in sorted(orpha.items(), key=lambda kv: int(kv[0])):
        if rx.search(v["name"] or ""):
            print(f"   {c}\t{v['name']}\t[{v['type']}]\tn={len(v['annotations'])}")
