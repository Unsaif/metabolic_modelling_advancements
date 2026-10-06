import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_v02")
import json
from collections import defaultdict
from common import get_all, REPO

d = get_all()
H = d["harvey"]
rxns, rn, gr = H["rxns"], H["rxnNames"], H.get("grRules")
protocol = {p["iem"]: p for p in json.load(open(f"{REPO}/data/iem/iem_protocol_v0.2.json"))}
for iem in sys.argv[1:]:
    p = protocol[iem]
    print(f"== {iem} include={p['include_patterns']} exclude={p['exclude_patterns']}")
    for pat in p["include_patterns"]:
        hits = [i for i, r in enumerate(rxns) if pat in r and not any(x in r for x in p["exclude_patterns"])]
        names = defaultdict(int)
        rules = defaultdict(int)
        for i in hits:
            names[rn[i]] += 1
            rules[gr[i] if isinstance(gr, list) else ""] += 1
        print(f"   {pat}: {len(hits)} rxns; names={dict(names)}")
        print(f"        grRules={dict(list(rules.items())[:3])}")
