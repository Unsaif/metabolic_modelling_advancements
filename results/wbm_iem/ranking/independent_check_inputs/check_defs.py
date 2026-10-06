import csv, re
SCR = "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad"
terms = {}; cur = None
for line in open(f"{SCR}/dl_hpo/hp.obo"):
    line = line.rstrip("\n")
    if line == "[Term]":
        cur = {"is_a": [], "syn": []}; continue
    if line.startswith("[") and line.endswith("]"):
        cur = None; continue
    if cur is None: continue
    k, _, v = line.partition(": ")
    if k == "id": terms[v] = cur; cur["id"] = v
    elif k == "name": cur["name"] = v
    elif k == "def": cur["def"] = v.split('" [')[0].strip('"')
    elif k == "is_a": cur["is_a"].append(v.split(" ! ")[0])
    elif k == "is_obsolete": cur["obsolete"] = True
    elif k == "replaced_by": cur["replaced_by"] = v
    elif k == "synonym": cur["syn"].append(v.split('"')[1])
rows = list(csv.DictReader(open("/home/claude/mma/data/iem/hpo_term_readout_map_v0.1.tsv"), delimiter="\t"))
for r in rows:
    t = terms.get(r["hpo_id"], {})
    flag = "" if t.get("name") == r["hpo_term"] else f"  LABEL DIFFERS: obo={t.get('name')!r}"
    if t.get("obsolete"): flag += "  OBSOLETE"
    print(f"{r['hpo_id']} [{r['decision']:7s}] {r['hpo_term']}{flag}\n      def: {t.get('def','(none)')}\n      -> {r['readout'] or r['reason']}")
print("\n### lost / heuristic-missed terms")
for h in ["HP:6000118", "HP:6000331", "HP:0001958", "HP:0001988", "HP:0003162", "HP:0003128", "HP:0002919", "HP:0008281", "HP:0002908", "HP:0005979", "HP:0001942"]:
    t = terms[h]; print(h, t["name"], "| is_a", [terms[p]["name"] for p in t["is_a"]], "\n      def:", t.get("def"))
