import json, re, csv, collections
H = json.load(open("/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_ranking/harvey_ids.json"))
fluid_bases = collections.defaultdict(set); hm = {}; nm = {}
for m, n, h in zip(H["mets"], H["metNames"], H["metHMDBID"]):
    mm = re.match(r"^([^\[]+)\[(bc|u|csf)\]$", m)
    if mm:
        b = mm.group(1); fluid_bases[b].add(mm.group(2)); hm[b] = h; nm[b] = n
by_hmdb = collections.defaultdict(set); by_name = collections.defaultdict(set)
for b in fluid_bases:
    if hm[b]: by_hmdb[hm[b]].add(b)
    by_name[re.sub(r"[^a-z0-9]", "", nm[b].lower())].add(b)
rows = [r for r in csv.DictReader(open("/home/claude/mma/data/iem/hpo_term_readout_map_v0.1.tsv"), delimiter="\t") if r["decision"] == "include"]
for r in rows:
    b = r["vmh_metabolite"]
    alts = (by_hmdb.get(hm.get(b), set()) | by_name.get(re.sub(r"[^a-z0-9]", "", nm.get(b, "").lower()), set())) - {b}
    if alts: print(r["hpo_id"], r["hpo_term"], b, "-> other fluid ids for same compound:", {a: (nm[a], sorted(fluid_bases[a])) for a in alts})
