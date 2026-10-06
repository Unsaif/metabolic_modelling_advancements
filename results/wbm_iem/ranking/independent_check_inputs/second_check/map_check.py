import sys; sys.path.insert(0, "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/independent_check_v02")
import csv
import json
import re
from collections import defaultdict
from common import get_all, ancestors, REPO

d = get_all()
orpha, terms, H = d["orpha"], d["terms"], d["harvey"]
tmap = list(csv.DictReader(open(f"{REPO}/data/iem/hpo_term_readout_map_v0.2.tsv"), delimiter="\t"))
mets = H["mets"]; names = H["metNames"]; rxns = set(H["rxns"])
hmdb = H.get("metHMDBID") if isinstance(H.get("metHMDBID"), list) else None
base_info = defaultdict(lambda: {"names": set(), "comps": set(), "hmdb": set()})
for i, m in enumerate(mets):
    mm = re.match(r"^(.*)\[(\w+)\]$", m)
    if not mm:
        continue
    full, comp = mm.groups()
    # organ prefix?
    base = full
    org = None
    if "_" in full:
        pre, rest = full.split("_", 1)
        if pre in {"Liver", "Brain", "Kidney", "Heart", "Muscle", "Lung", "Adipocytes", "Spleen", "Pancreas", "Stomach",
                   "Colon", "sIEC", "Gall", "RBC", "Skin", "Thyroidgland", "Agland", "Retina", "Urinarybladder", "Prostate",
                   "Testis", "Platelet", "Monocyte", "Nkcells", "CD4Tcells", "Bcells", "BBB", "Breast", "Cervix", "Ovary",
                   "Uterus", "Scord", "Parathyroidglands", "Adrenal"}:
            org, base = pre, rest
    base_info[base]["names"].add(names[i])
    base_info[base]["comps"].add((org or "") + ":" + comp)
    if hmdb:
        base_info[base]["hmdb"].add(hmdb[i])
mode = sys.argv[1] if len(sys.argv) > 1 else "include"
for r in tmap:
    if r["decision"] != mode:
        continue
    t = terms.get(r["hpo_id"], {})
    out = f"{r['hpo_id']} map='{r['hpo_term']}' obo='{t.get('name')}'\n   def={t.get('def')}\n"
    if mode == "include":
        met = r["vmh_metabolite"]
        bi = base_info.get(met)
        ro = r["readout"]
        ok = ro in rxns or (ro.startswith("DM_") and ro[3:] in set(mets))
        exp = {"blood": f"DM_{met}[bc]", "urine": f"EX_{met}[u]", "csf": f"DM_{met}[csf]"}[r["biofluid"]]
        out += f"   -> {met} {r['biofluid']} {r['direction']} {ro} exists={ok} convention={'OK' if exp == ro else 'BAD:' + exp}\n"
        if bi:
            comps = sorted(bi['comps'])
            out += f"   harvey names={sorted(bi['names'])} hmdb={sorted(bi['hmdb'])[:3]} ncomps={len(comps)} bc/u/csf={[c for c in comps if c.split(':')[1] in ('bc','u','csf') and c.split(':')[0]=='']}\n"
        else:
            out += "   harvey: NOT FOUND\n"
        out += f"   reason: {r['reason']}\n"
    else:
        out += f"   reason: {r['reason']}\n"
    print(out)
