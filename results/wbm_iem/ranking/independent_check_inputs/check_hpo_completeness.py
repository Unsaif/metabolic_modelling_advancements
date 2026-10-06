"""Compare, per protocol IEM, the full Orphanet HPO annotation set (product4) with what reaches the builder
(lab hpo_evidence terms + HPO-only rows). Uses the linker's own term heuristic for the 'metabolite' definition,
and a broader regex to flag metabolite-like terms the heuristic misses."""
import csv, json, re, sys, importlib.util, collections
import xml.etree.ElementTree as ET
SCR = "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad"
spec = importlib.util.spec_from_file_location("lk", "/home/claude/mma/scripts/link_iem_ground_truth.py")
lk = importlib.util.module_from_spec(spec); spec.loader.exec_module(lk)
tree = ET.parse(f"{SCR}/dl_orphadata/en_product4.xml")
orpha_hpo = collections.defaultdict(list)
for d in tree.getroot().iter("Disorder"):
    o = d.findtext("OrphaCode")
    for a in d.iter("HPODisorderAssociation"):
        orpha_hpo[o].append((a.findtext("HPO/HPOId"), a.findtext("HPO/HPOTerm"), a.findtext("HPOFrequency/Name") or ""))
prot = json.load(open("/home/claude/mma/data/iem/iem_protocol_v0.2.json"))
names = [p["iem"] for p in prot]
lab = list(csv.DictReader(open("/home/claude/mma/data/iem/iem_biomarkers_v0.1_linked.tsv"), delimiter="\t"))
only = list(csv.DictReader(open("/home/claude/mma/data/iem/hpo_metabolite_tuples_not_in_lab_table.tsv"), delimiter="\t"))
orpha = {}
for r in lab:
    orpha.setdefault(r["iem_abbr"], set()).add(r["orpha_code"])
print("IEMs with >1 orpha code:", {k: v for k, v in orpha.items() if len(v) > 1})
by_orpha = collections.defaultdict(list)
for i in names:
    for o in orpha[i]:
        if o: by_orpha[o].append(i)
print("protocol IEMs sharing an Orpha code:", {o: v for o, v in by_orpha.items() if len(v) > 1})
print("protocol IEMs without Orpha code:", [i for i in names if not any(orpha[i])])
# what reaches the builder: every HP id in evidence (any regex) + HPO-only rows
reach = collections.defaultdict(dict)
for r in lab:
    for hid, label, fr in re.findall(r"(HP:\d{7})\s+(.*?)\s*\[([^\]]*)\]", r["hpo_evidence"]):
        reach[r["iem_abbr"]].setdefault(hid, set()).add(fr)
    if r["hpo_evidence"] and not re.fullmatch(r"HP:\d{7} [^\[\]]+ \[[^\]]+\]", r["hpo_evidence"]):
        print("UNUSUAL evidence format:", r["iem_abbr"], r["hpo_evidence"])
for r in only:
    reach[r["iem_abbr"]].setdefault(r["hpo_id"], set()).add(r["frequency"])
BROAD = re.compile(r"(emia|aemia|uria|acidemia|aciduria|concentration|level|excretion|acidosis|ketosis|hypoglyc|hyperglyc)", re.I)
lost_total = 0
missed_broad = collections.Counter()
for i in names:
    o = next(iter(orpha[i]))
    if not o:
        continue
    anns = orpha_hpo.get(o, [])
    heur = {(h, t, f) for h, t, f in anns if lk.hpo_term_to_tuple(t)}
    lost = [(h, t, f) for h, t, f in sorted(heur) if h not in reach[i]]
    freq_mismatch = [(h, f, reach[i][h]) for h, t, f in heur if h in reach[i] and reach[i][h] != {f}]
    extra = [h for h in reach[i] if h not in {x[0] for x in heur}]
    broad = [(h, t, f) for h, t, f in anns if BROAD.search(t) and not lk.hpo_term_to_tuple(t)]
    if lost or freq_mismatch or extra:
        print(f"{i:7s} orpha {o}: LOST {lost}  FREQ-MISMATCH {freq_mismatch}  NOT-IN-ORPHANET {extra}")
    lost_total += len(lost)
    for b in broad:
        missed_broad[(b[0], b[1])] += 1
        print(f"   {i:7s} metabolite-like term not parsed by heuristic: {b}")
print("total lost annotations:", lost_total)
print("distinct heuristic-missed metabolite-like terms:", len(missed_broad))
for k, v in sorted(missed_broad.items(), key=lambda x: -x[1]): print("   ", k, v)
