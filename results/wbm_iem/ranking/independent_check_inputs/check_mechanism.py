import csv, json, re, importlib.util, collections
import xml.etree.ElementTree as ET
SCR = "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad"
spec = importlib.util.spec_from_file_location("lk", "/home/claude/mma/scripts/link_iem_ground_truth.py")
lk = importlib.util.module_from_spec(spec); spec.loader.exec_module(lk)
lab = list(csv.DictReader(open("/home/claude/mma/data/iem/iem_biomarkers_v0.1_linked.tsv"), delimiter="\t"))
vname = {r["vmh_metabolite"]: r["vmh_name"] for r in lab}
for vm, hp in [("thym", "Elevated urinary thymine level"), ("thym", "Elevated urinary dihydrothymine level"),
               ("56dura", "Elevated urinary dihydrouracil level"), ("56dura", "Uraciluria"),
               ("c5dc", "Glutaric aciduria"), ("4hpro_LT", "Prolinuria"), ("4hpro_LT", "Hydroxyprolinuria"), ("3ivcrn", "3-hydroxyisovaleric aciduria"),
               ("hcys_L", "Hyperhomocystinemia"), ("Lhcystin", "Hyperhomocystinemia")]:
    t = lk.hpo_term_to_tuple(hp)
    print(f"{vm:9s} ({vname.get(vm)!r:40s}) vs {hp!r:45s} parsed={t} match={lk.met_match(vname.get(vm) or vm, t[0]) if t else None}")
tree = ET.parse(f"{SCR}/dl_orphadata/en_product4.xml")
ann = collections.defaultdict(list); oname = {}
for d in tree.getroot().iter("Disorder"):
    o = d.findtext("OrphaCode"); oname[o] = d.findtext("Name")
    for a in d.iter("HPODisorderAssociation"):
        ann[o].append((a.findtext("HPO/HPOId"), a.findtext("HPO/HPOTerm"), a.findtext("HPOFrequency/Name")))
link = {}
for r in lab:
    link.setdefault(r["iem_abbr"], (r["orpha_code"], r["orpha_name"], r["gene_symbols"], r["orpha_match_score"]))
for i in ["GA1", "HYPRO1", "IVA", "ADSL", "AKGD", "CD", "DGK", "HMET", "HCYS"]:
    o = link[i][0]
    print(f"\n== {i} -> Orpha {o} {link[i][1]!r} genes={link[i][2]} score={link[i][3]}")
    for h, t, f in ann.get(o, []):
        if lk.hpo_term_to_tuple(t) or re.search(r"(emia|uria|acid|level|concentration|glyc|lact|keto)", t, re.I):
            print("    ", h, t, "|", f)
# Orphanet disorders by name for ADSL / Canavan
for o, n in oname.items():
    if re.search(r"adenylosuccin|canavan", n or "", re.I):
        print("\n## candidate", o, n)
        for h, t, f in ann.get(o, []):
            if lk.hpo_term_to_tuple(t) or re.search(r"(emia|uria|acid|level|concentration)", t, re.I):
                print("    ", h, t, "|", f, "| parsed:", lk.hpo_term_to_tuple(t))
