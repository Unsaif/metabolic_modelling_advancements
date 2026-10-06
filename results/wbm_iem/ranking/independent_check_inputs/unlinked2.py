import re, json, importlib.util, xml.etree.ElementTree as ET
SCR = "/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad"
spec = importlib.util.spec_from_file_location("lk", "/home/claude/mma/scripts/link_iem_ground_truth.py")
lk = importlib.util.module_from_spec(spec); spec.loader.exec_module(lk)
tree = ET.parse(f"{SCR}/dl_orphadata/en_product4.xml")
for d in tree.getroot().iter("Disorder"):
    n = d.findtext("Name") or ""
    if re.search(r"pentosuria|branched-chain amino acid transaminase|BCAT|leucine-isoleucin|LTC4|leukotriene", n, re.I):
        terms = [(a.findtext("HPO/HPOId"), a.findtext("HPO/HPOTerm"), a.findtext("HPOFrequency/Name")) for a in d.iter("HPODisorderAssociation")]
        print(d.findtext("OrphaCode"), n, "| parsed:", [t for t in terms if lk.hpo_term_to_tuple(t[1])], "| n terms", len(terms))
H = json.load(open(f"{SCR}/independent_check_ranking/harvey_ids.json"))
for m, nm in zip(H["mets"], H["metNames"]):
    if re.search(r"adip", nm, re.I) and re.match(r"^[^_]+\[(bc|u)\]$", m): print("  Harvey:", m, nm)
prot = json.load(open("/home/claude/mma/data/iem/iem_protocol_v0.2.json"))
for p in prot:
    if p["iem"] in ("2OAA", "STAR", "EP", "HYPVLI", "LTC4S"): print(p["iem"], p["biomarkers"])
