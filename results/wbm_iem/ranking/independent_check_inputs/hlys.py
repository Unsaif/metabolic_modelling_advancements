import re, xml.etree.ElementTree as ET
tree = ET.parse("/tmp/claude-0/-home-claude/5c6d4c6c-0648-59bb-a3ee-169f6ec3d723/scratchpad/dl_orphadata/en_product4.xml")
for d in tree.getroot().iter("Disorder"):
    n = d.findtext("Name") or ""
    if re.search(r"lysinemia|saccharopin", n, re.I):
        print(d.findtext("OrphaCode"), n, [(a.findtext("HPO/HPOTerm"), a.findtext("HPOFrequency/Name")) for a in d.iter("HPODisorderAssociation") if re.search(r"emia|uria|level|concentration", a.findtext("HPO/HPOTerm"), re.I)])
