"""Colleague copy of a Claude Docs Word export: drop the internal "Note for Tim" section and replace the byline.
Usage: python3 clean_for_colleague.py UNPACKED_DIR "Byline text" """
import sys
from lxml import etree
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
NS = {"w": W}
d, byline = sys.argv[1], sys.argv[2]
p = f"{d}/word/document.xml"
tree = etree.parse(p)
body = tree.getroot().find("w:body", NS)
kids = list(body)
def style(el):
    s = el.find("w:pPr/w:pStyle", NS)
    return s.get(f"{{{W}}}val") if s is not None else None
def text(el):
    return "".join(t.text or "" for t in el.iter(f"{{{W}}}t"))
# byline: the first paragraph after the title
title_i = next(i for i, el in enumerate(kids) if style(el) == "Heading1")
by = kids[title_i + 1]
assert "Tim Hulshof" in text(by), text(by)
runs = by.findall("w:r", NS)
for r in runs[1:]:
    by.remove(r)
runs[0].find("w:t", NS).text = byline
# drop the Note for Tim section
start = next(i for i, el in enumerate(kids) if style(el) == "Heading2" and text(el).strip() == "Note for Tim")
end = next(i for i in range(start + 1, len(kids)) if style(kids[i]) == "Heading2")
removed = [text(el)[:60] for el in kids[start:end]]
for el in kids[start:end]:
    body.remove(el)
tree.write(p, xml_declaration=True, encoding="UTF-8", standalone=True)
print(f"removed {len(removed)} blocks; next section: {text(kids[end])!r}")
