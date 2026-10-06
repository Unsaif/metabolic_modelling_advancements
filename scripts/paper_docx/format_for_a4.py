"""A4 page size and readable table column widths for the colleague copies (formatting only, no text changes,
except dropping Paper 2's duplicated figure title line, which repeats the caption beneath it).
Usage: python3 format_for_a4.py UNPACKED_DIR paper1|paper2"""
import sys
from lxml import etree
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"; NS = {"w": W}
q = lambda t: f"{{{W}}}{t}"
d, which = sys.argv[1], sys.argv[2]
TEXT_W = 11906 - 2 * 1440
WIDTHS = {
    "paper1": [[1290, 1840, 1085, 1185, 850, 755, 950, 1071],
               [2200, 1150, 850, 1250, 1100, 1150, 1326],
               [1900, 2300, 1500, 1300, 2026],
               [2526, 1300, 1300, 1300, 1300, 1300],
               [2200, 4626, 2200],
               [3200, 5826],
               [900, 4126, 2200, 1800]],
    "paper2": [[4026, 2500, 2500],
               [2500, 6526],
               [2500, 6526],
               [1526, 1500, 1500, 1500, 1500, 1500],
               [1100, 2000, 1800, 1500, 1226, 1400],
               [1700, 2600, 1700, 1400, 1626],
               [2826, 1500, 1500, 1500, 1700],
               [3200, 5826]],
}[which]
p = f"{d}/word/document.xml"
tree = etree.parse(p)
root = tree.getroot()
for sect in root.iter(q("sectPr")):
    pg = sect.find("w:pgSz", NS)
    pg.set(q("w"), "11906"); pg.set(q("h"), "16838")
tables = list(root.iter(q("tbl")))
assert len(tables) == len(WIDTHS), (len(tables), len(WIDTHS))
for tbl, widths in zip(tables, WIDTHS):
    assert sum(widths) == TEXT_W, (widths, sum(widths))
    grid = tbl.findall("w:tblGrid/w:gridCol", NS)
    assert len(grid) == len(widths), (len(grid), len(widths))
    for g, w in zip(grid, widths):
        g.set(q("w"), str(w))
    tw = tbl.find("w:tblPr/w:tblW", NS)
    tw.set(q("w"), str(TEXT_W)); tw.set(q("type"), "dxa")
    for tr in tbl.findall("w:tr", NS):
        col = 0
        for tc in tr.findall("w:tc", NS):
            span_el = tc.find("w:tcPr/w:gridSpan", NS)
            span = int(span_el.get(q("val"))) if span_el is not None else 1
            tcw = tc.find("w:tcPr/w:tcW", NS)
            tcw.set(q("w"), str(sum(widths[col:col + span]))); tcw.set(q("type"), "dxa")
            col += span
        assert col == len(widths), "row does not fill the grid"
# Paper 1's Table 1 has eight columns: 8.5 pt text so that species names and headers do not break mid-word
if which == "paper1":
    for r in tables[0].iter(q("r")):
        rpr = r.find("w:rPr", NS)
        if rpr is None:
            rpr = etree.SubElement(r, q("rPr")); r.remove(rpr); r.insert(0, rpr)
        for tag in ("sz", "szCs"):
            el = rpr.find(f"w:{tag}", NS)
            if el is None:
                el = etree.SubElement(rpr, q(tag))
            el.set(q("val"), "17")
removed = 0
if which == "paper2":
    body = root.find("w:body", NS)
    for el in list(body):
        txt = "".join(t.text or "" for t in el.iter(q("t")))
        if txt.startswith("Figure 1 · the runIEM_HH protocol and its five weak points") and el.find(".//w:drawing", NS) is None:
            body.remove(el); removed += 1
tree.write(p, xml_declaration=True, encoding="UTF-8", standalone=True)
print(which, "A4; tables re-gridded:", len(tables), "; duplicate figure titles removed:", removed)
