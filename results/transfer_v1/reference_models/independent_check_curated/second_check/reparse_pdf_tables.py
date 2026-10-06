"""Independent re-reading of the iAH991 supplement tables (Heinken et al. 2013, Gut Microbes 4:28-40).

The rebuild (scripts/rebuild_iAH991_from_pdf.py) used pdfplumber's extract_tables and line-joining rules. This
script reads the same pages differently, so that its errors would not be the rebuild's:

  * table cells are taken from the drawn cell borders (thin filled rectangles), not from extract_tables;
  * every character is assigned to a cell by its centre, so nothing crosses a column;
  * each cell is returned twice: its characters with every space removed ("despaced"), and its lines.

Comparisons against the model then use the despaced text, which side-steps the question of where a line break
falls inside a token. A row whose bottom border is missing (cut off by the page edge) runs to the page bottom and is
flagged.

The PDF is untrusted input: run with python -I, from this directory, with the PDF path as an argument.

Usage: python -I reparse_pdf_tables.py <pdf> <out.json> [first-last ...]   (default: S10a 36-262, S10b 263-306,
       S12 472-490, S8 and S3 pages are read by other scripts)
"""
from __future__ import annotations

import json
import sys

import pdfplumber


def borders(page):
    tall = sorted({round(r["x0"], 1) for r in page.rects if (r["bottom"] - r["top"]) > 3 and (r["x1"] - r["x0"]) < 2})
    wide = sorted({round(r["top"], 1) for r in page.rects if (r["x1"] - r["x0"]) > 3 and (r["bottom"] - r["top"]) < 2})
    # merge borders closer than 1 pt
    def merge(v):
        out = []
        for x in v:
            if not out or x - out[-1] > 1.0:
                out.append(x)
        return out
    return merge(tall), merge(wide)


def cell_text(chars):
    """(despaced, lines) for the characters of one cell."""
    chars = [c for c in chars if c["text"].strip()]
    if not chars:
        return "", []
    chars.sort(key=lambda c: (round(c["top"]), c["x0"]))
    lines, cur, cur_top = [], [], None
    for c in chars:
        if cur_top is None or abs(c["top"] - cur_top) > 2.0:
            if cur:
                lines.append(cur)
            cur, cur_top = [c], c["top"]
        else:
            cur.append(c)
    if cur:
        lines.append(cur)
    out_lines = []
    for ln in lines:
        ln.sort(key=lambda c: c["x0"])
        s, prev = "", None
        for c in ln:
            if prev is not None and c["x0"] - prev["x1"] > 0.18 * max(c["size"], 1.0):
                s += " "
            s += c["text"]
            prev = c
        out_lines.append(s)
    despaced = "".join("".join(c["text"] for c in sorted(ln, key=lambda c: c["x0"])) for ln in lines).replace(" ", "")
    return despaced, out_lines


def read_pages(pdf, first, last):
    rows = []
    for pn in range(first, last + 1):
        page = pdf.pages[pn - 1]
        xs, ys = borders(page)
        if len(xs) < 3 or len(ys) < 1:
            rows.append({"page": pn, "error": f"borders not found (x {len(xs)}, y {len(ys)})"})
            continue
        chars = [c for c in page.chars if xs[0] - 1 <= (c["x0"] + c["x1"]) / 2 <= xs[-1] + 1]
        below = [c for c in chars if (c["top"] + c["bottom"]) / 2 > ys[-1] and c["text"].strip()]
        bounds = list(ys) + ([page.height] if below else [])
        for i in range(len(bounds) - 1):
            top, bot = bounds[i], bounds[i + 1]
            rc = [c for c in chars if top < (c["top"] + c["bottom"]) / 2 < bot]
            if not rc:
                continue
            cells = []
            for j in range(len(xs) - 1):
                cc = [c for c in rc if xs[j] < (c["x0"] + c["x1"]) / 2 < xs[j + 1]]
                d, l = cell_text(cc)
                cells.append({"despaced": d, "lines": l})
            rows.append({"page": pn, "row_on_page": i, "top": top, "bottom": bot, "open_bottom": bot == page.height,
                         "column_borders": xs, "cells": cells})
    return rows


def main():
    pdf = pdfplumber.open(sys.argv[1])
    out = sys.argv[2]
    ranges = sys.argv[3:] or ["36-262", "263-306", "472-490"]
    res = {}
    for r in ranges:
        a, b = (int(x) for x in r.split("-"))
        res[r] = read_pages(pdf, a, b)
        print(r, "rows", len(res[r]), flush=True)
    with open(out, "w") as fh:
        json.dump(res, fh)


if __name__ == "__main__":
    main()
