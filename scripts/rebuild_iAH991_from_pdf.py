#!/usr/bin/env python3
"""Rebuild iAH991 from the supplementary PDF of Heinken, Sahoo, Fleming & Thiele 2013.

iAH991 is the genome-scale metabolic reconstruction of Bacteroides thetaiotaomicron
VPI-5482 (Gut Microbes 4:28-40, doi:10.4161/gmic.22370). No SBML file was published, so the
model is assembled from the supplementary PDF:

  * Table S10a (all 1,488 reactions; pp. 36-262) and Table S10b (all metabolites; pp. 263-306)
    are the model source.
  * Table S12 (the joint B. thetaiotaomicron/mouse model, pp. 472-491) lists every iAH991
    reaction a second time with a "BT" prefix. It is used (i) as an independent rendering to
    cross-check the S10a parse and (ii) only where S10a itself is unreadable; each such use is
    listed as a manual intervention.

The rebuilt model is then validated against the paper's own published model predictions
(Tables S1, S3, S8, S9, S14 and the gene-essentiality counts in the supplementary text). The
script never reads experimental data other than what is printed in the PDF.

Usage:
    python3 -I scripts/rebuild_iAH991_from_pdf.py <pdf> <outdir>

Outputs in <outdir>: iAH991_rebuilt.xml, rebuild_report.json, REBUILD.md
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import html
import itertools
import json
import math
import platform
import re
import sys
import time

import cobra
import pdfplumber
from cobra.flux_analysis import flux_variability_analysis, single_gene_deletion

# --------------------------------------------------------------------------------------
# Constants (all page numbers are 1-based PDF pages; every range is verified at run time)
# --------------------------------------------------------------------------------------
EXPECTED_SHA256 = "49502d765c7d049208aeb03a140aa7973d4569b7cddf22b9ca2ea4fa019975b9"
S10A_PAGES = (36, 262)
S10B_PAGES = (263, 306)
S10CG_PAGES = (307, 402)  # Tables S10c-S10g: only used as a word list for text columns
S12_PAGES = (472, 491)  # Table S12, B. thetaiotaomicron part (mouse rows start on p. 490)
S1_PAGES = (13, 14)
S3_PAGES = (17, 18)
S8_PAGES = {"S8a": 30, "S8b": 31, "S8c": 32, "S8d": 33}
S9_PAGE = 35
S14_PAGE = 21
TEXT_PAGE_ESSENTIALITY = 8
PROSE_PAGES = (2, 27)  # supplementary text, figure legends, small tables, references (word list)
# vitamins with an iAH991 exchange reaction that Table S8d (TYG) does not list; used only for an
# exploratory sensitivity table of the TYG essential-gene count
TYG_VITAMIN_PROBES = ["btn", "pnto-R", "nac", "fol"]

S10A_HEADER = ["SEED Rxn ID", "BIGG Rxn ID", "Function", "Formula", "Rev", "GA", "LB", "UB",
               "CS", "Subsystem", "Reference", "Notes", "EC\nNumber"]
S10B_HEADER = ["SEED Met ID", "BIGG Met ID", "Name", "Formula", "Charge", "KEGG ID"]

PUBLISHED_COUNTS = {"reactions": 1488, "metabolic_and_transport": 1213,
                    "exchange_and_demand": 275, "metabolites": 1152, "genes": 991}

GROWTH_EPS = 1e-6  # growth above this counts as "grows" (published values go down to 0.0012)
SOLVER = "glpk"

ARROWS = ("<=>", "=>", "->", "<=", "<-")
FORMULA_COMPLETE = re.compile(
    r"^(?:\+|<=>|=>|->|<=|<-|\(?\d+(?:\.\d+)?\)?|[^\s\[\]]+\[[a-z]\])$")
COEF_RE = re.compile(r"^\(?(\d+(?:\.\d+)?)\)?$")
MET_TOKEN_RE = re.compile(r"^([^\s\[\]]+)\[([a-z])\]$")
GENE_RE = re.compile(r"^BT_\d{4}$")
GPR_WORD_COMPLETE = re.compile(r"^(?:and|or|BT_\d{4})?$")
INT_RE = re.compile(r"^-?\d+$")
AMBIGUOUS_NUM_RE = re.compile(r"^-?\d{1,3}\.\d{3}$")  # "1.000": decimal point or 1000s separator?
WORD_RE = re.compile(r"[A-Za-z]+")
EC_RE = re.compile(r"^\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:n?\d+|-)$")


def log(msg):
    print(f"[rebuild] {msg}", file=sys.stderr, flush=True)


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------------------
# PDF table extraction
# --------------------------------------------------------------------------------------
def extract_page_tables(page):
    """Tables on one page.

    An explicit horizontal edge 0.25 pt above the page bottom closes table rows whose cell
    borders run past the page edge (rows clipped by the typesetting). pdfplumber drops such
    rows otherwise. The edge only intersects vertical cell borders that reach the page
    bottom, so ordinary pages are unaffected (checked in the report: row counts per page are
    compared with and without the edge on every clipped page).
    """
    return page.extract_tables({"explicit_horizontal_lines": [page.height - 0.25]})


def page_clip(page):
    """How far (pt) cell rectangles and characters extend below the page bottom."""
    rb = max((r["bottom"] for r in page.rects), default=0.0)
    cb = max((c["bottom"] for c in page.chars), default=0.0)
    return max(0.0, rb - page.height), max(0.0, cb - page.height)


def clean_cell(c):
    return "" if c is None else c


def extract_header_table(pdf, first, last, header, name):
    """Rows of a table whose header row is repeated on every page in [first, last]."""
    rows, pages_info = [], []
    for p in range(first, last + 1):
        page = pdf.pages[p - 1]
        tables = extract_page_tables(page)
        if len(tables) != 1:
            raise RuntimeError(f"{name}: page {p} has {len(tables)} tables, expected 1")
        t = [[clean_cell(c) for c in r] for r in tables[0]]
        if t[0] != header:
            raise RuntimeError(f"{name}: page {p} header {t[0]!r} differs from {header!r}")
        data = t[1:]
        rect_clip, char_clip = page_clip(page)
        info = {"page": p, "rows": len(data), "rect_clip_pt": round(rect_clip, 1),
                "char_clip_pt": round(char_clip, 1)}
        if rect_clip > 0.5:
            default = page.extract_tables()
            info["rows_without_bottom_edge"] = len(default[0]) - 1 if default else 0
        pages_info.append(info)
        for i, r in enumerate(data, start=1):
            if len(r) != len(header):
                raise RuntimeError(f"{name}: page {p} row {i} has {len(r)} cells")
            rows.append({"page": p, "row": i, "cells": r,
                         "clipped": rect_clip > 0.5 and i == len(data)})
    # the pages just outside the range must not carry the same header
    for p in (first - 1, last + 1):
        tables = extract_page_tables(pdf.pages[p - 1])
        if tables and [clean_cell(c) for c in tables[0][0]] == header:
            raise RuntimeError(f"{name}: page {p} outside the expected range has the header")
    return rows, pages_info


def extract_plain_tables(pdf, first, last):
    """All cells of all tables on pages [first, last] (used for the word list only)."""
    cells = []
    for p in range(first, last + 1):
        for t in extract_page_tables(pdf.pages[p - 1]):
            for r in t:
                cells.extend(clean_cell(c) for c in r)
    return cells


def extract_s12_bt_rows(pdf, first, last):
    """Table S12 rows that belong to B. thetaiotaomicron (Type column)."""
    first_text = pdf.pages[first - 1].extract_text() or ""
    if "Table S12" not in first_text:
        raise RuntimeError(f"S12: 'Table S12' not found on page {first}")
    rows, last_bt_page, mouse_seen_on = [], None, None
    for p in range(first, last + 1):
        page = pdf.pages[p - 1]
        rect_clip, _ = page_clip(page)
        if rect_clip > 0.5:
            raise RuntimeError(f"S12: page {p} has clipped rows; not handled")
        for t in extract_page_tables(page):
            for i, r in enumerate(t, start=1):
                r = [clean_cell(c) for c in r]
                if len(r) != 6:
                    continue
                typ = r[5].replace("\n", " ")
                if "thetaiotaomicron" in typ:
                    rows.append({"page": p, "row": i, "cells": r})
                    last_bt_page = p
                elif "Mouse" in typ or "mouse" in typ:
                    mouse_seen_on = mouse_seen_on or p
    if mouse_seen_on is None or last_bt_page is None or last_bt_page >= last:
        raise RuntimeError("S12: could not confirm that all B. theta rows were scanned")
    return rows, {"pages": [first, last_bt_page], "first_mouse_page": mouse_seen_on}


# --------------------------------------------------------------------------------------
# Joining wrapped cell fragments
# --------------------------------------------------------------------------------------
def split_frags(cell):
    return [f.strip() for f in (cell or "").split("\n") if f.strip()]


def join_nospace(cell):
    """Identifiers and numbers never contain spaces: a line break is always inside the token."""
    return "".join(split_frags(cell))


def join_formula(cell):
    """Join formula fragments. A fragment whose last token is complete (a metabolite with its
    [compartment], a stoichiometric coefficient, '+' or an arrow) was broken at a space;
    otherwise the break is inside a token and the fragments are joined without a space."""
    frags = split_frags(cell)
    breaks = []
    if not frags:
        return "", breaks
    out = frags[0]
    for nxt in frags[1:]:
        last = out.split()[-1]
        first = nxt.split()[0]
        if FORMULA_COMPLETE.match(last):
            breaks.append({"kind": "space", "left": last, "right": first})
            out = f"{out} {nxt}"
        else:
            breaks.append({"kind": "joined", "left": last, "right": first})
            out = out + nxt
    return out, breaks


def join_gpr(cell):
    """Join gene-rule fragments; same idea as join_formula with tokens BT_dddd/and/or/()."""
    frags = split_frags(cell)
    breaks = []
    if not frags:
        return "", breaks
    out = frags[0]
    for nxt in frags[1:]:
        last = out.split()[-1].strip("()")
        if GPR_WORD_COMPLETE.match(last):
            breaks.append({"kind": "space", "left": out.split()[-1], "right": nxt.split()[0]})
            out = f"{out} {nxt}"
        else:
            breaks.append({"kind": "joined", "left": out.split()[-1], "right": nxt.split()[0]})
            out = out + nxt
    return out, breaks


def interior_words(cell):
    """Words that cannot be a fragment of a word broken across lines."""
    frags = (cell or "").split("\n")
    n = len(frags)
    words = []
    for k, frag in enumerate(frags):
        ws = [(m.group(0), m.start(), m.end()) for m in WORD_RE.finditer(frag)]
        for w, s, e in ws:
            if k < n - 1 and e == len(frag.rstrip()) and frag.rstrip():
                continue  # touches the end of a non-final fragment
            if k > 0 and s == len(frag) - len(frag.lstrip()):
                continue  # touches the start of a non-first fragment
            words.append(w.lower())
    return words


class TextJoiner:
    """Join wrapped fragments of free-text cells (names, subsystems, notes, references).

    Rules: no space after a trailing '-', '/' or '(' or before ':', ')', ',' or '.';
    a word broken without a hyphen (w1 | w2) is re-joined when w1+w2 occurs as a word somewhere
    in Tables S10a-S10g or the supplementary text (pp. 2-27) and w1 or w2 is never seen as a
    complete word (i.e. away from a line break); otherwise a space is inserted.
    These columns are annotation only and do not affect the model's mathematics.
    """

    def __init__(self, cells, prose_lines=()):
        self.vocab = set()  # words seen away from any line break: certainly complete words
        self.seen = set()  # every word seen anywhere
        for c in cells:
            self.vocab.update(interior_words(c))
            self.seen.update(w.lower() for w in WORD_RE.findall(c or ""))
        for line in prose_lines:
            ws = [w.lower() for w in WORD_RE.findall(line)]
            self.seen.update(ws)
            self.vocab.update(ws[1:-1])
        self.stats = collections.Counter()

    def join(self, cell):
        frags = split_frags(cell)
        if not frags:
            return ""
        out = frags[0]
        for nxt in frags[1:]:
            if out.endswith(("-", "/", "(")) or nxt.startswith((":", ")", ",", ".", ";")):
                out += nxt
                self.stats["punctuation_join"] += 1
                continue
            m1 = re.search(r"[A-Za-z]+$", out)
            m2 = re.match(r"^[A-Za-z]+", nxt)
            if m1 and m2:
                w1, w2 = m1.group(0).lower(), m2.group(0).lower()
                if (w1 + w2) in self.seen and (w1 not in self.vocab or w2 not in self.vocab):
                    out += nxt
                    self.stats["word_rejoined"] += 1
                    continue
                if w1 not in self.vocab or w2 not in self.vocab:
                    self.stats["space_unverified"] += 1
                else:
                    self.stats["space"] += 1
            else:
                self.stats["space"] += 1
            out += " " + nxt
        return out


def join_ec(cell):
    frags = split_frags(cell)
    if not frags:
        return ""
    out = frags[0]
    for nxt in frags[1:]:
        out = out + nxt if out.endswith(".") else f"{out} {nxt}"
    return out


# --------------------------------------------------------------------------------------
# Formula / gene-rule parsing
# --------------------------------------------------------------------------------------
class ParseError(Exception):
    pass


def parse_side(tokens):
    terms, i, expect_term = [], 0, True
    while i < len(tokens):
        t = tokens[i]
        if expect_term:
            coef = 1.0
            m = COEF_RE.match(t)
            if m:
                coef = float(m.group(1))
                i += 1
                if i >= len(tokens):
                    raise ParseError(f"coefficient {t!r} without metabolite")
                t = tokens[i]
            mm = MET_TOKEN_RE.match(t)
            if not mm:
                raise ParseError(f"bad metabolite token {t!r}")
            terms.append((mm.group(1), mm.group(2), coef))
            i += 1
            expect_term = False
        else:
            if t != "+":
                raise ParseError(f"expected '+', got {t!r}")
            i += 1
            expect_term = True
    if terms and expect_term:
        raise ParseError("dangling '+'")
    return terms


def parse_formula(text):
    """Return (lhs, rhs, arrow); lhs/rhs are lists of (met, compartment, coefficient)."""
    tokens = text.split()
    idx = [i for i, t in enumerate(tokens) if t in ARROWS]
    if len(idx) != 1:
        raise ParseError(f"{len(idx)} arrows")
    i = idx[0]
    arrow = tokens[i]
    lhs, rhs = parse_side(tokens[:i]), parse_side(tokens[i + 1:])
    if arrow in ("<=", "<-"):
        lhs, rhs, arrow = rhs, lhs, "=>"
    if arrow == "->":
        arrow = "=>"
    return lhs, rhs, arrow


def stoichiometry(lhs, rhs):
    st = collections.defaultdict(float)
    for m, c, k in lhs:
        st[(m, c)] -= k
    for m, c, k in rhs:
        st[(m, c)] += k
    return dict(st)


def gpr_tokens(rule):
    return re.findall(r"\(|\)|[^\s()]+", rule)


def check_gpr(tokens):
    """Recursive-descent syntax check; returns an error string or None."""
    pos = 0

    def term():
        nonlocal pos
        if pos >= len(tokens):
            raise ParseError("unexpected end")
        t = tokens[pos]
        if t == "(":
            pos += 1
            expr()
            if pos >= len(tokens) or tokens[pos] != ")":
                raise ParseError("missing ')'")
            pos += 1
        elif GENE_RE.match(t):
            pos += 1
        else:
            raise ParseError(f"unexpected token {t!r}")

    def expr():
        nonlocal pos
        term()
        while pos < len(tokens) and tokens[pos] in ("and", "or"):
            pos += 1
            term()

    if not tokens:
        return None
    try:
        expr()
        if pos != len(tokens):
            raise ParseError(f"trailing token {tokens[pos]!r}")
    except ParseError as exc:
        return str(exc)
    return None


def normalize_gpr(tokens):
    return " ".join(tokens).replace("( ", "(").replace(" )", ")")


def gpr_dnf(tokens):
    """Set of frozensets of genes (disjunctive normal form) for logical comparison."""
    pos = 0

    def term():
        nonlocal pos
        t = tokens[pos]
        if t == "(":
            pos += 1
            r = expr()
            pos += 1
            return r
        pos += 1
        return {frozenset([t])}

    def expr():
        nonlocal pos
        result = term()
        while pos < len(tokens) and tokens[pos] in ("and", "or"):
            op = tokens[pos]
            pos += 1
            rhs = term()
            if op == "or":
                result = result | rhs
            else:
                result = {a | b for a in result for b in rhs}
        return result

    if not tokens:
        return set()
    dnf = expr()
    return {c for c in dnf if not any(o < c for o in dnf)}


def parse_number(txt):
    if INT_RE.match(txt):
        return float(int(txt)), None
    if AMBIGUOUS_NUM_RE.match(txt):
        return None, "ambiguous_number_format"
    try:
        return float(txt.replace(",", ".")), "non_integer"
    except ValueError:
        return None, "unparseable"


# --------------------------------------------------------------------------------------
# Identifier conversion
# --------------------------------------------------------------------------------------
def cobra_met_id(mid, comp):
    s = mid.replace("-", "__")
    s = re.sub(r"[^A-Za-z0-9_]", "_", s)
    return f"{s}_{comp}"


def cobra_rxn_id(rid):
    s = rid.replace("(e)", "_e").replace("-", "__")
    return re.sub(r"[^A-Za-z0-9_]", "_", s)


def esc(s):
    return html.escape(str(s), quote=False)


class OrderedGroup(cobra.core.Group):
    """cobra Group whose members iterate in id order.

    cobra stores group members in a set, so cobra.io.write_sbml_model would write them in a
    hash-dependent order and the SBML file would differ between runs.
    """

    @property
    def members(self):
        return sorted(self._members, key=lambda x: x.id)


# --------------------------------------------------------------------------------------
# Small tables printed in the PDF (simulation conditions and published predictions)
# --------------------------------------------------------------------------------------
NUM_TOKEN = r"\d+(?:[.,]\d+)?"


def parse_rate_table(text, known_ids):
    """Lines '<name> <BIGG id> <rate>' (rate may use a decimal comma)."""
    entries = []
    for line in text.split("\n"):
        m = re.match(rf"^(?:(?P<name>.*?)\s+)?(?P<id>\S+)\s+(?P<rate>{NUM_TOKEN})\s*$",
                     line.strip())
        if m and m.group("id") in known_ids:
            entries.append({"name": (m.group("name") or "").strip(), "id": m.group("id"),
                            "rate": float(m.group("rate").replace(",", "."))})
    return entries


def parse_s8(pdf, known_ids):
    out = {}
    expected_titles = {"S8a": "a) Essential minerals and cofactors",
                       "S8b": "b) Uptake rates for calculations of iAH991 growth",
                       "S8c": "c) Simulation constraints for gene essentiality analysis of iAH991 "
                              "corresponding to glucose minimal medium",
                       "S8d": "d) Simulation constraints for gene essentiality analysis of iAH991 "
                              "corresponding to anaerobic tryptone-yeast extract-glucose medium"}
    for key, p in S8_PAGES.items():
        text = pdf.pages[p - 1].extract_text() or ""
        flat = re.sub(r"\s+", " ", text)
        if expected_titles[key] not in flat:
            raise RuntimeError(f"{key}: title not found on page {p}")
        out[key] = {"page": p, "entries": parse_rate_table(text, known_ids)}
    return out


def parse_s9_western(pdf, known_ids, exchange_ids):
    """Table S9 (diets of the joint model): name, BIGG id, then 5 rates (HP, HC, HF, KD, W)."""
    p = S9_PAGE
    text = pdf.pages[p - 1].extract_text() or ""
    if "Table S9" not in text or "Western diet" not in text:
        raise RuntimeError("S9: title not found")
    entries = []
    for line in text.split("\n"):
        toks = line.strip().split()
        if len(toks) < 7:
            continue
        nums = toks[-5:]
        if not all(re.fullmatch(NUM_TOKEN, t) for t in nums):
            continue
        bid = toks[-6]
        if bid not in known_ids and bid not in exchange_ids:
            continue
        rates = [float(t.replace(",", ".")) for t in nums]
        entries.append({"name": " ".join(toks[:-6]), "id": bid, "western": rates[4],
                        "all_diets": dict(zip(["HP", "HC", "HF", "KD", "W"], rates))})
    return {"page": p, "entries": entries}


def column_words(pdf, pages, x_lo, x_hi, y_start=None, y_stop_word=None):
    """Words whose x0 lies in [x_lo, x_hi) in reading order over several pages."""
    out = []
    for p in pages:
        page = pdf.pages[p - 1]
        words = page.extract_words()
        stop_y = None
        if y_stop_word:
            for w in words:
                if w["text"] == y_stop_word[0] and w["x0"] < 100:
                    stop_y = w["top"]
                    break
        for w in sorted(words, key=lambda w: (round(w["top"]), w["x0"])):
            if y_start and p == pages[0] and w["top"] < y_start:
                continue
            if stop_y is not None and w["top"] >= stop_y:
                continue
            if x_lo <= w["x0"] < x_hi:
                out.append((p, round(w["top"]), w["text"]))
    return out


# Table S3a (pp. 17-18): carbon source -> S8b uptake ids, published iAH991 ("in silico") growth.
# The transcription is checked against the in-silico column of the PDF at run time.
S3_ROWS = [
    ("Arabinose", {"arab-L": None}, "0.23"),
    ("Fructose", {"fru": None}, "0.24"),
    ("Fucose", {"fuc-L": None}, "0.0012"),
    ("Galactose", {"gal": None}, "0.24"),
    ("Galacturonate", {"galur": None}, "0.14"),
    ("Glucosamine", {"gam": None}, "0.24"),
    ("Glucose", {"glc-D": None}, "0.24"),
    ("Glucose + 1 mmol/gDW/h oxygen", {"glc-D": None, "o2": 1.0}, "0.27"),
    ("Glucuronate", {"glcur": None}, "0.14"),
    ("Mannose", {"man": None}, "0.16"),
    ("N-acetylgalactosamine", {"acgal": None}, "0.17"),
    ("N-acetylglucosamine", {"acgam": None}, "0.18"),
    ("N-acetylneuraminate", {"acnam": None}, "0"),
    ("Rhamnose", {"rmn": None}, "0.0012"),
    ("Ribose", {"rib-D": None}, "0.23"),
    ("Xylose", {"xyl-D": None}, "0.23"),
    ("Arabinan", {"arabinan101": None}, "0.23"),
    ("Arabinogalactan", {"arabinogal": None}, "0.23"),
    ("Homogalacturonan", {"homogal": None}, "0.12"),
    ("Pectic galactan", {"pecticgal": None}, "0.22"),
    ("Rhamnogalacturonan I", {"rhamnogalurI": None}, "0.13"),
    ("Rhamnogalacturonan II", {"rhamnogalurII": None}, "0.08"),
    ("Amylopectin", {"strch1": None}, "0.25"),
    ("Pullulan", {"pullulan1200": None}, "0.25"),
    ("Dextran", {"dextran40": None}, "0.26"),
    ("Inulin", {"inulin": None}, "0.24"),
    ("Levan", {"levan1000": None}, "0.24"),
    ("F1alpha (mucin O-glycan)", {"f1a": None}, "0.15"),
    ("Core 4", {"core4": None}, "0.12"),
    ("Core 5", {"core5": None}, "0.11"),
    ("Core 7", {"core7": None}, "0.11"),
    ("Core 8", {"core8": None}, "0.17"),
    ("GlcNAc-alpha-1,4-Core 1", {"gncore1": None}, "0.15"),
    ("GlcNAc-alpha-1,4-Core 2", {"gncore2": None}, "0.14"),
    ("Sialyl-Tn antigen", {"sTn_antigen": None}, "0"),
    ("Disialyl-T antigen", {"dsT_antigen": None}, "0.01"),
    ("alpha-mannan", {"amannan140": None}, "0.16"),
    ("Chondroitin sulfate A", {"cspg_a": None}, "0.05"),
    ("Chondroitin sulfate B", {"cspg_b": None}, "0.02"),
    ("Chondroitin sulfate (C)", {"cspg_c": None}, "0.05"),
    ("Hyaluronan", {"ha": None}, "0.13"),
    ("Heparin", {"hspg": None}, "0.11"),
    ("Glycogen", {"glycogen1500": None}, "0.25"),
    ("Laminarin", {"lmn30": None}, "0.23"),
    ("Cellobiose", {"cellb": None}, "0"),
]
# S3 has a single "Amylopectin" row but S8b lists two starch substrates (starch1200 and strch1);
# the alternative is reported as an extra row.
S3_ALTERNATIVES = [("Amylopectin (alternative mapping: starch1200)", {"starch1200": None}, "0.25")]

# Table S3b (p. 18): in silico secretion at 10 mmol/gDW/h substrate uptake (reported as ranges).
S3B_ROWS = [
    ("Glucose", "glc-D", {"ac": "5.56", "ppa": "0 to 8.95", "succ": "0 to 8.95",
                          "co2": "-1.53 to 7.41", "h2": "2.10"}),
    ("N-acetylglucosamine", "acgam", {"ac": "8.10 to 12.63", "ppa": "0 to 14.56",
                                      "succ": "0 to 14.56", "co2": "-5.82 to 8.75",
                                      "h2": "< 0.001"}),
    ("Glucuronate", "glcur", {"ac": "7.72 to 8.16", "ppa": "0 to 7.66", "succ": "0 to 7.66",
                              "co2": "7.60 to 14.66", "h2": "< 0.001"}),
    ("Fructose", "fru", {"ac": "5.56", "ppa": "0 to 8.95", "succ": "0 to 8.95",
                         "co2": "-1.53 to 7.41", "h2": "2.10"}),
]

# Table S14 (p. 21): knockout strains and the published iAH991 ("in silico") growth call.
# Each case: deleted gene sets (each set is one strain), substrates of the condition, call.
S14_ROWS = [
    ("BT_3763-BT_3767", [["BT_3763", "BT_3764", "BT_3765", "BT_3766", "BT_3767"]],
     "L-rhamnose as sole carbon source", [{"rmn": None}], "no"),
    ("BT_3763-BT_3767", [["BT_3763", "BT_3764", "BT_3765", "BT_3766", "BT_3767"]],
     "Glucose minimal medium", [{"glc-D": None}], "yes"),
    ("BT_3767", [["BT_3767"]], "L-rhamnose as sole carbon source", [{"rmn": None}], "no"),
    ("BT_3717", [["BT_3717"]], "Glucose minimal medium", [{"glc-D": None}], "no"),
    ("BT_3717", [["BT_3717"]], "Glucose minimal medium with L-arginine supply",
     [{"glc-D": None, "arg-L": 10.0}], "yes"),
    ("SusC (BT_3702) or SusD (BT_3701) or SusG (BT_3698)",
     [["BT_3702"], ["BT_3701"], ["BT_3698"]],
     "Maltose or maltotriose as sole carbon source", [{"malt": None}, {"malttr": None}], "yes"),
    ("SusC (BT_3702) or SusD (BT_3701) or SusG (BT_3698)",
     [["BT_3702"], ["BT_3701"], ["BT_3698"]],
     "Starch as sole carbon source", [{"starch1200": None}], "no"),
    ("SusG (BT_3698)", [["BT_3698"]], "Maltoheptaose as sole carbon source",
     [{"malthp": None}], "yes"),
    ("SusC (BT_3702) or SusD (BT_3701)", [["BT_3702"], ["BT_3701"]],
     "Maltoheptaose as sole carbon source", [{"malthp": None}], "no"),
    ("BT_1663 or BT_4689", [["BT_1663"], ["BT_4689"]], "Pullulan as sole carbon source",
     [{"pullulan1200": None}], "yes"),
    ("BT_0238", [["BT_0238"]], "Chondroitin sulfate A, B or C as sole carbon source",
     [{"cspg_a": None}, {"cspg_b": None}, {"cspg_c": None}], "no"),
    ("BT_0238", [["BT_0238"]], "Glucose minimal medium", [{"glc-D": None}], "yes"),
    ("csuF (BT_3332)", [["BT_3332"]], "Chondroitin sulfate A, B or C as sole carbon source",
     [{"cspg_a": None}, {"cspg_b": None}, {"cspg_c": None}], "no"),
    ("csuF (BT_3332)", [["BT_3332"]], "Sulfated disaccharides as sole carbon source",
     [{"cspg_a_degr": None}, {"cspg_b_degr": None}, {"cspg_c_degr": None}], "yes"),
    ("BT_2970", [["BT_2970"]], "Glucose minimal medium", [{"glc-D": None}], "yes"),
    ("BT_1760", [["BT_1760"]], "Levan as sole carbon source", [{"levan1000": None}], "no"),
    ("BT_1760", [["BT_1760"]], "Glucose, fructose, sucrose or oligofructose as sole carbon source",
     [{"glc-D": None}, {"fru": None}, {"sucr": None}, {"oligofru4": None}], "yes"),
]

S1_SECTION_HEADERS = {"Monosaccharides", "Di-and oligosaccharides", "Polysaccharides"}


def parse_s1(pdf):
    """Table S1 from word positions: name (x<190), BIGG id(s) (190-285), BΘ_Seed_v1 (285-330),
    iAH991 (360-400). A row closes on the line carrying the iAH991 yes/no; id and name words
    since the previous row belong to it (group rows list several ids on preceding lines)."""
    rows, ids, names = [], [], []
    for p in range(S1_PAGES[0], S1_PAGES[1] + 1):
        page = pdf.pages[p - 1]
        words = page.extract_words()
        header = [w for w in words if w["text"] == "iAH991" and 350 < w["x0"] < 400]
        y0 = header[0]["bottom"] if header and p == S1_PAGES[0] else 0
        lines = collections.defaultdict(list)
        for w in words:
            if w["top"] > y0 and w["top"] < page.height - 60:
                lines[round(w["top"])].append(w)
        for y in sorted(lines):
            ws = sorted(lines[y], key=lambda w: w["x0"])
            name_ws = [w["text"] for w in ws if w["x0"] < 190]
            id_ws = [w["text"] for w in ws if 190 <= w["x0"] < 285]
            seed_ws = [w["text"] for w in ws if 285 <= w["x0"] < 330]
            iah_ws = [w["text"] for w in ws if 360 <= w["x0"] < 400]
            line_name = " ".join(name_ws)
            if line_name in S1_SECTION_HEADERS and not id_ws and not iah_ws:
                continue
            names.extend(name_ws)
            for t in id_ws:
                ids.extend(x for x in t.split(",") if x)
            if iah_ws:
                rows.append({"page": p, "name": " ".join(names), "ids": ids,
                             "seed_v1": " ".join(seed_ws), "iAH991": iah_ws[0]})
                ids, names = [], []
    return rows


# --------------------------------------------------------------------------------------
# Simulation helpers
# --------------------------------------------------------------------------------------
class Sim:
    def __init__(self, model, ex_by_met, mets_info):
        self.model = model
        self.ex_by_met = ex_by_met
        self.ex_rxns = [model.reactions.get_by_id(r) for r in sorted(set(ex_by_met.values()))]
        self.mets_info = mets_info

    def apply(self, uptake):
        """Close every exchange uptake, then open the listed ones (lower bound = -rate)."""
        missing = [m for m in uptake if m not in self.ex_by_met]
        for r in self.ex_rxns:
            r.lower_bound = 0.0
        for m, rate in uptake.items():
            if m in self.ex_by_met:
                self.model.reactions.get_by_id(self.ex_by_met[m]).lower_bound = -float(rate)
        return missing

    def growth(self, uptake, knockouts=()):
        with self.model:
            missing = self.apply(uptake)
            for g in knockouts:
                if g in self.model.genes:
                    self.model.genes.get_by_id(g).knock_out()
            val = self.model.slim_optimize(error_value=float("nan"))
        return (0.0 if (val is None or math.isnan(val)) else float(val)), missing


def hexose_equivalent_rate(met_formula):
    """Uptake giving 60 mmol C/gDW/h (10 mmol hexose units), the scaling used in Table S3."""
    m = re.findall(r"C(\d*)(?![a-z])", met_formula or "")
    c = sum(int(x) if x else 1 for x in m)
    return (60.0 / c) if c else None


def fmt(x, nd=4):
    if x is None:
        return ""
    if isinstance(x, float) and math.isnan(x):
        return "nan"
    if abs(x) < 0.5 * 10 ** (-nd):
        x = 0.0  # no "-0.0000"
    return f"{x:.{nd}f}"


def published_decimals(s):
    return len(s.split(".")[1]) if "." in s else 0


def num(x, nd=6):
    """Round a solver value; -0.0 and values that round to zero become 0.0 (stable output)."""
    v = round(float(x), nd)
    return 0.0 if v == 0 else v


def rounds_to(g, pub):
    """True if g, rounded to the number of decimals printed in pub, equals pub."""
    nd = published_decimals(pub)
    p = float(pub)
    return abs(g - p) <= 0.5 * 10 ** (-nd) + 1e-9


def md_table(headers, rows):
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c).replace("|", "\\|") for c in r) + " |")
    return "\n".join(out)


# --------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("pdf")
    ap.add_argument("outdir")
    ap.add_argument("--processes", type=int, default=2)
    args = ap.parse_args()
    t_start = time.time()

    report = collections.OrderedDict()
    report["tool"] = {"script": "scripts/rebuild_iAH991_from_pdf.py",
                      "python": platform.python_version(), "cobra": cobra.__version__,
                      "pdfplumber": pdfplumber.__version__, "solver": SOLVER}
    pdf_sha = sha256_of(args.pdf)
    report["source"] = {
        "pdf_file": args.pdf.split("/")[-1], "sha256": pdf_sha,
        "sha256_matches_expected": pdf_sha == EXPECTED_SHA256,
        "citation": "Heinken A, Sahoo S, Fleming RMT, Thiele I. Systems-level characterization "
                    "of a host-microbe metabolic symbiosis in the mammalian gut. Gut Microbes "
                    "2013;4(1):28-40. doi:10.4161/gmic.22370 (supplementary material)",
    }
    if pdf_sha != EXPECTED_SHA256:
        log(f"WARNING: sha256 {pdf_sha} differs from the expected {EXPECTED_SHA256}")

    manual = []  # manual interventions
    parse_failures = []

    with pdfplumber.open(args.pdf) as pdf:
        report["source"]["pages"] = len(pdf.pages)
        log("extracting Table S10a")
        s10a_raw, s10a_pages = extract_header_table(pdf, *S10A_PAGES, S10A_HEADER, "S10a")
        log("extracting Table S10b")
        s10b_raw, s10b_pages = extract_header_table(pdf, *S10B_PAGES, S10B_HEADER, "S10b")
        log("extracting Tables S10c-S10g (word list for text columns)")
        s10cg_cells = extract_plain_tables(pdf, *S10CG_PAGES)
        log("extracting Table S12 (B. thetaiotaomicron part)")
        s12_raw, s12_info = extract_s12_bt_rows(pdf, *S12_PAGES)
        log("reading condition tables S1, S3, S8, S9, S14 and the supplementary text")
        prose_lines = []
        for p in range(PROSE_PAGES[0], PROSE_PAGES[1] + 1):
            prose_lines.extend((pdf.pages[p - 1].extract_text() or "").split("\n"))
        text_p8 = pdf.pages[TEXT_PAGE_ESSENTIALITY - 1].extract_text() or ""
        s10c_text = re.sub(r"\s+", " ", pdf.pages[S10CG_PAGES[0] - 1].extract_text() or "")
        s3_col = column_words(pdf, list(range(S3_PAGES[0], S3_PAGES[1] + 1)), 265, 300,
                              y_start=190, y_stop_word=("b)",))
        s14_col = column_words(pdf, [S14_PAGE], 435, 470, y_start=120)
        s1_rows = parse_s1(pdf)
        s10b_ids_for_tables = {join_nospace(r["cells"][1]) for r in s10b_raw}
        s8 = parse_s8(pdf, s10b_ids_for_tables)

        report["extraction"] = {
            "S10a": {"pages": list(S10A_PAGES), "rows": len(s10a_raw),
                     "pages_with_clipped_rows": [i for i in s10a_pages if i["rect_clip_pt"] > 0.5]},
            "S10b": {"pages": list(S10B_PAGES), "rows": len(s10b_raw),
                     "pages_with_clipped_rows": [i for i in s10b_pages if i["rect_clip_pt"] > 0.5]},
            "S12_bt": {**s12_info, "rows": len(s12_raw)},
            "table_settings": "pdfplumber extract_tables with explicit_horizontal_lines="
                              "[page.height - 0.25] (closes rows clipped at the page bottom)",
        }

        # ------------------------------------------------------------------ text joiner
        vocab_cells = [c for r in s10a_raw for k, c in enumerate(r["cells"]) if k in (2, 9, 10, 11)]
        vocab_cells += [r["cells"][2] for r in s10b_raw] + s10cg_cells
        joiner = TextJoiner(vocab_cells, prose_lines)

        # ------------------------------------------------------------------ S10b
        mets_pub = collections.OrderedDict()
        for r in s10b_raw:
            seed, mid, name, formula, charge, kegg = r["cells"]
            mid = join_nospace(mid)
            ch_txt = join_nospace(charge)
            ch = int(ch_txt) if INT_RE.match(ch_txt) else None
            if ch is None:
                parse_failures.append({"table": "S10b", "page": r["page"], "row": r["row"],
                                       "id": mid, "problem": f"charge {ch_txt!r} not an integer"})
            if mid in mets_pub:
                parse_failures.append({"table": "S10b", "page": r["page"], "row": r["row"],
                                       "id": mid, "problem": "duplicate metabolite id"})
            mets_pub[mid] = {"seed": join_nospace(seed), "name": joiner.join(name),
                             "formula": join_nospace(formula), "charge": ch,
                             "kegg": join_nospace(kegg), "page": r["page"], "row": r["row"]}
        known_mets = set(mets_pub)

        # ------------------------------------------------------------------ S10a
        rxns = []
        unreadable = []
        formula_breaks = collections.Counter()
        gpr_breaks = collections.Counter()
        numeric_breaks = []
        naive_bad_tokens = 0
        for r in s10a_raw:
            c = r["cells"]
            rec = {"page": r["page"], "row": r["row"], "clipped": r["clipped"],
                   "seed": join_nospace(c[0]), "id": join_nospace(c[1]),
                   "name": joiner.join(c[2]), "rev": join_nospace(c[4]),
                   "lb_txt": join_nospace(c[6]), "ub_txt": join_nospace(c[7]),
                   "cs": join_nospace(c[8]), "subsystem": joiner.join(c[9]),
                   "references": joiner.join(c[10]), "notes": joiner.join(c[11]),
                   "ec": join_ec(c[12]), "source": "S10a"}
            rec["formula"], fb = join_formula(c[3])
            for b in fb:
                formula_breaks[b["kind"]] += 1
                if b["kind"] == "space" and COEF_RE.match(b["left"]):
                    numeric_breaks.append({"page": r["page"], "row": r["row"], "id": rec["id"],
                                           **b})
            # how many tokens a naive space join would have broken (for the report only)
            for tok in " ".join(split_frags(c[3])).split():
                if MET_TOKEN_RE.match(tok) is None and not FORMULA_COMPLETE.match(tok):
                    naive_bad_tokens += 1
                elif MET_TOKEN_RE.match(tok) and MET_TOKEN_RE.match(tok).group(1) not in known_mets:
                    naive_bad_tokens += 1
            rec["gpr_txt"], gb = join_gpr(c[5])
            for b in gb:
                gpr_breaks[b["kind"]] += 1
            if not rec["id"] or not any(a in rec["formula"].split() for a in ARROWS):
                unreadable.append(rec)
                continue
            rxns.append(rec)

        # ambiguous numeric breaks: would "<number><next token>" also be a valid metabolite?
        ambiguous_numeric = []
        for b in numeric_breaks:
            alt = MET_TOKEN_RE.match(b["left"] + b["right"])
            if alt and alt.group(1) in known_mets:
                ambiguous_numeric.append(b)

        # ------------------------------------------------------------------ S12 parse
        s12 = collections.OrderedDict()
        s12_exch = collections.OrderedDict()
        s12_failures = []
        for r in s12_raw:
            c = r["cells"]
            rid = join_nospace(c[0])
            ftxt, _ = join_formula(c[1])
            gtxt, _ = join_gpr(c[4])
            typ = joiner.join(c[5])
            entry = {"page": r["page"], "row": r["row"], "raw_id": rid, "formula": ftxt,
                     "lb_txt": join_nospace(c[2]), "ub_txt": join_nospace(c[3]),
                     "gpr_txt": gtxt, "type": typ}
            try:
                lhs, rhs, arrow = parse_formula(ftxt)
            except ParseError as exc:
                s12_failures.append({**entry, "problem": str(exc)})
                continue
            if "Internal exchange" in typ:
                m = re.match(r"^BT(.+)\[u\]tr$", rid)
                if not m or len(lhs) != 1 or len(rhs) != 1:
                    s12_failures.append({**entry, "problem": "unexpected exchange format"})
                    continue
                s12_exch[m.group(1)] = entry
                continue
            if not rid.startswith("BT"):
                s12_failures.append({**entry, "problem": "id without BT prefix"})
                continue

            def strip(terms):
                out = []
                for mid, comp, k in terms:
                    if not mid.startswith("BT"):
                        raise ParseError(f"metabolite {mid} without BT prefix")
                    out.append((mid[2:], comp, k))
                return out

            try:
                lhs, rhs = strip(lhs), strip(rhs)
            except ParseError as exc:
                s12_failures.append({**entry, "problem": str(exc)})
                continue
            entry.update({"lhs": lhs, "rhs": rhs, "arrow": arrow,
                          "tokens": [t[2:] if MET_TOKEN_RE.match(t) and t.startswith("BT") else t
                                     for t in ftxt.split()]})
            s12[rid[2:]] = entry

        # ------------------------------------------------------------------ unreadable S10a rows
        s10a_ids = {r["id"] for r in rxns}
        for u in unreadable:
            vis = u["formula"].split()
            cands = [k for k, e in s12.items()
                     if k not in s10a_ids and e["tokens"][:len(vis)] == vis]
            fix = {"table": "S10a", "page": u["page"], "row": u["row"],
                   "visible_text": {"BIGG Rxn ID": u["id"], "Formula": u["formula"],
                                    "Rev": u["rev"], "LB": u["lb_txt"], "UB": u["ub_txt"]},
                   "clipped_at_page_bottom": u["clipped"]}
            if len(cands) == 1 and u["clipped"]:
                e = s12[cands[0]]
                rec = dict(u)
                rec.update({"id": cands[0], "formula": " ".join(e["tokens"]),
                            "lb_txt": e["lb_txt"], "ub_txt": e["ub_txt"],
                            "gpr_txt": e["gpr_txt"], "rev": "0" if e["arrow"] == "=>" else "1",
                            "source": f"S12 p.{e['page']} (S10a p.{u['page']} row clipped)",
                            "name": u["name"] or cands[0]})
                rxns.append(rec)
                fix.update({"action": "row taken from Table S12",
                            "reaction": cands[0], "s12_page": e["page"],
                            "visible_tokens_matched": len(vis),
                            "s12_tokens": len(e["tokens"]),
                            "s12_lb": e["lb_txt"], "s12_ub": e["ub_txt"],
                            "s12_gene_rule": e["gpr_txt"],
                            "reason": ("The S10a cell is taller than the page and is clipped by "
                                       "the page edge (cell border at y=%s pt on a %s pt page); "
                                       "only the first %d formula tokens are printed and the "
                                       "BIGG id, Rev, LB, UB and GA cells are not visible. "
                                       "Table S12 is the only complete rendering in the PDF; its "
                                       "formula (BT prefix removed) starts with exactly the "
                                       "visible S10a tokens and its id is the only S12 "
                                       "reaction missing from S10a." % (
                                           next(i["rect_clip_pt"] + 595.3 for i in s10a_pages
                                                if i["page"] == u["page"]),
                                           "595.3", len(vis)))})
                manual.append({"kind": "clipped_row_replaced_from_S12", **fix})
            else:
                fix.update({"action": "none possible", "candidates": cands})
                parse_failures.append({"table": "S10a", "page": u["page"], "row": u["row"],
                                       "id": u["id"], "problem": "unreadable row", **fix})
        rxns.sort(key=lambda r: (r["page"], r["row"]))

        # rows recovered by the bottom edge that were otherwise complete
        recovered = []
        for i in s10a_pages:
            if i["rect_clip_pt"] > 0.5:
                recovered.append({"page": i["page"],
                                  "rows_with_bottom_edge": i["rows"],
                                  "rows_default_settings": i.get("rows_without_bottom_edge")})
        report["extraction"]["S10a"]["recovered_rows"] = recovered

        # ------------------------------------------------------------------ parse + validate S10a
        for rec in rxns:
            try:
                lhs, rhs, arrow = parse_formula(rec["formula"])
            except ParseError as exc:
                parse_failures.append({"table": "S10a", "page": rec["page"], "row": rec["row"],
                                       "id": rec["id"], "problem": f"formula: {exc}",
                                       "formula": rec["formula"]})
                rec["bad"] = True
                continue
            rec["lhs"], rec["rhs"], rec["arrow"] = lhs, rhs, arrow
            bad_mets = sorted({f"{m}[{c}]" for m, c, _ in lhs + rhs if m not in known_mets})
            bad_comp = sorted({f"{m}[{c}]" for m, c, _ in lhs + rhs if c not in ("c", "e")})
            if bad_mets or bad_comp:
                parse_failures.append({"table": "S10a", "page": rec["page"], "row": rec["row"],
                                       "id": rec["id"], "problem": "metabolite not in S10b"
                                       if bad_mets else "unknown compartment",
                                       "tokens": bad_mets + bad_comp, "formula": rec["formula"]})
                rec["bad"] = True
            both = {m for m, c, _ in lhs} & {m for m, c, _ in rhs}
            both = {(m, c) for m, c, _ in lhs if (m, c) in {(x, y) for x, y, _ in rhs}}
            if both:
                rec["same_met_both_sides"] = sorted(f"{m}[{c}]" for m, c in both)
            toks = gpr_tokens(rec["gpr_txt"])
            bad_g = [t for t in toks if t not in ("(", ")", "and", "or") and not GENE_RE.match(t)]
            syn = check_gpr(toks)
            if bad_g or syn:
                parse_failures.append({"table": "S10a", "page": rec["page"], "row": rec["row"],
                                       "id": rec["id"], "problem": "gene rule",
                                       "bad_tokens": bad_g, "syntax": syn, "GA": rec["gpr_txt"]})
                rec["bad"] = True
            rec["gpr_tokens"] = toks
            rec["gpr"] = normalize_gpr(toks)
            if not re.match(r"^[A-Za-z0-9_()\-]+$", rec["id"]):
                parse_failures.append({"table": "S10a", "page": rec["page"], "row": rec["row"],
                                       "id": rec["id"], "problem": "reaction id characters"})

        # bounds
        for rec in rxns:
            for key in ("lb", "ub"):
                val, flag = parse_number(rec[f"{key}_txt"])
                rec[key] = val
                rec[f"{key}_flag"] = flag
                if flag is not None:
                    e = s12.get(rec["id"])
                    s12_val = parse_number(e[f"{key}_txt"])[0] if e else None
                    if flag == "non_integer" and rec["source"] != "S10a":
                        continue  # value already taken from S12 for a clipped row
                    item = {"kind": "bound_cell_number_format", "table": "S10a",
                            "page": rec["page"], "row": rec["row"], "reaction": rec["id"],
                            "column": key.upper(), "printed": rec[f"{key}_txt"],
                            "readings": ({"decimal_point": float(rec[f"{key}_txt"]),
                                          "thousands_separator": float(
                                              rec[f"{key}_txt"].replace(".", ""))}
                                         if flag == "ambiguous_number_format" else None),
                            "s12_value": s12_val,
                            "s12_printed": e[f"{key}_txt"] if e else None,
                            "s12_page": e["page"] if e else None}
                    if s12_val is not None:
                        rec[key] = s12_val
                        item["action"] = f"used the Table S12 value {s12_val:g}"
                        item["reason"] = (
                            "Every other S10a bound cell is an integer (-1000, 0 or 1000); this "
                            "cell is printed with three decimals, so it is either a decimal "
                            "number or a number with a thousands separator, and the printing "
                            "alone does not say which. Table S12 lists the same reaction's "
                            "bound in an unambiguous format and was used.")
                        manual.append(item)
                    else:
                        parse_failures.append({**item, "problem": "unresolved bound"})

        # Rev vs bounds vs arrow
        rev_issues = []
        for rec in rxns:
            if rec.get("lb") is None or rec.get("ub") is None or "arrow" not in rec:
                continue
            bounds_rev = rec["lb"] < 0 < rec["ub"]
            arrow_rev = rec["arrow"] == "<=>"
            rev_flag = rec["rev"] == "1" if rec["rev"] in ("0", "1") else None
            if rev_flag is None or rev_flag != bounds_rev or arrow_rev != bounds_rev:
                rev_issues.append({"reaction": rec["id"], "page": rec["page"], "row": rec["row"],
                                   "Rev": rec["rev"], "arrow": rec["arrow"], "LB": rec["lb"],
                                   "UB": rec["ub"], "note": "bounds kept"})

        # ------------------------------------------------------------------ S12 cross-check
        cross = {"reactions_compared": 0, "stoichiometry_mismatch": [], "arrow_mismatch": [],
                 "bound_mismatch": [], "gpr_text_differs_logic_same": [],
                 "gpr_logic_mismatch": [], "missing_in_s12": [], "missing_in_s10a": [],
                 "s12_parse_failures": s12_failures, "exchange_metabolites": {}}
        by_id = {r["id"]: r for r in rxns}
        nonex = [r for r in rxns if not r["id"].startswith("EX_")]
        for rec in nonex:
            e = s12.get(rec["id"])
            if e is None:
                cross["missing_in_s12"].append(rec["id"])
                continue
            if rec["source"] != "S10a" or rec.get("bad"):
                continue
            cross["reactions_compared"] += 1
            st_a = stoichiometry(rec["lhs"], rec["rhs"])
            st_b = stoichiometry(e["lhs"], e["rhs"])
            keys = set(st_a) | set(st_b)
            if any(abs(st_a.get(k, 0) - st_b.get(k, 0)) > 1e-9 for k in keys):
                cross["stoichiometry_mismatch"].append(
                    {"reaction": rec["id"], "s10a_page": rec["page"], "s12_page": e["page"],
                     "s10a": rec["formula"], "s12": " ".join(e["tokens"])})
            if rec["arrow"] != e["arrow"]:
                cross["arrow_mismatch"].append({"reaction": rec["id"], "s10a": rec["arrow"],
                                                "s12": e["arrow"]})
            lb12, ub12 = parse_number(e["lb_txt"])[0], parse_number(e["ub_txt"])[0]
            if lb12 != rec["lb"] or ub12 != rec["ub"]:
                cross["bound_mismatch"].append({"reaction": rec["id"], "s10a": [rec["lb"], rec["ub"]],
                                                "s12": [lb12, ub12],
                                                "s10a_printed": [rec["lb_txt"], rec["ub_txt"]]})
            t12 = gpr_tokens(e["gpr_txt"])
            if normalize_gpr(t12) != rec["gpr"]:
                syn12 = check_gpr(t12)
                bad12 = [t for t in t12 if t not in ("(", ")", "and", "or") and not GENE_RE.match(t)]
                if syn12 is None and not bad12 and gpr_dnf(t12) == gpr_dnf(rec["gpr_tokens"]):
                    cross["gpr_text_differs_logic_same"].append(rec["id"])
                else:
                    if not t12:
                        problem = "S12 rule empty"
                    elif bad12:
                        problem = f"S12 rule has invalid token(s) {bad12}"
                    elif syn12:
                        problem = f"S12 rule malformed ({syn12})"
                    else:
                        problem = "logic differs"
                    cross["gpr_logic_mismatch"].append({
                        "reaction": rec["id"], "s10a": rec["gpr"], "s12": normalize_gpr(t12),
                        "s10a_rule_valid": check_gpr(rec["gpr_tokens"]) is None,
                        "s12_problem": problem, "s12_page": e["page"]})
        cross["missing_in_s10a"] = sorted(set(s12) - {r["id"] for r in nonex})
        ex_mets_a = {}
        for rec in rxns:
            if rec["id"].startswith("EX_") and "lhs" in rec:
                for m, c, k in rec["lhs"] + rec["rhs"]:
                    ex_mets_a[m] = rec["id"]
        cross["exchange_metabolites"] = {
            "s10a": len(ex_mets_a), "s12": len(s12_exch),
            "only_in_s10a": sorted(set(ex_mets_a) - set(s12_exch)),
            "only_in_s12": sorted(set(s12_exch) - set(ex_mets_a))}

        # ------------------------------------------------------------------ exchange orientation
        ex_issues = []
        for rec in rxns:
            if not rec["id"].startswith("EX_") or "lhs" not in rec:
                continue
            terms = rec["lhs"] + rec["rhs"]
            if len(terms) != 1 or rec["rhs"] or terms[0][1] != "e" or terms[0][2] != 1.0:
                ex_issues.append({"reaction": rec["id"], "formula": rec["formula"]})

        # ------------------------------------------------------------------ build the cobra model
        log("building the cobra model")
        model = cobra.Model("iAH991")
        model.name = ("iAH991 Bacteroides thetaiotaomicron VPI-5482 (Heinken et al. 2013), "
                      "rebuilt from supplementary Table S10")
        model.compartments = {"c": "cytosol", "e": "extracellular space"}
        used = collections.OrderedDict()
        for rec in rxns:
            for m, c, k in rec.get("lhs", []) + rec.get("rhs", []):
                used[(m, c)] = True
        met_objs = {}
        cid_seen = {}
        id_collisions = []
        for (m, c) in sorted(used, key=lambda x: (x[0], x[1])):
            cid = cobra_met_id(m, c)
            if cid in cid_seen:
                id_collisions.append({"cobra_id": cid, "published": [cid_seen[cid], f"{m}[{c}]"]})
            cid_seen[cid] = f"{m}[{c}]"
            info = mets_pub.get(m, {})
            met = cobra.Metabolite(cid, name=info.get("name", m), formula=info.get("formula") or None,
                                   charge=info.get("charge"), compartment=c)
            met.notes = {"published_id": esc(f"{m}[{c}]"),
                         "source": esc(f"Table S10b p.{info.get('page')} row {info.get('row')}")}
            ann = {}
            if re.fullmatch(r"cpd\d{5}", info.get("seed", "")):
                ann["seed.compound"] = info["seed"]
            if re.fullmatch(r"C\d{5}", info.get("kegg", "")):
                ann["kegg.compound"] = info["kegg"]
            elif re.fullmatch(r"G\d{5}", info.get("kegg", "")):
                ann["kegg.glycan"] = info["kegg"]
            met.annotation = ann
            met_objs[(m, c)] = met
        model.add_metabolites(list(met_objs.values()))

        rid_seen = {}
        cobra_rxns = []
        ex_by_met = {}
        for rec in rxns:
            cid = cobra_rxn_id(rec["id"])
            if cid in rid_seen:
                id_collisions.append({"cobra_id": cid, "published": [rid_seen[cid], rec["id"]]})
            rid_seen[cid] = rec["id"]
            rec["cobra_id"] = cid
            rxn = cobra.Reaction(cid, name=rec["name"] or rec["id"], subsystem=rec["subsystem"],
                                 lower_bound=rec["lb"], upper_bound=rec["ub"])
            st = stoichiometry(rec["lhs"], rec["rhs"])
            rxn.add_metabolites({met_objs[k]: v for k, v in st.items() if v != 0})
            rxn.gene_reaction_rule = rec["gpr"]
            notes = {"published_id": rec["id"], "published_formula": rec["formula"],
                     "published_Rev": rec["rev"], "source": rec["source"] if rec["source"] != "S10a"
                     else f"Table S10a p.{rec['page']} row {rec['row']}"}
            for k, label in (("seed", "SEED_id"), ("cs", "confidence_score"), ("ec", "EC_number"),
                             ("references", "references"), ("notes", "curation_notes")):
                if rec.get(k):
                    notes[label] = rec[k]
            rxn.notes = {k: esc(v) for k, v in notes.items()}
            ann = {}
            if re.fullmatch(r"rxn\d{5}", rec["seed"]):
                ann["seed.reaction"] = rec["seed"]
            ecs = [x for x in re.split(r"[,\s]+|\bor\b", rec["ec"]) if EC_RE.match(x)]
            if ecs:
                ann["ec-code"] = ecs
            rxn.annotation = ann
            cobra_rxns.append(rxn)
            if rec["id"].startswith("EX_"):
                ex_by_met[rec["lhs"][0][0]] = cid
        model.add_reactions(cobra_rxns)
        # cobra collects genes from the rules through sets, whose order depends on string hashing
        # (randomised per process); sort them so that the SBML file is byte-identical across runs
        model.genes.sort(key=lambda g: g.id)

        # subsystems as SBML groups (cobra keeps reaction.subsystem only via groups)
        groups = []
        for sub in sorted({r.subsystem for r in model.reactions if r.subsystem}):
            gid = "subsystem_" + re.sub(r"[^A-Za-z0-9_]", "_", sub)
            g = OrderedGroup(gid, name=sub, kind="partonomy")
            g.add_members([r for r in model.reactions if r.subsystem == sub])
            groups.append(g)
        model.add_groups(groups)

        biomass = [r for r in rxns if "biomass" in r["id"].lower() or "biomass" in r["name"].lower()]
        if len(biomass) != 1:
            raise RuntimeError(f"expected one biomass reaction, found {[b['id'] for b in biomass]}")
        objective_id = biomass[0]["cobra_id"]
        model.objective = objective_id
        model.notes = {"description": esc(
            "Rebuilt from Table S10a/S10b of the supplementary PDF (sha256 " + pdf_sha + ") by "
            "scripts/rebuild_iAH991_from_pdf.py; see REBUILD.md for every manual intervention.")}

        # S12-only search for other biomass reactions (in the whole S12 B. theta part)
        s12_biomass = [k for k in s12 if "biomass" in k.lower()]

        out_xml = f"{args.outdir}/iAH991_rebuilt.xml"
        log(f"writing {out_xml}")
        cobra.io.write_sbml_model(model, out_xml)

        # ------------------------------------------------------------------ round trip
        m2 = cobra.io.read_sbml_model(out_xml)
        rt = {"reactions": len(m2.reactions) == len(model.reactions),
              "metabolites": len(m2.metabolites) == len(model.metabolites),
              "genes": len(m2.genes) == len(model.genes), "differences": []}
        for r in model.reactions:
            r2 = m2.reactions.get_by_id(r.id)
            st1 = {m.id: v for m, v in r.metabolites.items()}
            st2 = {m.id: v for m, v in r2.metabolites.items()}
            if (st1 != st2 or r.bounds != r2.bounds
                    or r.gene_reaction_rule.replace("(", "").replace(")", "") !=
                    r2.gene_reaction_rule.replace("(", "").replace(")", "")
                    or r.subsystem != r2.subsystem):
                rt["differences"].append(r.id)
        rt["objective"] = [r.id for r in m2.reactions if r.objective_coefficient != 0]
        rt["ok"] = (rt["reactions"] and rt["metabolites"] and rt["genes"]
                    and not rt["differences"] and rt["objective"] == [objective_id])
        _, sbml_errors = cobra.io.validate_sbml_model(out_xml)
        rt["libsbml_validation"] = {k: len(v) for k, v in sbml_errors.items()}
        rt["libsbml_messages_sample"] = {k: [str(x)[:300] for x in v[:3]]
                                         for k, v in sbml_errors.items() if v}

        # ------------------------------------------------------------------ counts
        n_ex = sum(r.id.startswith("EX_") for r in m2.reactions)
        n_dm = sum(r.id.startswith("DM_") for r in m2.reactions)
        n_sink = sum(r.id.startswith("sink_") for r in m2.reactions)
        counts = {"reactions": len(m2.reactions), "exchange": n_ex, "demand": n_dm, "sink": n_sink,
                  "exchange_and_demand": n_ex + n_dm + n_sink,
                  "metabolic_and_transport": len(m2.reactions) - n_ex - n_dm - n_sink,
                  "metabolites": len(m2.metabolites),
                  "unique_metabolites_without_compartment": len({m for m, c in used}),
                  "s10b_metabolites": len(mets_pub),
                  "s10b_metabolites_unused": sorted(set(mets_pub) - {m for m, c in used}),
                  "genes": len(m2.genes)}

        # ------------------------------------------------------------------ validation set-up
        m2.solver = SOLVER
        mets_info = {k: v for k, v in mets_pub.items()}
        sim = Sim(m2, ex_by_met, mets_info)
        s8a = {e["id"]: e["rate"] for e in s8["S8a"]["entries"]}
        s8b = {e["id"]: e["rate"] for e in s8["S8b"]["entries"]}
        s8c = {e["id"]: e["rate"] for e in s8["S8c"]["entries"]}
        s8d = {e["id"]: e["rate"] for e in s8["S8d"]["entries"]}

        def rate_for(mid):
            if mid in s8b:
                return s8b[mid], "S8b"
            r_ = hexose_equivalent_rate(mets_pub.get(mid, {}).get("formula"))
            return r_, "60 mmol C/gDW/h (S3 scaling)"

        def condition(subs):
            up = dict(s8a)
            used_rates = {}
            for mid, rate in subs.items():
                if rate is None:
                    rate, how = rate_for(mid)
                else:
                    how = "fixed"
                up[mid] = rate
                used_rates[mid] = {"rate": rate, "source": how}
            return up, used_rates

        validation = collections.OrderedDict()

        # S3 transcription check against the PDF column
        col_vals = [t for (_, _, t) in s3_col if re.fullmatch(r"\d+(?:\.\d+)?", t)]
        s3_check = col_vals == [row[2] for row in S3_ROWS]

        log("validation: Table S3 growth rates")
        s3_results = []
        for label, subs, pub in S3_ROWS + S3_ALTERNATIVES:
            up, used_rates = condition(subs)
            g, missing = sim.growth(up)
            s3_results.append({"condition": label, "uptake": used_rates, "published": pub,
                               "rebuilt": num(g), "abs_diff": num(abs(g - float(pub))),
                               "rounds_to_published": rounds_to(g, pub),
                               "missing_exchange": missing,
                               "alternative_mapping": label.startswith("Amylopectin (alt")})
        main_rows = [r for r in s3_results if not r["alternative_mapping"]]
        validation["table_S3_growth"] = {
            "transcription_matches_pdf_column": s3_check,
            "pdf_column_values": col_vals,
            "setup": "Table S8a minerals/cofactors + Table S8b carbon source; every other "
                     "exchange lower bound set to 0; upper bounds as published; sinks/demands "
                     "as published; objective Biomass_BT_v2",
            "rows": s3_results,
            "max_abs_diff": max(r["abs_diff"] for r in main_rows),
            "n_rounding_to_published": sum(r["rounds_to_published"] for r in main_rows),
            "n": len(main_rows)}

        # Diagnostics for S3 rows that do not reproduce. They do not change the comparison above.
        # (1) Is each printed S8b rate consistent with the scaling used for its peers?
        def n_carbon(mid):
            f = mets_pub.get(mid, {}).get("formula") or ""
            return sum(int(x) if x else 1 for x in re.findall(r"C(\d*)(?![a-z])", f))

        s8b_check = []
        for e in s8["S8b"]["entries"]:
            c_ = n_carbon(e["id"])
            r60 = 60.0 / c_ if c_ else None
            s8b_check.append({"id": e["id"], "printed_rate": e["rate"], "n_carbon_S10b": c_,
                              "rate_60_mmol_C": round(r60, 5) if r60 else None,
                              "ratio_printed_to_60C": round(e["rate"] / r60, 3) if r60 else None,
                              "deviates_gt_2pct": bool(r60 and abs(e["rate"] / r60 - 1) > 0.02)})
        s3_by_sub = {next(iter(subs)): (label, pub) for label, subs, pub in S3_ROWS if len(subs) == 1}
        deviating = []
        for chk in s8b_check:
            if not chk["deviates_gt_2pct"] or chk["id"] not in s3_by_sub:
                continue
            label, pub = s3_by_sub[chk["id"]]
            g_print, _ = sim.growth(condition({chk["id"]: chk["printed_rate"]})[0])
            g_60, _ = sim.growth(condition({chk["id"]: chk["rate_60_mmol_C"]})[0])
            deviating.append({"id": chk["id"], "s3_row": label, "published": pub,
                              "printed_rate": chk["printed_rate"],
                              "growth_at_printed_rate": num(g_print),
                              "rate_60_mmol_C": chk["rate_60_mmol_C"],
                              "growth_at_60_mmol_C": num(g_60),
                              "printed_rate_reproduces": rounds_to(g_print, pub),
                              "60C_rate_reproduces": rounds_to(g_60, pub)})
        # (2) alternative substrate mapping for the Amylopectin row
        diag_rows = []
        for row in main_rows:
            if row["rounds_to_published"]:
                continue
            alts = [{"what": d["id"] + " at " + str(d["rate_60_mmol_C"]) + " (60 mmol C/gDW/h) "
                             "instead of the printed " + str(d["printed_rate"]),
                     "rebuilt": d["growth_at_60_mmol_C"],
                     "rounds_to_published": d["60C_rate_reproduces"]}
                    for d in deviating if d["s3_row"] == row["condition"]]
            alts += [{"what": a["condition"], "rebuilt": a["rebuilt"],
                      "rounds_to_published": a["rounds_to_published"]}
                     for a in s3_results if a["alternative_mapping"]
                     and a["condition"].split(" (")[0] == row["condition"]]
            diag_rows.append({"condition": row["condition"], "published": row["published"],
                              "rebuilt": row["rebuilt"], "alternatives": alts})
        all_alt = {d["condition"]: min((a for a in d["alternatives"] if a["rounds_to_published"]),
                                       key=lambda a: abs(a["rebuilt"] - float(d["published"])),
                                       default=None) for d in diag_rows}
        diffs_alt = []
        for row in main_rows:
            a = all_alt.get(row["condition"])
            g = a["rebuilt"] if a else row["rebuilt"]
            diffs_alt.append(abs(g - float(row["published"])))
        validation["table_S3_growth"]["diagnostics"] = {
            "s8b_rate_vs_60_mmol_C": s8b_check,
            "s8b_entries_deviating_from_60C_with_S3_row": deviating,
            "non_reproducing_rows": diag_rows,
            "max_abs_diff_if_alternatives_accepted": num(max(diffs_alt)),
            "n_rows_reproduced_if_alternatives_accepted": sum(
                1 for row in main_rows if row["rounds_to_published"] or all_alt.get(row["condition"]))}

        # ATPM: growth on glucose minimal medium (S8c) under the possible readings of its LB
        atpm_diag = {}
        atpm_rxn = m2.reactions.get_by_id("ATPM")
        for val in (0.0, 8.43, 84.3):
            with m2:
                sim.apply(s8c)
                atpm_rxn.lower_bound = val
                g = m2.slim_optimize(error_value=float("nan"))
            atpm_diag[f"LB={val:g}"] = "infeasible" if math.isnan(g) else num(g)
        atpm_diag["LB=84300"] = "invalid (lower bound above the upper bound of 1000)"
        validation["atpm_lower_bound_readings_glucose_S8c"] = atpm_diag

        # sensitivity: sinks closed (findExcRxns in the COBRA Toolbox would also close them)
        sens = []
        sinks = [r for r in m2.reactions if r.id.startswith("sink_")]
        with m2:
            for r in sinks:
                r.lower_bound = 0.0
            for label, subs, pub in S3_ROWS:
                up, _ = condition(subs)
                g, _ = sim.growth(up)
                sens.append({"condition": label, "published": pub, "rebuilt_sinks_closed": num(g)})
        sink_role = {}
        with m2:
            sim.apply(s8c)
            sol = m2.optimize()
            fluxes = {r.id: num(sol.fluxes[r.id]) for r in sinks}
        for r in sinks:
            with m2:
                sim.apply(s8c)
                r.lower_bound = 0.0
                g_closed = m2.slim_optimize(error_value=float("nan"))
            met = next(iter(r.metabolites))
            sink_role[r.id] = {"flux_at_optimum_glucose": fluxes[r.id],
                               "growth_glucose_if_uptake_closed": 0.0 if math.isnan(g_closed)
                               else num(max(g_closed, 0.0)),
                               "other_reactions_of_metabolite": sorted(
                                   x.id for x in met.reactions if x.id != r.id)}
        validation["table_S3_growth"]["sensitivity_sinks_closed"] = {
            "max_abs_diff": max(abs(s["rebuilt_sinks_closed"] - float(s["published"])) for s in sens),
            "n_published_growth_lost": sum(1 for s in sens if s["rebuilt_sinks_closed"] < GROWTH_EPS
                                           and float(s["published"]) > 0),
            "sink_roles": sink_role,
            "rows": sens}

        # S3b secretion (FVA at 100 % of optimum)
        log("validation: Table S3b secretion ranges (FVA)")
        s3b = []
        for label, mid, pub in S3B_ROWS:
            up = dict(s8a)
            up[mid] = 10.0
            with m2:
                sim.apply(up)
                ids = [ex_by_met[x] for x in ("ac", "ppa", "succ", "co2", "h2")]
                fva = flux_variability_analysis(m2, reaction_list=ids, fraction_of_optimum=1.0,
                                                processes=1)
            row = {"substrate": label, "published": pub, "rebuilt": {}, "matches": {}}
            for x in ("ac", "ppa", "succ", "co2", "h2"):
                lo, hi = fva.loc[ex_by_met[x], "minimum"], fva.loc[ex_by_met[x], "maximum"]
                row["rebuilt"][x] = [num(lo, 4), num(hi, 4)]
                p_ = pub[x]
                if p_.startswith("<"):
                    ok = abs(lo) < float(p_[1:]) and abs(hi) < float(p_[1:])
                else:
                    parts = [float(t) for t in p_.split(" to ")]
                    plo, phi = (parts[0], parts[-1])
                    ok = rounds_to(float(lo), f"{plo:.2f}") and rounds_to(float(hi), f"{phi:.2f}")
                row["matches"][x] = ok
            s3b.append(row)
        validation["table_S3b_secretion"] = {
            "method_assumed": "FVA of the exchange fluxes at 100% of maximal growth, Table S8a + "
                              "10 mmol/gDW/h substrate (the paper does not state the method)",
            "rows": s3b,
            "n_values_matching": sum(v_ for r_ in s3b for v_ in r_["matches"].values()),
            "n_values": sum(len(r_["matches"]) for r_ in s3b)}

        # S1
        log("validation: Table S1 carbon sources")
        s1_alias = {"ara-L": "arab-L"}
        s1_results = []
        for row in s1_rows:
            members = []
            for mid in dict.fromkeys(row["ids"]):
                tested = s1_alias.get(mid, mid)
                up, used_rates = condition({tested: None})
                g, missing = sim.growth(up)
                members.append({"id": mid, "tested_as": tested, "rate": used_rates[tested],
                                "growth": num(g), "grows": g > GROWTH_EPS,
                                "has_exchange": not missing})
            grows_any = any(m["grows"] for m in members)
            s1_results.append({"carbohydrate": row["name"], "ids": row["ids"],
                               "published_iAH991": row["iAH991"],
                               "published_BTheta_Seed_v1": row["seed_v1"],
                               "members": members, "rebuilt": "yes" if grows_any else "no",
                               "agrees": ("yes" if grows_any else "no") == row["iAH991"]})
        extra_core8 = None
        if any("core7" in r["ids"] and r["ids"].count("core7") > 1 for r in s1_rows):
            up, used_rates = condition({"core8": None})
            g, _ = sim.growth(up)
            extra_core8 = {"id": "core8", "rate": used_rates["core8"], "growth": num(g)}
        all_members = [m for r in s1_results for m in r["members"]]
        s1_diag = []
        for mem in all_members:
            if mem["grows"]:
                continue
            row_ = {"id": mem["id"], "growth_at_set_up_rate": mem["growth"]}
            for rate in (10.0, 100.0):
                g, _ = sim.growth(condition({mem["tested_as"]: rate})[0])
                row_[f"growth_at_{rate:g}"] = num(g)
            ex_rxn = m2.reactions.get_by_id(ex_by_met[mem["tested_as"]]) \
                if mem["tested_as"] in ex_by_met else None
            if ex_rxn is not None:
                met_e = next(iter(ex_rxn.metabolites))
                row_["reactions_using_it"] = sorted(r.id for r in met_e.reactions if r.id != ex_rxn.id)
            s1_diag.append(row_)
        validation["table_S1_carbon_sources"] = {
            "setup": "Table S8a + the listed carbohydrate as the only carbon source, at the "
                     "Table S8b rate when listed there, otherwise at 60 mmol C/gDW/h "
                     "(10 mmol hexose units, the Table S3 scaling); grows = growth > 1e-6/h",
            "rows": s1_results,
            "rows_agreeing": sum(r["agrees"] for r in s1_results), "n_rows": len(s1_results),
            "ids_growing": sum(m["grows"] for m in all_members), "n_ids": len(all_members),
            "ids_not_growing": [m["id"] for m in all_members if not m["grows"]],
            "diagnostics_non_growing_ids": s1_diag,
            "id_aliases_used": s1_alias, "core8_extra": extra_core8}

        # gene essentiality
        log("validation: gene essentiality")
        ess_text = re.sub(r"\s+", " ", text_p8)
        pub_ess = {}
        mm = re.search(r"In glucose minimal medium, (\d+)/(\d+)", ess_text)
        if mm:
            pub_ess["glucose_minimal_S8c"] = int(mm.group(1))
        mm = re.search(r"In rich medium, only (\d+)/(\d+)", ess_text)
        if mm:
            pub_ess["rich"] = int(mm.group(1))
        mm = re.search(r"(\d+) of the (\d+) in silico essential and (\d+) ?of the (\d+) "
                       r"non-essential", ess_text)
        if mm:
            pub_ess["TYG_S8d"] = int(mm.group(2))
            pub_ess["TYG_S8d_nonessential"] = int(mm.group(4))
        s9 = parse_s9_western(pdf, known_mets, {k for k in ex_by_met})

    # PDF closed from here on --------------------------------------------------------------

    s9_up = {}
    s9_unmapped = []
    ex_id_by_name = {}
    for rec in rxns:
        if rec["id"].startswith("EX_"):
            m_ = re.match(r"^EX_(.+)\(e\)$", rec["id"])
            if m_:
                ex_id_by_name[m_.group(1)] = rec["lhs"][0][0]
    for e in s9["entries"]:
        mid = e["id"]
        if mid in ex_by_met:
            s9_up[mid] = e["western"]
        elif mid.replace("-", "_") in ex_id_by_name:
            s9_up[ex_id_by_name[mid.replace("-", "_")]] = e["western"]
        else:
            s9_unmapped.append(mid)

    rich = {m: 1000.0 for m in ex_by_met}
    media = collections.OrderedDict([
        ("glucose_minimal_S8c", s8c), ("TYG_S8d", s8d), ("rich", rich),
        ("western_diet_S9", s9_up)])
    gene_ids = sorted(g.id for g in m2.genes)

    def knockout_growth(up):
        with m2:
            missing_ = sim.apply(up)
            wt_ = m2.slim_optimize(error_value=float("nan"))
            df = single_gene_deletion(m2, gene_list=gene_ids, processes=args.processes)
        ko_ = {}
        for ids, gr in zip(df["ids"], df["growth"]):
            ko_[next(iter(ids))] = 0.0 if (gr is None or math.isnan(gr)) else float(gr)
        return (0.0 if math.isnan(wt_) else wt_), ko_, missing_

    ess_results = collections.OrderedDict()
    for name, up in media.items():
        t0 = time.time()
        wt, ko, missing = knockout_growth(up)
        ess_abs = sorted(g for g, v in ko.items() if v < GROWTH_EPS)
        res = {"wild_type_growth": num(wt), "missing_exchange_for": missing,
               "essential_growth_below_1e-6": len(ess_abs), "essential_genes": ess_abs}
        for frac in (0.01, 0.05, 0.1, 0.5):
            res[f"essential_below_{int(frac * 100)}pct_of_wt"] = sum(
                v < frac * wt for v in ko.values()) if wt > 0 else None
        res["published"] = pub_ess.get(name) if name != "western_diet_S9" else 160
        ess_results[name] = res
        log(f"  {name}: WT {wt:.4f}, essential {len(ess_abs)} (published {res['published']}), "
            f"{time.time() - t0:.0f} s")

    # Cross-reference with Table S14: genes whose deletion the paper reports as viable in silico
    # on glucose minimal medium, but which are essential in the rebuild.
    s14_glc_viable = sorted({g for label, gsets, cond, subs, pub in S14_ROWS
                             if cond == "Glucose minimal medium" and pub == "yes"
                             for gs in gsets for g in gs})
    glc_ess = set(ess_results["glucose_minimal_S8c"]["essential_genes"])
    conflict = sorted(glc_ess & set(s14_glc_viable))
    ess_results["glucose_minimal_S8c"]["s14_reports_viable_but_essential_here"] = conflict
    ess_results["glucose_minimal_S8c"]["count_without_those"] = len(glc_ess - set(conflict))

    # Exploratory: TYG count when vitamins that iAH991 can take up but S8d does not list are added.
    probes = [v_ for v_ in TYG_VITAMIN_PROBES if v_ in ex_by_met and v_ not in s8d]
    tyg_base = set(ess_results["TYG_S8d"]["essential_genes"])
    tyg_probe = []
    for k in range(1, len(probes) + 1):
        for combo in itertools.combinations(probes, k):
            up = dict(s8d)
            up.update({v_: 10.0 for v_ in combo})
            wt_, ko_, _ = knockout_growth(up)
            ess_ = {g for g, v_ in ko_.items() if v_ < GROWTH_EPS}
            tyg_probe.append({"added_at_10": list(combo), "wild_type_growth": num(wt_),
                              "essential": len(ess_),
                              "no_longer_essential": sorted(tyg_base - ess_)})
    ess_results["TYG_S8d"]["exploratory_vitamin_additions"] = tyg_probe
    validation["gene_essentiality"] = {
        "definition": "single-gene deletion (cobra single_gene_deletion, FBA); a gene is "
                      "essential if growth after deletion is < 1e-6/h; counts at relative "
                      "thresholds are given for sensitivity",
        "published_counts_source": f"supplementary text p.{TEXT_PAGE_ESSENTIALITY} (glucose "
                                   "minimal 204/991, rich 61/991, TYG 116 essential/875 "
                                   "non-essential); Western diet 160 from the main text",
        "media": ess_results,
        "western_diet_note": "Table S9 defines the Western diet for the joint mouse/B. theta "
                             "model (lumen uptake rates); applying it to iAH991 alone is an "
                             "approximation of the paper's setting.",
        "western_diet_unmapped_ids": s9_unmapped}

    # S14
    log("validation: Table S14 knockouts")
    s14_pub_col = [t for (_, _, t) in s14_col if t in ("yes", "no")]
    s14_results = []
    for label, gene_sets, cond, substrate_sets, pub in S14_ROWS:
        cases = []
        for gs in gene_sets:
            for subs in substrate_sets:
                up, used_rates = condition(subs)
                g, missing = sim.growth(up, knockouts=gs)
                cases.append({"deleted": [x for x in gs if x in m2.genes],
                              "deleted_not_in_model": [x for x in gs if x not in m2.genes],
                              "substrates": used_rates, "growth": num(g),
                              "grows": g > GROWTH_EPS, "missing_exchange": missing})
        testable = [c_ for c_ in cases if not c_["missing_exchange"]]
        calls = {("yes" if c_["grows"] else "no") for c_ in testable}
        rebuilt = calls.pop() if len(calls) == 1 else ("mixed" if calls else "not testable")
        explanation = []
        if rebuilt != pub:
            for gs in gene_sets:
                for subs in substrate_sets:
                    up, _ = condition(subs)
                    with m2:
                        sim.apply(up)
                        for g in gs:
                            if g in m2.genes:
                                m2.genes.get_by_id(g).knock_out()
                        sol = m2.optimize()
                        grows = sol.status == "optimal" and sol.objective_value > GROWTH_EPS
                        still_active = []
                        if grows:
                            for g in gs:
                                for r in m2.genes.get_by_id(g).reactions:
                                    if r.bounds != (0.0, 0.0) and abs(sol.fluxes[r.id]) > 1e-9:
                                        still_active.append({"reaction": r.id,
                                                             "flux": num(sol.fluxes[r.id]),
                                                             "gene_rule": r.gene_reaction_rule})
                    lethal = []
                    if not grows:
                        for g in gs:
                            if g not in m2.genes:
                                continue
                            g_alone, _ = sim.growth(up, knockouts=[g])
                            if g_alone > GROWTH_EPS:
                                continue
                            lethal_rxns = []
                            for r in m2.genes.get_by_id(g).reactions:
                                with m2:
                                    sim.apply(up)
                                    m2.genes.get_by_id(g).knock_out()
                                    disabled = r.bounds == (0.0, 0.0)
                                if not disabled:
                                    continue
                                with m2:
                                    sim.apply(up)
                                    r.knock_out()
                                    gr = m2.slim_optimize(error_value=float("nan"))
                                if math.isnan(gr) or gr < GROWTH_EPS:
                                    lethal_rxns.append({"reaction": r.id, "equation": r.reaction,
                                                        "gene_rule": r.gene_reaction_rule})
                            lethal.append({"gene": g, "lethal_reactions": lethal_rxns})
                    if grows and pub == "no":
                        explanation.append({"deleted": gs, "substrates": list(subs),
                                            "reactions_of_deleted_genes_still_carrying_flux":
                                                still_active})
                    elif not grows and pub == "yes":
                        explanation.append({"deleted": gs, "substrates": list(subs),
                                            "single_deletions_that_are_lethal": lethal})
        s14_results.append({"deleted": label, "condition": cond, "published_in_silico": pub,
                            "rebuilt": rebuilt, "agrees": rebuilt == pub, "cases": cases,
                            "explanation": explanation})
    validation["table_S14_knockouts"] = {
        "transcription_matches_pdf_column": s14_pub_col == [r[4] for r in S14_ROWS],
        "rows": s14_results,
        "n_agree": sum(r["agrees"] for r in s14_results), "n": len(s14_results)}

    # ------------------------------------------------------------------ assemble report
    report["counts"] = {"rebuilt": counts, "published": PUBLISHED_COUNTS,
                        "match": {k: counts.get(k) == v for k, v in PUBLISHED_COUNTS.items()}}
    report["joins"] = {
        "formula_line_breaks": dict(formula_breaks),
        "formula_breaks_after_a_number": len(numeric_breaks),
        "formula_breaks_after_a_number_also_valid_as_metabolite": ambiguous_numeric,
        "gene_rule_line_breaks": dict(gpr_breaks),
        "text_columns": dict(joiner.stats),
        "tokens_a_naive_space_join_would_break": naive_bad_tokens,
    }
    report["parse_failures"] = parse_failures
    report["manual_interventions"] = manual
    report["rev_vs_bounds_disagreements"] = rev_issues
    report["exchange_orientation_issues"] = ex_issues
    report["same_metabolite_on_both_sides"] = [
        {"reaction": r["id"], "metabolites": r["same_met_both_sides"]} for r in rxns
        if r.get("same_met_both_sides")]
    report["id_collisions"] = id_collisions
    report["renamed_reaction_ids"] = {r["id"]: r["cobra_id"] for r in rxns
                                      if r["id"] != r["cobra_id"] and not r["id"].startswith("EX_")}
    report["s12_crosscheck"] = cross
    report["biomass"] = {
        "s10a_biomass_reactions": [b["id"] for b in biomass],
        "s12_bt_biomass_reactions": s12_biomass,
        "objective": objective_id,
        "evidence_checked": {
            "s10a_rows_named_biomass": len(biomass),
            "s10c_p307_mentions_bio01632_and_Biomass_BT_v2":
                "bio01632" in s10c_text and "Biomass_BT_v2" in s10c_text,
            "s12_bt_biomass_rows": len(s12_biomass)},
        "evidence": [
            "Table S10a contains a single biomass row (p. 89; clipped, see manual interventions).",
            "Table S10c (p. 307) lists the Model SEED draft biomass bio01632 as rejected and "
            "'Replaced with Biomass_BT_v2' ('The biomass reaction was modified based on "
            "literature').",
            "Table S12 (joint model) contains exactly one B. thetaiotaomicron biomass reaction, "
            "BTBiomass_BT_v2 (p. 475).",
        ],
        "formula": next(r["formula"] for r in rxns if r["cobra_id"] == objective_id)}
    report["sbml_round_trip"] = rt
    report["validation"] = validation
    report["s8_tables"] = {k: v for k, v in s8.items()}
    report["s9_western"] = s9

    with open(f"{args.outdir}/rebuild_report.json", "w") as fh:
        json.dump(report, fh, indent=1, default=str)
    write_markdown(report, f"{args.outdir}/REBUILD.md")
    log(f"done in {time.time() - t_start:.0f} s")


# --------------------------------------------------------------------------------------
# REBUILD.md
# --------------------------------------------------------------------------------------
def write_markdown(rep, path):
    c = rep["counts"]
    v = rep["validation"]
    rb = c["rebuilt"]
    s3 = v["table_S3_growth"]
    s3d = s3["diagnostics"]
    s3b = v["table_S3b_secretion"]
    s1 = v["table_S1_carbon_sources"]
    ge = v["gene_essentiality"]
    s14 = v["table_S14_knockouts"]
    cx = rep["s12_crosscheck"]
    j = rep["joins"]
    L = []
    L.append("# iAH991 rebuilt from the supplementary PDF\n")
    L.append("Model: iAH991, *Bacteroides thetaiotaomicron* VPI-5482. "
             f"Source: {rep['source']['citation']}.\n")
    L.append(f"Source PDF: `source/{rep['source']['pdf_file']}` "
             f"({rep['source']['pages']} pages), sha256 `{rep['source']['sha256']}`.\n")
    L.append("Files: `iAH991_rebuilt.xml` (SBML L3 FBC v2, written by cobra), "
             "`rebuild_report.json` (every number in this file, plus per-row details), this file.\n")
    L.append("Rebuilt by `scripts/rebuild_iAH991_from_pdf.py` "
             f"(python {rep['tool']['python']}, cobra {rep['tool']['cobra']}, "
             f"pdfplumber {rep['tool']['pdfplumber']}, solver {rep['tool']['solver']}). "
             "Re-run from the repository root with\n"
             "`python3 -I scripts/rebuild_iAH991_from_pdf.py "
             "models/curated/iAH991/source/2012GUTMICROBES0043R-Sup.pdf models/curated/iAH991`. "
             "This file is generated by that script.\n")

    # ---------------------------------------------------------------- method
    L.append("## 1. Method\n")
    L.append(
        "1. **Extraction.** Every page of Table S10a (pp. 36-262) and Table S10b (pp. 263-306) is "
        "read with pdfplumber `extract_tables`; the header row is checked on every page, and the "
        "pages just outside each range are checked not to carry it. An explicit horizontal edge "
        "0.25 pt above the page bottom is added so that rows whose cell borders run past the page "
        "edge are still returned (pdfplumber silently drops them otherwise; this affects pp. 89 "
        "and 179 only, and no other page changes).\n"
        "2. **Wrapped cells.** Formulas: a line break is a space when the fragment before it ends "
        "with a complete token (a metabolite with its `[compartment]`, a coefficient such as `2`, "
        "`0.5` or `(0.5)`, `+`, `=>` or `<=>`); otherwise it falls inside a token (`ala-` / "
        "`L[c]`, `10m3hddcaACP[` / `c]`) and the fragments are joined without a space. Gene "
        "rules: the same with tokens `BT_dddd`, `and`, `or`, parentheses. Reaction ids, SEED ids "
        "and numbers never contain spaces and are joined without one (`2FUCLAC_FUCA` / `SEe`). "
        "Free text (names, subsystems, notes) is joined with a space, except after a trailing "
        "hyphen or before closing punctuation, and a word broken without a hyphen "
        "(`Glycerophosph` / `olipid`) is re-joined when the whole word occurs elsewhere in the "
        "PDF and one of the two pieces never occurs as a complete word.\n"
        "3. **Checks.** Every formula token is checked against the S10b metabolite list and every "
        "gene token against `^BT_\\d{4}$` (plus a syntax check of the Boolean rule). For every "
        "line break that follows a number, the alternative reading (the number glued to the next "
        "token as one metabolite id) is also checked against S10b.\n"
        "4. **Independent cross-check.** Table S12 (the joint mouse/*B. theta* model, pp. 472-490) "
        "prints every iAH991 reaction a second time with a `BT` prefix and with different column "
        "widths, so its line breaks fall elsewhere. Every non-exchange reaction's stoichiometry, "
        "direction, bounds and gene rule is compared with the S10a parse; exchange reactions are "
        "compared by the exchanged metabolite (their S12 bounds describe the lumen, not iAH991).\n"
        "5. **Identifiers.** Metabolites `ala-L[c]` -> `ala__L_c` (`-` -> `__`, `[x]` -> `_x`); "
        "reactions: published BiGG id with `(e)` -> `_e`, `-` -> `__`, other illegal characters "
        "-> `_` (only `PFK(ppi)` -> `PFK_ppi_` needed it). Published ids and formulas, Rev, SEED "
        "ids, confidence scores, EC numbers, references, curation notes and the page/row of "
        "origin are kept in the SBML notes (XML-escaped); SEED/KEGG ids and EC numbers are also "
        "written as annotations. Metabolite names, formulas and charges come from S10b. "
        "Subsystems are written as SBML groups.\n"
        "6. **Bounds and gene rules** are as published in S10a; Rev is only compared with them. "
        "Every S10a exchange reaction is written `met[e] <=>`, i.e. it consumes its metabolite "
        "(coefficient -1), so none needed reorienting.\n"
        "7. **Validation set-up** (section 5): every exchange lower bound is set to 0, then the "
        "medium's uptakes are opened (lower bound = -rate) exactly as listed in Table S8; upper "
        "bounds, the 3 demand and 2 sink reactions stay as published; objective `Biomass_BT_v2`; "
        "GLPK. 'Grows' means growth > 1e-6 /h (the smallest published non-zero rate is 0.0012).\n")

    # ---------------------------------------------------------------- manual interventions
    L.append("## 2. Manual interventions\n")
    rows = []
    for m in rep["manual_interventions"]:
        if m["kind"] == "clipped_row_replaced_from_S12":
            rows.append([f"S10a p.{m['page']} row {m['row']}", m["reaction"],
                         f"The cell is taller than the page and is clipped by the page edge: only "
                         f"the first {m['visible_tokens_matched']} formula tokens are printed, and "
                         "the id, Rev, LB, UB and GA cells (vertically centred) are not visible. "
                         f"The whole row (formula of {m['s12_tokens']} tokens, LB {m['s12_lb']}, "
                         f"UB {m['s12_ub']}, gene rule "
                         f"{'`' + m['s12_gene_rule'] + '`' if m['s12_gene_rule'] else 'none'}) was "
                         f"taken from Table S12 p.{m['s12_page']} with the `BT` prefix removed. "
                         "The visible S10a tokens are an exact prefix of the S12 formula, and this "
                         "is the only S12 reaction absent from S10a."])
        else:
            rd = m.get("readings") or {}
            rows.append([f"S10a p.{m['page']} row {m['row']}", m["reaction"],
                         f"{m['column']} printed `{m['printed']}`; every other S10a bound cell is "
                         f"an integer (-1000, 0 or 1000). It reads as {rd.get('decimal_point'):g} "
                         f"(decimal point) or {rd.get('thousands_separator'):g} (thousands "
                         f"separator). Table S12 p.{m['s12_page']} prints `{m['s12_printed']}` for "
                         f"the same reaction; {m['action']}."])
    L.append(md_table(["Where", "Reaction", "What and why"], rows) + "\n" if rows else "None.\n")
    atpm = v["atpm_lower_bound_readings_glucose_S8c"]
    L.append("Effect of the ATPM reading on growth on glucose minimal medium (Table S8c): "
             + ", ".join(f"{k}: {val}" for k, val in atpm.items())
             + ". The published Table S3 glucose rate is 0.24. The ARAI reading (-1 or -1000) does "
               "not change any growth rate reported here, because the reaction runs forward "
               "(arabinose -> ribulose).\n")
    L.append("No other reaction, bound or gene rule was edited, and no row needed a hand-made "
             "join: every line-break join comes from the general rules in section 1. HYD4 (p. 179) "
             "is the other clipped row; it is recovered complete (only the end of its free-text "
             "note is cut off) and agrees with S12.\n")

    # ---------------------------------------------------------------- parse checks
    L.append("## 3. Parse checks\n")
    L.append(f"- S10a rows read: {rep['extraction']['S10a']['rows']}; without the bottom edge "
             "pdfplumber returns "
             + ", ".join(f"{x['rows_default_settings']} instead of {x['rows_with_bottom_edge']} "
                         f"row(s) on p.{x['page']}" for x in rep["extraction"]["S10a"]["recovered_rows"])
             + ". S10b rows read: " + str(rep["extraction"]["S10b"]["rows"]) + ".")
    L.append(f"- Formula line breaks: {j['formula_line_breaks'].get('space', 0)} at a space and "
             f"{j['formula_line_breaks'].get('joined', 0)} inside a token. Joining every break "
             f"with a space would have produced {j['tokens_a_naive_space_join_would_break']} "
             "invalid tokens.")
    amb = j["formula_breaks_after_a_number_also_valid_as_metabolite"]
    L.append(f"- Breaks right after a number: {j['formula_breaks_after_a_number']}; "
             f"{len(amb)} of them would also give a valid metabolite id if the number were glued "
             "to the next token" + (": " + ", ".join(f"{b['id']}" for b in amb) if amb else "")
             + ".")
    L.append(f"- Gene-rule line breaks: {j['gene_rule_line_breaks'].get('space', 0)} at a space, "
             f"{j['gene_rule_line_breaks'].get('joined', 0)} inside a token.")
    L.append(f"- Free-text joins: {j['text_columns'].get('word_rejoined', 0)} broken words "
             f"re-joined, {j['text_columns'].get('punctuation_join', 0)} joins at hyphens or "
             f"punctuation, {j['text_columns'].get('space', 0) + j['text_columns'].get('space_unverified', 0)} "
             f"spaces ({j['text_columns'].get('space_unverified', 0)} of them next to a word not seen "
             "elsewhere, i.e. not independently confirmed).")
    pf = rep["parse_failures"]
    L.append(f"- **Rows failing validation** (metabolite not in S10b, invalid gene token or rule "
             f"syntax, unreadable row, unresolved bound): **{len(pf)}**"
             + ("." if not pf else ":"))
    for f in pf:
        L.append(f"  - {json.dumps(f, default=str)}")
    L.append(f"- Exchange reactions not written `met[e] <=>`: {len(rep['exchange_orientation_issues'])}. "
             f"Identifier collisions after conversion: {len(rep['id_collisions'])}. Reactions with "
             "the same metabolite on both sides: "
             + (", ".join(x["reaction"] for x in rep["same_metabolite_on_both_sides"]) or "none")
             + ".")
    ri = rep["rev_vs_bounds_disagreements"]
    L.append(f"- Rev / arrow / bounds disagreements: {len(ri)}"
             + ("" if not ri else " (" + "; ".join(
                 f"{x['reaction']}: Rev {x['Rev']}, `{x['arrow']}`, LB {x['LB']:g}, UB {x['UB']:g}"
                 for x in ri) + ")") + "; the published bounds are kept (S12 prints the same).")
    rt = rep["sbml_round_trip"]
    L.append(f"- SBML round trip (write, read back, compare stoichiometry, bounds, gene rule and "
             f"subsystem of every reaction): {'OK' if rt['ok'] else 'FAILED'}; libsbml validation "
             f"messages: {rt.get('libsbml_validation')}.\n")

    L.append("### Cross-check against Table S12\n")
    L.append(md_table(["Check", "Result"], [
        ["Non-exchange reactions compared", cx["reactions_compared"]],
        ["Stoichiometry differs", len(cx["stoichiometry_mismatch"])],
        ["Direction differs", len(cx["arrow_mismatch"])],
        ["Bounds differ", len(cx["bound_mismatch"])],
        ["Gene rule differs", len(cx["gpr_logic_mismatch"])],
        ["Gene rule written differently, same logic", len(cx["gpr_text_differs_logic_same"])],
        ["S10a reactions missing from S12 / S12 reactions missing from S10a",
         f"{len(cx['missing_in_s12'])} / {len(cx['missing_in_s10a'])}"],
        ["Exchanged metabolites in S10a / S12 (identical sets)",
         f"{cx['exchange_metabolites']['s10a']} / {cx['exchange_metabolites']['s12']} "
         f"({'yes' if not cx['exchange_metabolites']['only_in_s10a'] and not cx['exchange_metabolites']['only_in_s12'] else 'no'})"],
        ["S12 rows that failed to parse", len(cx["s12_parse_failures"])],
    ]) + "\n")
    if cx["gpr_logic_mismatch"]:
        L.append("Gene-rule differences (in every case the S10a rule is well formed and is used):\n")
        L.append(md_table(["Reaction", "Problem in S12", "S10a rule valid"],
                          [[x["reaction"], x["s12_problem"], "yes" if x["s10a_rule_valid"] else "no"]
                           for x in cx["gpr_logic_mismatch"]]) + "\n")
    for key in ("stoichiometry_mismatch", "arrow_mismatch", "bound_mismatch"):
        for x in cx[key]:
            L.append(f"- {key}: {json.dumps(x, default=str)}")

    # ---------------------------------------------------------------- biomass
    L.append("## 4. Biomass reaction and objective\n")
    b = rep["biomass"]
    L.append(f"There is one biomass reaction, `{b['objective']}`, and it is the objective. "
             f"Biomass rows in S10a: {b['s10a_biomass_reactions']}; in the S12 *B. theta* part: "
             f"{b['s12_bt_biomass_reactions']}. Evidence:\n")
    for e in b["evidence"]:
        L.append(f"- {e}")
    L.append("- The supplementary text (p. 3, 'Curation of the biomass reaction') describes a "
             "single curated biomass reaction, and the Table S3 growth rates are reproduced with "
             "this objective (section 5b).")
    L.append(f"- Checks run by the script: {b['evidence_checked']}.\n")

    # ---------------------------------------------------------------- validation
    L.append("## 5. Validation against the paper's own model predictions\n")
    L.append("### 5a. Counts\n")
    L.append(md_table(["Quantity", "Paper", "Rebuilt"], [
        ["Reactions", c["published"]["reactions"], rb["reactions"]],
        ["Metabolic and transport", c["published"]["metabolic_and_transport"],
         rb["metabolic_and_transport"]],
        ["Exchange and demand (EX_ + DM_ + sink_)", c["published"]["exchange_and_demand"],
         f"{rb['exchange_and_demand']} ({rb['exchange']} + {rb['demand']} + {rb['sink']})"],
        ["Metabolites (compartment-specific)", c["published"]["metabolites"], rb["metabolites"]],
        ["Genes", c["published"]["genes"], rb["genes"]],
    ]) + "\n")
    L.append(f"S10b lists {rb['s10b_metabolites']} compartment-free metabolites; "
             f"{rb['unique_metabolites_without_compartment']} are used by reactions"
             + (f" (unused: {', '.join(rb['s10b_metabolites_unused'])})"
                if rb["s10b_metabolites_unused"] else "") + ".\n")

    L.append("### 5b. Table S3 growth rates\n")
    L.append("Published = the 'growth rate in silico' column of Table S3 (iAH991's own prediction), "
             "transcribed and checked by the script against the words in that column of the PDF "
             f"({'match' if s3['transcription_matches_pdf_column'] else 'MISMATCH'}). Medium: Table "
             "S8a + the Table S8b carbon source at the printed rate, all other uptakes closed.\n")
    rows = []
    for r in s3["rows"]:
        up = "; ".join(f"{k} {x['rate']:g}" for k, x in r["uptake"].items())
        rows.append([r["condition"], up, r["published"], fmt(r["rebuilt"]), fmt(r["abs_diff"]),
                     "yes" if r["rounds_to_published"] else "**no**"])
    L.append(md_table(["Condition", "Uptake (mmol/gDW/h)", "Published", "Rebuilt", "abs diff",
                       "Rounds to published"], rows) + "\n")
    worst = max((r for r in s3["rows"] if not r["alternative_mapping"]), key=lambda r: r["abs_diff"])
    L.append(f"**{s3['n_rounding_to_published']}/{s3['n']}** published rows are reproduced to the "
             f"printed precision; the largest absolute difference is **{fmt(s3['max_abs_diff'])} "
             f"/h** ({worst['condition']}). The last row is not one of the {s3['n']}: it is the "
             "second possible substrate for the 'Amylopectin' row.\n")
    L.append("Why the three remaining rows differ (diagnostics only; the table above is not "
             "changed by them):\n")
    dv = s3d["s8b_entries_deviating_from_60C_with_S3_row"]
    L.append("- Table S8b doses monosaccharides per molecule (hexoses 10, pentoses 12 mmol/gDW/h) "
             "and oligo-/polysaccharides and glycans at 60 mmol carbon/gDW/h (10 hexose "
             "equivalents). The script checks every S8b rate against 60/(carbon atoms in the S10b "
             "formula). Entries off by more than 2% that have an S3 row:\n")
    L.append(md_table(["Substrate", "S3 row", "Published", "Printed rate", "Growth at printed rate",
                       "60 mmol C rate", "Growth at 60 mmol C rate"],
                      [[d["id"], d["s3_row"], d["published"], f"{d['printed_rate']:g}",
                        fmt(d["growth_at_printed_rate"]), f"{d['rate_60_mmol_C']:g}",
                        fmt(d["growth_at_60_mmol_C"])] for d in dv]) + "\n")
    L.append("  The amino sugars and sialic acid are dosed per molecule like the other "
             "monosaccharides, and the printed rate reproduces the published value. Two entries "
             "break the rule of their own group: L-rhamnose (a 6-carbon deoxy sugar like L-fucose, "
             "which is dosed at 10) is printed at the pentose rate 12, and GlcNAc-Core 2 (C30) is "
             "printed at 2.7273, the rate of GlcNAc-Core 1 (C22), whereas every other mucin glycan "
             "follows 60/C exactly. At the group-consistent rates the rebuilt model gives "
             + "; ".join(f"{d['s3_row']} {fmt(d['growth_at_60_mmol_C'])} (published {d['published']})"
                         for d in dv if d["60C_rate_reproduces"] and not d["printed_rate_reproduces"])
             + ". Growth this close to zero is very sensitive to the uptake rate because the ATP "
               "maintenance requirement (ATPM 8.43) has to be met first.")
    amyl = [d for d in s3d["non_reproducing_rows"] if d["condition"] == "Amylopectin"]
    if amyl:
        a = amyl[0]
        alt = "; ".join(f"{x['what']}: {fmt(x['rebuilt'])}" for x in a["alternatives"])
        L.append(f"- Amylopectin: S3 has one 'Amylopectin' row (in vivo column '0.075/0.0056') "
                 "but S8b lists two starch substrates, strch1 (a branched 11-glucose structure, "
                 "used here) and starch1200 (potato starch, 900 amylopectin + 300 amylose units). "
                 f"strch1 gives {fmt(a['rebuilt'])}; {alt} (published {a['published']}).")
    L.append(f"- Accepting these readings, {s3d['n_rows_reproduced_if_alternatives_accepted']}/"
             f"{s3['n']} rows reproduce and the largest difference is "
             f"{fmt(s3d['max_abs_diff_if_alternatives_accepted'])} /h.")
    sens = s3["sensitivity_sinks_closed"]
    roles = "; ".join(
        f"`{k}` (flux {x['flux_at_optimum_glucose']:g} at the glucose optimum; other reactions of "
        f"its metabolite: {', '.join(x['other_reactions_of_metabolite'])}; growth on glucose with "
        f"its uptake closed: {fmt(x['growth_glucose_if_uptake_closed'])})"
        for k, x in sens["sink_roles"].items())
    L.append(f"- Sinks: the two sink reactions are not exchange reactions, so Table S8 does not "
             f"cover them, and they are left open as published. They supply required precursors: "
             f"{roles}. If their uptake is also closed, {sens['n_published_growth_lost']} of the "
             f"{sum(1 for r in sens['rows'] if float(r['published']) > 0)} conditions with "
             f"published growth no longer grow (largest difference "
             f"{fmt(sens['max_abs_diff'])} /h), so the paper's simulations must have kept them "
             "open.\n")

    L.append("Table S3b secretion ranges (supplementary check; " + s3b["method_assumed"] + "):\n")
    rows = []
    for r in s3b["rows"]:
        for x in ("ac", "ppa", "succ", "co2", "h2"):
            lo, hi = r["rebuilt"][x]
            val = f"{lo:.2f}" if abs(lo - hi) < 0.005 else f"{lo:.2f} to {hi:.2f}"
            rows.append([r["substrate"], x, r["published"][x], val,
                         "yes" if r["matches"][x] else "**no**"])
    L.append(md_table(["Substrate", "Product", "Published", "Rebuilt", "Match"], rows) + "\n")
    miss = [f"{r['substrate']} {x} (published {r['published'][x]}, rebuilt "
            f"{r['rebuilt'][x][0]:.2f} to {r['rebuilt'][x][1]:.2f})"
            for r in s3b["rows"] for x, ok in r["matches"].items() if not ok]
    L.append(f"{s3b['n_values_matching']}/{s3b['n_values']} values match to two decimals."
             + (f" Not matching: {'; '.join(miss)}. This is unexplained: the growth rate and the "
                "other secretion ranges of the same condition do match, so a misprint is possible "
                "but cannot be shown." if miss else "") + "\n")

    L.append("### 5c. Table S1 sole carbon sources\n")
    L.append(f"Set-up: {s1['setup']}. Id alias: {s1['id_aliases_used']} (S1 prints `ara-L`, which "
             "is not in S10b; S8b and S10b use `arab-L`).\n")
    rows = []
    for r in s1["rows"]:
        mem = ", ".join(f"{m['id']} {fmt(m['growth'])}" for m in r["members"])
        rows.append([r["carbohydrate"], mem, r["published_iAH991"], r["rebuilt"],
                     "yes" if r["agrees"] else "**no**"])
    L.append(md_table(["Carbohydrate", "Id and rebuilt growth (/h)", "Published iAH991",
                       "Rebuilt", "Agrees"], rows) + "\n")
    L.append(f"Rows agreeing: **{s1['rows_agreeing']}/{s1['n_rows']}** (a row listing several ids "
             f"counts as growth if any of them supports growth). Individual ids supporting growth: "
             f"{s1['ids_growing']}/{s1['n_ids']}."
             + (f" S1 prints `core7` twice in the mucin row; core8 (probably meant) gives "
                f"{fmt(s1['core8_extra']['growth'])} /h." if s1.get("core8_extra") else ""))
    for d in s1["diagnostics_non_growing_ids"]:
        L.append(f"- `{d['id']}` does not grow at the set-up rate; growth at 10 mmol/gDW/h: "
                 f"{fmt(d['growth_at_10'])}, at 100: {fmt(d['growth_at_100'])}; reactions using it: "
                 f"{', '.join(d.get('reactions_using_it', []))}.")
    L.append("")

    L.append("### 5d. Gene essentiality\n")
    L.append(f"Definition: {ge['definition']}. Published counts: {ge['published_counts_source']}.\n")
    rows = []
    for name, r in ge["media"].items():
        rows.append([name, fmt(r["wild_type_growth"]), r["published"],
                     r["essential_growth_below_1e-6"], r["essential_below_1pct_of_wt"],
                     r["essential_below_5pct_of_wt"], r["essential_below_10pct_of_wt"],
                     r["essential_below_50pct_of_wt"]])
    L.append(md_table(["Medium", "WT growth (/h)", "Published essential", "Rebuilt (< 1e-6)",
                       "< 1% WT", "< 5% WT", "< 10% WT", "< 50% WT"], rows) + "\n")
    g = ge["media"]["glucose_minimal_S8c"]
    conf = g["s14_reports_viable_but_essential_here"]
    L.append(f"- Glucose minimal medium: {g['essential_growth_below_1e-6']} vs {g['published']}. "
             "Table S14 states that the BT_3763-BT_3767 deletion grows in silico on glucose "
             "minimal medium, but in the rebuild "
             + (f"{', '.join(conf)} {'is' if len(conf) == 1 else 'are'} essential there (see 5e). "
                f"Without {'it' if len(conf) == 1 else 'them'} the count would be "
                f"{g['count_without_those']}." if conf else "none of those genes is essential."))
    t = ge["media"]["TYG_S8d"]
    hits = [x for x in t["exploratory_vitamin_additions"] if x["essential"] == t["published"]]
    L.append(f"- TYG: {t['essential_growth_below_1e-6']} vs {t['published']}; the count does not "
             "depend on the threshold. The printed Table S8d list has no biotin, pantothenate, "
             "nicotinate or folate, although iAH991 has exchange reactions for them, and the S10a "
             "note to BTNt2 says biotin is taken up from yeast extract. Exploratory only, not "
             "adopted: essential-gene counts when these are added at 10 mmol/gDW/h:\n")
    L.append(md_table(["Added to S8d", "Essential", "No longer essential"],
                      [[", ".join(x["added_at_10"]), x["essential"],
                        ", ".join(x["no_longer_essential"])] for x in t["exploratory_vitamin_additions"]])
             + "\n")
    L.append(f"  {len(hits)} of these {len(t['exploratory_vitamin_additions'])} additions happen "
             "to give the published count. The paper does not describe the medium beyond Table "
             "S8d, so the TYG difference is left unexplained.")
    r_ = ge["media"]["rich"]
    L.append(f"- Rich medium (every exchange open at -1000): {r_['essential_growth_below_1e-6']} vs "
             f"{r_['published']}.")
    w = ge["media"]["western_diet_S9"]
    L.append(f"- Western diet: {w['essential_growth_below_1e-6']} vs {w['published']} (main text). "
             f"{ge['western_diet_note']} Table S9 ids without an iAH991 exchange (ignored): "
             f"{', '.join(ge['western_diet_unmapped_ids']) or 'none'}.\n")

    L.append("### 5e. Table S14 knockout predictions\n")
    L.append("Same carbon-source set-up as 5c; 'glucose minimal medium' = Table S8c; L-arginine "
             "supplied at 10 mmol/gDW/h; 'sulfated disaccharides' = cspg_a_degr, cspg_b_degr, "
             "cspg_c_degr; 'starch' = starch1200. The published in-silico column was checked against "
             f"the PDF ({'match' if s14['transcription_matches_pdf_column'] else 'MISMATCH'}).\n")
    rows = []
    for r in s14["rows"]:
        detail = "; ".join(
            f"{'+'.join(cc['deleted']) or '-'} on {'/'.join(cc['substrates'])}: {fmt(cc['growth'])}"
            + (" (no exchange)" if cc["missing_exchange"] else "") for cc in r["cases"])
        rows.append([r["deleted"], r["condition"], r["published_in_silico"], r["rebuilt"],
                     "yes" if r["agrees"] else "**no**", detail])
    L.append(md_table(["Deleted", "Condition", "Published in silico", "Rebuilt", "Agrees",
                       "Growth per case (/h)"], rows) + "\n")
    L.append(f"Agreement: **{s14['n_agree']}/{s14['n']}** rows. Disagreements:\n")
    for r in s14["rows"]:
        if r["agrees"]:
            continue
        for e in r["explanation"]:
            if "single_deletions_that_are_lethal" in e:
                parts = []
                for le in e["single_deletions_that_are_lethal"]:
                    rx = "; ".join(f"{x['reaction']} (`{x['equation']}`)" for x in le["lethal_reactions"])
                    parts.append(f"{le['gene']} alone is lethal through {rx or 'no single reaction'}")
                L.append(f"- {r['deleted']} on {'/'.join(e['substrates'])}: published grows, rebuilt "
                         f"does not. {'; '.join(parts)}.")
            elif "reactions_of_deleted_genes_still_carrying_flux" in e:
                parts = [f"{x['reaction']} (gene rule `{x['gene_rule']}`)"
                         for x in e["reactions_of_deleted_genes_still_carrying_flux"]]
                L.append(f"- {'+'.join(e['deleted'])} on {'/'.join(e['substrates'])}: published "
                         f"no growth, rebuilt grows; still active: {'; '.join(parts)}.")
    L.append("")

    # ---------------------------------------------------------------- conclusions
    L.append("## 6. Outcome\n")
    L.append(
        f"- The model content is read without a single token-validation failure, and the "
        f"independent S12 rendering agrees with S10a on the stoichiometry, direction and bounds of "
        f"all {cx['reactions_compared']} reactions compared; the only differences are "
        f"{len(cx['gpr_logic_mismatch'])} gene rules that are malformed or empty in S12. Three values "
        "rest on S12 instead: the clipped biomass row, and the ARAI and ATPM bound cells, for which "
        "S12 was the tie-breaker.")
    L.append(f"- Counts match the paper exactly ({rb['reactions']} reactions, {rb['metabolites']} "
             f"metabolites, {rb['genes']} genes).")
    L.append(f"- Reproduced: Table S3 growth {s3['n_rounding_to_published']}/{s3['n']} as set up "
             f"(largest difference {fmt(s3['max_abs_diff'])} /h; "
             f"{s3d['n_rows_reproduced_if_alternatives_accepted']}/{s3['n']} and "
             f"{fmt(s3d['max_abs_diff_if_alternatives_accepted'])} /h with the readings in 5b); "
             f"Table S3b {s3b['n_values_matching']}/{s3b['n_values']}; Table S1 "
             f"{s1['rows_agreeing']}/{s1['n_rows']} rows; Table S14 {s14['n_agree']}/{s14['n']} "
             f"rows; essential genes in rich medium {r_['essential_growth_below_1e-6']} vs "
             f"{r_['published']}.")
    open_items = []
    for d in s3d["non_reproducing_rows"]:
        how = ("S8b uptake rate" if any("60 mmol C" in a["what"] for a in d["alternatives"])
               else "substrate mapping")
        ok = any(a["rounds_to_published"] for a in d["alternatives"])
        open_items.append(f"Table S3 {d['condition']}: {fmt(d['rebuilt'])} vs {d['published']} "
                          + (f"(reproduced with the alternative {how}, 5b)" if ok
                             else "(unexplained)"))
    for r in s3b["rows"]:
        for x, ok in r["matches"].items():
            if not ok:
                open_items.append(f"Table S3b {r['substrate']} {x} maximum: "
                                  f"{r['rebuilt'][x][1]:.2f} vs {r['published'][x]} (unexplained)")
    for r in s1["rows"]:
        if not r["agrees"]:
            open_items.append(f"Table S1 {r['carbohydrate']}: no growth at the set-up rate "
                              "(grows at 10 mmol/gDW/h, see 5c)")
    for r in s14["rows"]:
        if not r["agrees"]:
            open_items.append(f"Table S14 {r['deleted']} / {r['condition']}: published "
                              f"{r['published_in_silico']}, rebuilt {r['rebuilt']} (5e)")
    open_items.append(f"essential genes on glucose minimal medium {g['essential_growth_below_1e-6']} "
                      f"vs {g['published']} ("
                      + (f"{', '.join(conf)} is essential here although Table S14 implies it is "
                         f"not; without it the count would be {g['count_without_those']}"
                         if conf else "no candidate gene identified") + ")")
    open_items.append(f"essential genes on TYG {t['essential_growth_below_1e-6']} vs "
                      f"{t['published']} (unexplained; possibly the TYG medium, 5d)")
    open_items.append(f"essential genes on the Western diet {w['essential_growth_below_1e-6']} vs "
                      f"{w['published']} (approximate set-up)")
    L.append("- Not reproduced:")
    for item in open_items:
        L.append(f"  - {item}")
    L.append("- None of these is caused by a parse error that the S12 cross-check could see: the "
             "two renderings print the same network. Where an explanation was found it lies either "
             "in the condition tables (two S8b rates, the Amylopectin row) or in the printed network "
             "itself: LCAR2 (BT_3767) is the only consumer of the glycolaldehyde made in thiamine "
             "and folate synthesis, and the printed AMYe rule lets BT_0773 or BT_4305 replace SusG. "
             "That suggests the authors simulated a slightly different network or medium than the "
             "one printed, but the PDF cannot settle this.\n")

    L.append("## 7. Known limitations\n")
    L.append(LIMITATIONS)
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


LIMITATIONS = """\
- The biomass reaction cannot be read from Table S10a (the row is clipped by the page edge on
  p. 89). Its full formula and bounds come from Table S12 of the same PDF (joint model). The
  visible S10a part agrees token for token, but the rest of the formula is attested only by S12.
  Its S10a name, subsystem and notes are not visible and are left empty (name = id).
- The ATPM lower bound is printed `84.300` in S10a (and also `84.300` for the mouse model in S11a)
  but `8,43` in S12. 8.43 is used. 84.3 makes the model infeasible on glucose minimal medium and
  84300 exceeds the upper bound, while 8.43 reproduces Table S3; still, the S10a cell itself is
  unreadable.
- HYD4's free-text note (p. 179) ends mid-sentence because that row is also clipped by the page
  edge; its formula, gene rule and bounds are complete and agree with S12.
- Names, subsystems and curation notes are free text joined by heuristic rules. A word broken
  across lines may still carry a spurious space where the whole word does not occur elsewhere in
  the PDF. This does not affect the model's mathematics.
- Metabolite formulas and charges are copied from S10b; mass and charge balance were not checked.
- FORt is printed reversible (Rev 1, `<=>`) with LB 0 in both S10a and S12; the bounds are kept.
- Table S8 only addresses exchange reactions, so the two sink reactions stay open as published
  (-1000..1000) in every simulated medium. They supply elemental sulfur for biotin synthase and
  choline sulfate (the only choline source), and closing their uptake stops growth in the
  Table S3 conditions (section 5b), so the paper must have kept them open too. A literal reading
  of "all other uptakes closed" that included the sinks would make every published growth rate
  unreachable.
- The paper does not state its gene-essentiality threshold. Counts are given at several
  thresholds and barely change.
- The Western-diet count applies the joint-model diet of Table S9 to iAH991 alone, which is not
  the paper's setting.
- Table S3 has one "Amylopectin" row but Table S8b lists two starch substrates; strch1 was chosen
  before running, and starch1200 is reported as the alternative.
- Gene rules are as printed, including the duplicated `(BT_1683 and BT_1683)` in
  RHAMNOGALURASEe_I/II.
"""


if __name__ == "__main__":
    main()
