"""Check 1b: 30 reactions of the rebuilt iAH991 read a third way, with poppler's pdftotext (a different PDF text
engine from pdfplumber, which both the rebuild and reparse_pdf_tables.py use).

For each reaction the row's Formula and GA cells are cut out of the page (cell borders from the reparse JSON) and
extracted with `pdftotext -layout -x -y -W -H`. Every space is removed and the result compared with the model's
reaction written in the paper's notation (also with every space removed, metabolites in the model's order). Gene
rules are compared the same way after removing parentheses-only differences (as Boolean functions).

Selection: the 10 reactions with the most metabolites, the 5 with the longest gene rules, 5 named cases (HYD4,
RHAMNOGALURASEe_I, AMYe, GLYASNt, HSK) and 10 at random (seed 20261006), excluding the biomass row.

Usage: python -I spot_check_poppler.py <pdf> <reparse.json>   (writes spot_check_poppler.json)
"""
from __future__ import annotations

import json
import os
import random
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT  # noqa: E402
from compare_rebuild_to_pdf import dnf_tokens, gpr_tokens, model_dnf, model_id  # noqa: E402

import cobra  # noqa: E402

OUT = os.path.dirname(os.path.abspath(__file__))


def crop_text(pdf, page, x0, x1, top, bottom):
    cmd = ["pdftotext", "-f", str(page), "-l", str(page), "-layout", "-x", str(int(x0) + 1), "-y", str(int(top)),
           "-W", str(int(x1 - x0) - 2), "-H", str(int(bottom - top) + 3), pdf, "-"]
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def printed(model, r):
    def term(m, c):
        base, comp = m.id.rsplit("_", 1)
        name = base.replace("__", "-")
        c = abs(c)
        coef = "" if c == 1 else (f"{c:g}" if c != int(c) else str(int(c)))
        return f"{coef}{name}[{comp}]"
    lhs = "+".join(term(m, c) for m, c in r.metabolites.items() if c < 0)
    rhs = "+".join(term(m, c) for m, c in r.metabolites.items() if c > 0)
    arrow = "<=>" if r.lower_bound < 0 < r.upper_bound else "=>"
    return f"{lhs}{arrow}{rhs}"


def main():
    pdf, rp = sys.argv[1], json.load(open(sys.argv[2]))
    rows = [x for x in rp["36-262"] if "cells" in x and not x["cells"][0]["despaced"].startswith("SEED") and x["cells"][1]["despaced"]]
    by_id = {model_id(x["cells"][1]["despaced"]): x for x in rows}
    model = cobra.io.read_sbml_model(os.path.join(ROOT, "models/curated/iAH991/iAH991_rebuilt.xml"))
    rxns = [r for r in model.reactions if r.id in by_id]
    most_mets = sorted(rxns, key=lambda r: -len(r.metabolites))[:10]
    longest_gpr = sorted([r for r in rxns if r not in most_mets], key=lambda r: -len(r.gene_reaction_rule))[:5]
    named = [model.reactions.get_by_id(i) for i in ("HYD4", "RHAMNOGALURASEe_I", "AMYe", "GLYASNt", "HSK")]
    rest = [r for r in rxns if r not in most_mets + longest_gpr + named]
    rnd = random.Random(20261006).sample(rest, 10)
    out = []
    for group, lst in (("most metabolites", most_mets), ("longest gene rule", longest_gpr), ("named", named), ("random", rnd)):
        for r in lst:
            x = by_id[r.id]
            xs = x["column_borders"]
            ftxt = crop_text(pdf, x["page"], xs[3], xs[4], x["top"], x["bottom"])
            gtxt = crop_text(pdf, x["page"], xs[5], xs[6], x["top"], x["bottom"])
            fdesp = re.sub(r"\s+", "", ftxt)
            gdesp = re.sub(r"\s+", "", gtxt)
            mform = printed(model, r)
            toks = gpr_tokens(gdesp)
            try:
                gequal = toks is not None and dnf_tokens(toks) == model_dnf(r)
            except (ValueError, IndexError):
                gequal = False
            out.append({"group": group, "reaction": r.id, "page": x["page"], "n_metabolites": len(r.metabolites),
                        "formula_equal": fdesp == mform, "gene_rule_equal": gequal,
                        "pdf_formula": fdesp, "model_formula": mform, "pdf_gene_rule": gdesp, "model_gene_rule": r.gene_reaction_rule})
    summary = {"n": len(out), "formula_equal": sum(o["formula_equal"] for o in out), "gene_rule_equal": sum(o["gene_rule_equal"] for o in out)}
    with open(os.path.join(OUT, "spot_check_poppler.json"), "w") as fh:
        json.dump({"summary": summary, "reactions": out}, fh, indent=1)
    print(summary)
    for o in out:
        if not (o["formula_equal"] and o["gene_rule_equal"]):
            print("DIFF", o["reaction"], o["page"], "\n  pdf  ", o["pdf_formula"][:400], "\n  model", o["model_formula"][:400],
                  "\n  pdf gpr  ", o["pdf_gene_rule"][:300], "\n  model gpr", o["model_gene_rule"][:300])
    print([(o["group"], o["reaction"], o["n_metabolites"]) for o in out])


if __name__ == "__main__":
    main()
