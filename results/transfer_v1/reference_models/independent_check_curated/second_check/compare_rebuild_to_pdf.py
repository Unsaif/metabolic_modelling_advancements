"""Check 1a: every reaction of models/curated/iAH991/iAH991_rebuilt.xml against this check's own reading of
Tables S10a, S10b and S12 (reparse_pdf_tables.py), not only a sample.

For every S10a row: identifier, stoichiometry (from the despaced formula, split into terms with the S10b
metabolite list), direction, lower and upper bound, and gene rule (as a Boolean function) are compared with the
rebuilt model. The biomass row (cut off on p. 89) is compared with S12, and its visible S10a part with the model.
Every S12 B. theta reaction is compared too (formula, bounds), to see whether S12 could have hidden S10a errors.

Usage: python -I compare_rebuild_to_pdf.py <reparse.json>   (writes compare_rebuild_to_pdf.json)
"""
from __future__ import annotations

import ast
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import ROOT  # noqa: E402

import cobra  # noqa: E402

OUT = os.path.dirname(os.path.abspath(__file__))
TOL = 1e-6


def model_id(published):
    s = published.replace("(e)", "_e").replace("-", "__")
    return re.sub(r"[^A-Za-z0-9_]", "_", s)


def met_model_id(printed):
    m = re.match(r"^(.*)\[(\w+)\]$", printed)
    base, comp = m.group(1), m.group(2)
    return f"{base.replace('-', '__')}_{comp}"


def split_terms(side):
    return [t for t in side.split("+") if t] if side else []


NUM = r"\(?\d+(?:\.\d+)?(?:[eE]-?\d+)?\)?"


def parse_formula(desp, mets, prefix=""):
    """desp: formula with every space removed. mets: set of printed metabolite ids (no compartment).
    Returns (stoich {printed_id_with_comp: coef}, arrow, ambiguities, errors)."""
    for arrow in ("<=>", "=>", "->", "<=", "<-"):
        if arrow in desp:
            lhs, rhs = desp.split(arrow, 1)
            break
    else:
        return None, None, [], [f"no arrow in {desp!r}"]
    st, amb, err = {}, [], []
    for sign, side in ((-1.0, lhs), (1.0, rhs)):
        for term in split_terms(side):
            m = re.match(r"^(.*)\[(\w+)\]$", term)
            if not m:
                err.append(f"term without compartment: {term!r}")
                continue
            body, comp = m.group(1), m.group(2)
            cands = []
            for i in range(len(body)):
                pre, met = body[:i], body[i:]
                if prefix:
                    if not met.startswith(prefix):
                        continue
                    met = met[len(prefix):]
                if met in mets and (pre == "" or re.fullmatch(NUM, pre)):
                    cands.append((pre, met))
            if not cands:
                err.append(f"cannot split {term!r}")
                continue
            if len(cands) > 1:
                amb.append({"term": term, "readings": cands})
            pre, met = cands[0]
            coef = float(pre.strip("()")) if pre else 1.0
            key = f"{met}[{comp}]"
            st[key] = st.get(key, 0.0) + sign * coef
    return st, arrow, amb, err


def gpr_tokens(desp):
    toks = re.findall(r"BT_\d{4}|and|or|\(|\)", desp)
    return toks if "".join(toks) == desp else None


def dnf_tokens(toks):
    """Disjunctive normal form of a token list (and binds tighter than or)."""
    pos = [0]

    def expr():
        terms = [term()]
        while pos[0] < len(toks) and toks[pos[0]] == "or":
            pos[0] += 1
            terms.append(term())
        out = set()
        for t in terms:
            out |= t
        return out

    def term():
        facs = [factor()]
        while pos[0] < len(toks) and toks[pos[0]] == "and":
            pos[0] += 1
            facs.append(factor())
        out = {frozenset()}
        for f in facs:
            out = {a | b for a in out for b in f}
        return out

    def factor():
        t = toks[pos[0]]
        if t == "(":
            pos[0] += 1
            v = expr()
            if pos[0] >= len(toks) or toks[pos[0]] != ")":
                raise ValueError("unbalanced")
            pos[0] += 1
            return v
        pos[0] += 1
        return {frozenset([t])}

    if not toks:
        return frozenset()
    v = expr()
    if pos[0] != len(toks):
        raise ValueError("trailing tokens")
    return frozenset(s for s in v if not any(o < s for o in v))


def dnf_ast(node):
    if node is None:
        return frozenset()
    if isinstance(node, ast.Name):
        return frozenset([frozenset([node.id])])
    if isinstance(node, ast.BoolOp):
        parts = [dnf_ast(v) for v in node.values]
        if isinstance(node.op, ast.Or):
            out = set()
            for s in parts:
                out |= s
            return frozenset(out)
        out = {frozenset()}
        for s in parts:
            out = {a | b for a in out for b in s}
        return frozenset(out)
    raise TypeError(type(node))


def model_dnf(r):
    d = dnf_ast(getattr(r.gpr, "body", None))
    return frozenset(s for s in d if not any(o < s for o in d))


def fnum(s):
    s = s.replace(",", ".")
    try:
        return float(s)
    except ValueError:
        return None


def main():
    rp = json.load(open(sys.argv[1]))
    s10a = [x for x in rp["36-262"] if "cells" in x and not x["cells"][0]["despaced"].startswith("SEED")]
    s10b = [x for x in rp["263-306"] if "cells" in x and not x["cells"][0]["despaced"].startswith("SEED")]
    s12 = [x for x in rp["472-490"] if "cells" in x and x["cells"][0]["despaced"] != "ReactionID"]
    printed_mets = {x["cells"][1]["despaced"]: {"seed": x["cells"][0]["despaced"], "name": " ".join(x["cells"][2]["lines"]),
                                                "formula": x["cells"][3]["despaced"], "charge": x["cells"][4]["despaced"]}
                    for x in s10b}
    mets = set(printed_mets)
    model = cobra.io.read_sbml_model(os.path.join(ROOT, "models/curated/iAH991/iAH991_rebuilt.xml"))
    rep = {"s10a_rows": len(s10a), "s10b_rows": len(s10b), "s12_rows": len(s12),
           "model_counts": {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)}}
    # metabolite id map
    base_of = {}
    for m in model.metabolites:
        mm = re.match(r"^(.*)_(\w)$", m.id)
        base_of[m.id] = mm.group(1)
    conv = {p: p.replace("-", "__") for p in mets}
    inv = {}
    for p, b in conv.items():
        inv.setdefault(b, []).append(p)
    unmatched_model = sorted({b for b in base_of.values() if b not in inv})
    rep["metabolites"] = {"s10b_ids": len(mets), "model_bases": len(set(base_of.values())),
                          "model_bases_without_s10b_id": unmatched_model,
                          "s10b_ids_unused_by_model": sorted(p for p, b in conv.items() if b not in set(base_of.values())),
                          "s10b_ids_with_plus": sorted(p for p in mets if "+" in p),
                          "non_injective": {b: v for b, v in inv.items() if len(v) > 1}}

    def to_printed(stoich):
        out = {}
        for mid, c in stoich.items():
            b, comp = base_of[mid], mid.rsplit("_", 1)[1]
            out[f"{inv[b][0]}[{comp}]"] = c
        return out

    probs = {"id_not_in_model": [], "formula_parse": [], "stoichiometry": [], "direction": [], "bounds": [],
             "gene_rule_parse": [], "gene_rule": [], "ambiguous_split": []}
    seen = set()
    biomass_row = None
    for x in s10a:
        c = x["cells"]
        pid, form, rev, ga, lb, ub = (c[1]["despaced"], c[3]["despaced"], c[4]["despaced"], c[5]["despaced"],
                                      c[6]["despaced"], c[7]["despaced"])
        where = f"p.{x['page']} row {x['row_on_page']}"
        if not pid:
            biomass_row = {"where": where, "visible_formula": form, "open_bottom": x["open_bottom"]}
            continue
        rid = model_id(pid)
        if rid not in model.reactions:
            probs["id_not_in_model"].append({"where": where, "published": pid, "converted": rid})
            continue
        seen.add(rid)
        r = model.reactions.get_by_id(rid)
        st, arrow, amb, err = parse_formula(form, mets)
        if err:
            probs["formula_parse"].append({"where": where, "id": pid, "formula": form, "errors": err})
            continue
        if amb:
            probs["ambiguous_split"].append({"where": where, "id": pid, "ambiguities": amb})
        mst = to_printed({m.id: v for m, v in r.metabolites.items()})
        if set(st) != set(mst) or any(abs(st[k] - mst[k]) > TOL * max(1, abs(st[k])) for k in st):
            probs["stoichiometry"].append({"where": where, "id": pid, "pdf": st, "model": mst})
        rev_pdf = arrow == "<=>"
        if rev_pdf != (r.lower_bound < 0 < r.upper_bound) or rev != ("1" if rev_pdf else "0"):
            probs["direction"].append({"where": where, "id": pid, "arrow": arrow, "rev_cell": rev, "model_bounds": r.bounds})
        lbv, ubv = fnum(lb), fnum(ub)
        if lbv is None or ubv is None or abs(lbv - r.lower_bound) > TOL or abs(ubv - r.upper_bound) > TOL:
            probs["bounds"].append({"where": where, "id": pid, "pdf": [lb, ub], "model": list(r.bounds)})
        toks = gpr_tokens(ga)
        if toks is None:
            probs["gene_rule_parse"].append({"where": where, "id": pid, "ga": ga})
        else:
            try:
                d = dnf_tokens(toks)
            except ValueError as e:
                probs["gene_rule_parse"].append({"where": where, "id": pid, "ga": ga, "error": str(e)})
            else:
                if d != model_dnf(r):
                    probs["gene_rule"].append({"where": where, "id": pid, "pdf": ga, "model": r.gene_reaction_rule})
    rep["s10a_rows_compared"] = len(seen)
    rep["model_reactions_without_s10a_row"] = sorted(r.id for r in model.reactions if r.id not in seen)
    rep["problems"] = probs
    rep["problem_counts"] = {k: len(v) for k, v in probs.items()}
    # biomass
    bm = [r for r in model.reactions if r.id not in seen]
    rep["biomass"] = {"s10a_row": biomass_row}
    if biomass_row and len(bm) == 1:
        r = bm[0]
        mst = to_printed({m.id: v for m, v in r.metabolites.items()})
        # S12 biomass row
        s12bm = [x for x in s12 if x["cells"][0]["despaced"] in ("BTBiomass_BT_v2",)]
        rep["biomass"]["model_reaction"] = r.id
        rep["biomass"]["model_n_metabolites"] = len(r.metabolites)
        rep["biomass"]["model_bounds"] = list(r.bounds)
        rep["biomass"]["model_gene_rule"] = r.gene_reaction_rule
        if s12bm:
            f12 = s12bm[0]["cells"][1]["despaced"]
            st12, arrow12, amb12, err12 = parse_formula(f12, mets, prefix="BT")
            rep["biomass"]["s12"] = {"page": s12bm[0]["page"], "lb_ub": [s12bm[0]["cells"][2]["despaced"], s12bm[0]["cells"][3]["despaced"]],
                                     "gpr": s12bm[0]["cells"][4]["despaced"], "arrow": arrow12, "parse_errors": err12,
                                     "n_terms": len(st12 or {}),
                                     "equal_to_model": bool(st12) and set(st12) == set(mst) and all(abs(st12[k] - mst[k]) <= 1e-9 for k in st12)}
            if st12 and not rep["biomass"]["s12"]["equal_to_model"]:
                rep["biomass"]["s12"]["diff"] = {k: [st12.get(k), mst.get(k)] for k in set(st12) | set(mst)
                                                 if abs((st12.get(k) or 0) - (mst.get(k) or 0)) > 1e-9}
        # visible S10a prefix: rebuild the despaced model formula in the model's metabolite order and compare prefix
        vis = biomass_row["visible_formula"]
        st_vis, _, _, err_vis = parse_formula(vis + "=>", mets) if "=>" not in vis else parse_formula(vis, mets)
        # the visible part has no arrow if cut before it; parse the reactant terms directly
        terms = [t for t in vis.split("=>")[0].split("+") if t]
        last_complete = terms if vis.endswith("]") else terms[:-1]
        ok, bad = 0, []
        for t in last_complete:
            s1, _, _, e1 = parse_formula(t + "=>", mets)
            if e1 or not s1:
                bad.append(t)
                continue
            (k, v), = s1.items()
            if k in mst and abs(mst[k] - v) <= 1e-9:
                ok += 1
            else:
                bad.append({"term": t, "pdf": v, "model": mst.get(k)})
        rep["biomass"]["s10a_visible_terms_checked"] = ok
        rep["biomass"]["s10a_visible_terms_disagreeing"] = bad
    # S12 against the model (non-exchange reactions; exchange bounds in S12 describe the lumen)
    s12cmp = {"compared": 0, "stoichiometry": [], "bounds": [], "not_in_model": [], "parse": []}
    for x in s12:
        rid12 = x["cells"][0]["despaced"]
        if not rid12.startswith("BT"):
            continue
        rid = model_id(rid12[2:])
        f12 = x["cells"][1]["despaced"]
        if rid not in model.reactions:
            s12cmp["not_in_model"].append(rid12)
            continue
        r = model.reactions.get_by_id(rid)
        if r.boundary:
            continue
        st12, arrow12, amb12, err12 = parse_formula(f12, mets, prefix="BT")
        if err12:
            s12cmp["parse"].append({"id": rid12, "errors": err12})
            continue
        s12cmp["compared"] += 1
        mst = to_printed({m.id: v for m, v in r.metabolites.items()})
        if set(st12) != set(mst) or any(abs(st12[k] - mst[k]) > TOL * max(1, abs(st12[k])) for k in st12):
            s12cmp["stoichiometry"].append({"id": rid12, "s12": st12, "model": mst})
        lb12, ub12 = fnum(x["cells"][2]["despaced"]), fnum(x["cells"][3]["despaced"])
        if lb12 is None or ub12 is None or abs(lb12 - r.lower_bound) > TOL or abs(ub12 - r.upper_bound) > TOL:
            s12cmp["bounds"].append({"id": rid12, "s12": [x["cells"][2]["despaced"], x["cells"][3]["despaced"]], "model": list(r.bounds)})
    rep["s12_vs_model"] = {k: (v if not isinstance(v, list) else {"n": len(v), "items": v[:20]}) for k, v in s12cmp.items()}
    with open(os.path.join(OUT, "compare_rebuild_to_pdf.json"), "w") as fh:
        json.dump(rep, fh, indent=1, default=str)
    print(json.dumps({k: rep[k] for k in ("s10a_rows", "s10b_rows", "s10a_rows_compared", "model_counts", "problem_counts",
                                          "model_reactions_without_s10a_row")}, indent=1))
    print("metabolites:", {k: (v if not isinstance(v, (list, dict)) else (len(v), list(v)[:10])) for k, v in rep["metabolites"].items()})
    print("biomass:", json.dumps({k: v for k, v in rep["biomass"].items() if k != "s10a_row"}, default=str)[:1500])
    print("s12:", json.dumps({k: (v if not isinstance(v, dict) else v["n"]) for k, v in rep["s12_vs_model"].items()}))


if __name__ == "__main__":
    main()
