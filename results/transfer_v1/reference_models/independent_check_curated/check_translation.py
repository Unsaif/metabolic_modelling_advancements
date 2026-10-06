"""Check 1: is the translation relabelling only?

For iSO783 and iGD1575, compare the published SBML (as read by cobrapy) with the BiGG view reaction by reaction,
mapping identifiers back with the renames in *_translation.json:
  * stoichiometry, bounds, gene rules (as Boolean functions), objective coefficients;
  * iGD1575's two-step exchanges (X_e0 <=> X_b, X_b <=>) against the single X_e <=> of the view, with the
    orientation and the bounds checked (and that the dropped boundary species took part in nothing else);
  * every published reaction, metabolite and gene accounted for;
  * optimal growth with both versions under: published bounds; every boundary open; each mapped carbon-source
    condition of the organism's medium applied as the protocol applies it (no medium completion); 30 random media;
  * 60 single-gene deletions in one condition.
The published iGD1575 gene rules are also re-read from the SBML notes and compared with what cobrapy made of them.

Usage: python -I check_translation.py   (writes check_translation.json)
"""
from __future__ import annotations

import ast
import math
import os
import random
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import MODELS, back_maps, dump, fb_conditions, medium_bounds, p, read_model, translation  # noqa: E402

TOL = 1e-9


def dnf(node):
    """Disjunctive normal form of a cobrapy GPR AST as a set of frozensets of gene ids."""
    if node is None:
        return frozenset()
    if isinstance(node, ast.Name):
        return frozenset([frozenset([node.id])])
    if isinstance(node, ast.BoolOp):
        parts = [dnf(v) for v in node.values]
        if isinstance(node.op, ast.Or):
            out = set()
            for s in parts:
                out |= s
            return frozenset(out)
        out = {frozenset()}
        for s in parts:
            out = {a | b for a in out for b in s}
        return frozenset(out)
    if isinstance(node, ast.Expression):
        return dnf(node.body)
    raise TypeError(type(node))


def rule_dnf(r):
    body = getattr(r.gpr, "body", None)
    d = dnf(body)
    # minimise (drop supersets) so that equivalent rules compare equal
    return frozenset(s for s in d if not any(o < s for o in d))


def sbml_boundary_species(path):
    with open(p(path)) as fh:
        text = fh.read()
    return {m.group(1) for m in re.finditer(r'<species\s+id="([^"]+)"[^>]*boundaryCondition="true"', text)}


def notes_rules(path):
    """reaction id -> GENE_ASSOCIATION text from SBML notes (Level 2 KBase export)."""
    with open(p(path)) as fh:
        text = fh.read()
    out = {}
    for m in re.finditer(r'<reaction\s+id="([^"]+)"(.*?)</reaction>', text, re.S):
        g = re.search(r"GENE_ASSOCIATION:([^<]*)<", m.group(2))
        out[m.group(1)] = g.group(1).strip() if g else None
    return out


def compare_structure(label):
    cfg = MODELS[label]
    orig, view = read_model(cfg["orig"]), read_model(cfg["view"])
    met_back, rxn_back = back_maps(label)
    tr = translation(label)
    bspecies = sbml_boundary_species(cfg["orig"])
    repairs = {r["reaction"]: r for r in tr["gpr_repairs"]}
    res = {"label": label, "n_reactions": [len(orig.reactions), len(view.reactions)],
           "n_metabolites": [len(orig.metabolites), len(view.metabolites)], "n_genes": [len(orig.genes), len(view.genes)],
           "sbml_boundary_species": len(bspecies)}
    mism = {"missing_in_published": [], "stoichiometry": [], "bounds": [], "gene_rule": [], "objective": [], "collapse_not_clean": []}
    flips, collapsed, used_orig, dropped_species, dropped_rxns = {}, 0, set(), set(), set()
    for vr in view.reactions:
        oid = rxn_back.get(vr.id, vr.id)
        if oid not in orig.reactions:
            mism["missing_in_published"].append(vr.id)
            continue
        orr = orig.reactions.get_by_id(oid)
        used_orig.add(oid)
        vst = {met_back.get(m.id, m.id): c for m, c in vr.metabolites.items()}
        ost = {m.id: c for m, c in orr.metabolites.items()}
        gone = [m for m in ost if m not in vst]
        flip = False
        if gone:
            # must be exactly one SBML boundary species that only this reaction and cobrapy's EX_<b> use
            ok = len(gone) == 1 and gone[0] in bspecies
            if ok:
                b = orig.metabolites.get_by_id(gone[0])
                others = [r for r in b.reactions if r.id != oid]
                ok = (len(others) == 1 and others[0].boundary and len(others[0].metabolites) == 1
                      and others[0].bounds == (-1000.0, 1000.0) and abs(ost[gone[0]]) == 1)
                if ok:
                    dropped_species.add(b.id)
                    dropped_rxns.add(others[0].id)
            if not ok:
                mism["collapse_not_clean"].append({"view": vr.id, "published": oid, "gone": gone})
                continue
            collapsed += 1
            rest = {k: v for k, v in ost.items() if k not in gone}
            if set(rest) == set(vst) and all(abs(vst[k] - rest[k]) <= TOL for k in rest):
                flip = False
            elif set(rest) == set(vst) and all(abs(vst[k] + rest[k]) <= TOL for k in rest):
                flip = True
            else:
                mism["stoichiometry"].append({"view": vr.id, "published": oid, "view_st": vst, "pub_st": ost})
                continue
            flips[vr.id] = flip
        else:
            if set(ost) != set(vst) or any(abs(vst[k] - ost[k]) > TOL for k in ost):
                mism["stoichiometry"].append({"view": vr.id, "published": oid, "view_st": vst, "pub_st": ost})
        exp_b = (-orr.upper_bound, -orr.lower_bound) if flip else (orr.lower_bound, orr.upper_bound)
        if any(abs(a - b) > TOL for a, b in zip(vr.bounds, exp_b)):
            mism["bounds"].append({"view": vr.id, "published": oid, "view": vr.bounds, "published_bounds": orr.bounds, "flip": flip})
        if rule_dnf(vr) != rule_dnf(orr):
            mism["gene_rule"].append({"view": vr.id, "published": oid, "view_rule": vr.gene_reaction_rule,
                                      "published_rule_as_read": orr.gene_reaction_rule,
                                      "documented_repair": oid in repairs})
        oc = orr.objective_coefficient * (-1 if flip else 1)
        if abs(vr.objective_coefficient - oc) > TOL:
            mism["objective"].append({"view": vr.id, "published": oid, "coef": [vr.objective_coefficient, orr.objective_coefficient]})
    unaccounted_r = sorted(set(r.id for r in orig.reactions) - used_orig - dropped_rxns)
    view_mets_back = {met_back.get(m.id, m.id) for m in view.metabolites}
    unaccounted_m = sorted(set(m.id for m in orig.metabolites) - view_mets_back - dropped_species)
    extra_m = sorted(view_mets_back - set(m.id for m in orig.metabolites))
    res.update({"collapsed_two_step_exchanges_verified": collapsed, "collapsed_flipped": sum(flips.values()),
                "published_reactions_unaccounted": unaccounted_r, "published_metabolites_unaccounted": unaccounted_m,
                "view_metabolites_not_in_published": extra_m,
                "genes_identical": sorted(g.id for g in orig.genes) == sorted(g.id for g in view.genes),
                "objective_direction": [orig.objective_direction, view.objective_direction],
                "mismatches": {k: v for k, v in mism.items()}, "mismatch_counts": {k: len(v) for k, v in mism.items()}})
    uncollapsed_b = sorted(bspecies - dropped_species)
    res["boundary_species_not_collapsed"] = [{"id": b, "reactions": [f"{r.id}: {r.reaction} {r.bounds}" for r in orig.metabolites.get_by_id(b).reactions]}
                                             for b in uncollapsed_b if b in orig.metabolites]
    return orig, view, flips, res


def to_published(view, orig, constraints, rxn_back, flips):
    """Apply {view reaction id: (lb, ub)} to the published model (with orientation flips)."""
    for vid, (lb, ub) in constraints.items():
        r = orig.reactions.get_by_id(rxn_back.get(vid, vid))
        r.bounds = (-ub, -lb) if flips.get(vid) else (lb, ub)


def grow(m):
    v = m.slim_optimize(error_value=float("nan"))
    return float(v) if v is not None else float("nan")


def growth_tests(label, orig, view, flips, seed=1):
    cfg = MODELS[label]
    _, rxn_back = back_maps(label)
    vb = [r for r in view.reactions if r.boundary]
    ex_ids = {r.id for r in view.exchanges}          # what the protocol treats as exchanges (closed, then medium)
    out = []

    def run(name, cons):
        with view, orig:
            for vid, b in cons.items():
                view.reactions.get_by_id(vid).bounds = b
            to_published(view, orig, cons, rxn_back, flips)
            a, b = grow(orig), grow(view)
        out.append({"test": name, "published": a, "view": b, "equal": (math.isnan(a) and math.isnan(b)) or abs(a - b) <= 1e-6 * max(1.0, abs(a))})

    run("published bounds", {})
    run("every boundary reaction open (-1000, 1000)", {r.id: (-1000.0, 1000.0) for r in vb})
    med = medium_bounds(cfg["medium"])
    for c in fb_conditions(cfg["org"]):
        if not c["bigg_ids"]:
            continue
        cons = {rid: (0.0, 1000.0) for rid in ex_ids}
        for rid, lb in med.items():
            if rid in ex_ids:
                cons[rid] = (lb, 1000.0)
        for b in c["bigg_ids"]:
            rid = f"EX_{b}_e"
            if rid in ex_ids:
                cons[rid] = (-10.0, 1000.0)
        run(f"protocol medium, no completion: {c['name']}", cons)
    rng = random.Random(seed)
    for k in range(30):
        cons = {r.id: ((0.0, 1000.0) if rng.random() < 0.3 else (-rng.choice([1.0, 10.0, 1000.0]), 1000.0)) for r in vb}
        run(f"random medium {k}", cons)
    # single-gene deletions in one protocol condition
    first = [c for c in fb_conditions(cfg["org"]) if c["bigg_ids"]][0]
    cons = {rid: (0.0, 1000.0) for rid in ex_ids}
    for rid, lb in med.items():
        if rid in ex_ids:
            cons[rid] = (lb, 1000.0)
    for b in first["bigg_ids"]:
        if f"EX_{b}_e" in ex_ids:
            cons[f"EX_{b}_e"] = (-10.0, 1000.0)
    genes = sorted(g.id for g in view.genes)
    rng.shuffle(genes)
    dels = []
    with view, orig:
        for vid, b in cons.items():
            view.reactions.get_by_id(vid).bounds = b
        to_published(view, orig, cons, rxn_back, flips)
        for gid in genes[:60]:
            with view, orig:
                view.genes.get_by_id(gid).knock_out()
                orig.genes.get_by_id(gid).knock_out()
                a, b = grow(orig), grow(view)
            dels.append({"gene": gid, "published": a, "view": b, "equal": (math.isnan(a) and math.isnan(b)) or abs(a - b) <= 1e-6 * max(1.0, abs(a))})
    return out, {"condition": first["name"], "deletions": dels}


def notes_vs_cobra(label, orig):
    """Published gene-association text (SBML notes) against the rules cobrapy built (iGD1575 only)."""
    rules = notes_rules(MODELS[label]["orig"])
    diffs = []
    for rid, text in rules.items():
        if text is None or rid not in orig.reactions:
            continue
        genes_text = set(re.findall(r"[A-Za-z0-9_.\-]+", text)) - {"and", "or", "AND", "OR"}
        genes_cobra = {g.id for g in orig.reactions.get_by_id(rid).genes}
        if genes_text != genes_cobra:
            diffs.append({"reaction": rid, "notes": text, "cobra_rule": orig.reactions.get_by_id(rid).gene_reaction_rule})
    gene_ids = sorted(g.id for g in orig.genes)
    odd = [g for g in gene_ids if not re.match(r"^sm[abc]_?\d+$", g, re.I)]
    orphan = [g.id for g in orig.genes if len(g.reactions) == 0]
    return {"reactions_with_notes": sum(1 for v in rules.values() if v is not None), "gene_set_differences": diffs,
            "n_genes_read": len(gene_ids), "genes_not_shaped_like_sm_locus_tags": odd, "genes_in_no_reaction": orphan}


def main():
    report = {}
    for label in MODELS:
        orig, view, flips, res = compare_structure(label)
        tests, dels = growth_tests(label, orig, view, flips)
        res["growth_tests"] = tests
        res["growth_tests_all_equal"] = all(t["equal"] for t in tests)
        res["growth_tests_summary"] = {"n": len(tests), "n_equal": sum(t["equal"] for t in tests),
                                       "n_positive_growth": sum(1 for t in tests if t["view"] > 1e-6)}
        res["single_gene_deletions"] = dels
        res["single_gene_deletions_all_equal"] = all(d["equal"] for d in dels["deletions"])
        if label == "iGD1575":
            res["published_notes_vs_cobrapy_rules"] = notes_vs_cobra(label, orig)
        report[label] = res
        print(label, {k: res[k] for k in ("n_reactions", "n_metabolites", "n_genes", "collapsed_two_step_exchanges_verified",
                                         "collapsed_flipped", "mismatch_counts", "genes_identical",
                                         "published_reactions_unaccounted", "published_metabolites_unaccounted",
                                         "growth_tests_summary", "single_gene_deletions_all_equal")}, flush=True)
    dump(report, "check_translation.json")


if __name__ == "__main__":
    main()
