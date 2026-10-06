"""Check 1c: recompute some of the rebuild's validation numbers from this check's own reading of the supplement.

Medium (Table S8a, read here from the PDF text): ca2, cbl1, cobalt2, fe2, fe3, h2s, k, mg2, mn2, na1, pheme, pi, zn2 at
10; nh4, co2, h2o at 100; every other exchange lower bound 0, upper bounds as published, the 2 sinks and 3 demands as
published. Carbon sources at the Table S8b rates. Published in-silico growth (Table S3): glucose 0.24, arabinan 0.23,
dextran 0.26, levan 0.24, pullulan 0.25, glycogen 0.25, inulin 0.24, rhamnose 0.0012, GlcNAc-alpha-1,4-Core 2 0.14,
amylopectin 0.25.
Rich medium (supplement p. 7: every exchange -1000..1000): number of genes whose deletion stops growth (< 1e-6).
Also: growth on glucose with ATPM at 0, 8.43 and 84.3, and the three rows the rebuild reports as not reproduced, at
the printed and the alternative rates.

Usage: python -I validate_rebuild.py   (writes validate_rebuild.json)
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from common import ROOT  # noqa: E402

import cobra  # noqa: E402

S8A = {"ca2": 10, "cbl1": 10, "cobalt2": 10, "fe2": 10, "fe3": 10, "h2s": 10, "k": 10, "mg2": 10, "mn2": 10, "na1": 10,
       "pheme": 10, "pi": 10, "zn2": 10, "nh4": 100, "co2": 100, "h2o": 100}
TESTS = [("Glucose", "glc__D", 10, 0.24), ("Arabinan", "arabinan101", 0.1179, 0.23), ("Dextran", "dextran40", 0.125, 0.26),
         ("Levan", "levan1000", 0.01, 0.24), ("Pullulan", "pullulan1200", 0.00833, 0.25), ("Glycogen", "glycogen1500", 0.0067, 0.25),
         ("Inulin", "inulin", 0.33333, 0.24), ("Rhamnose (printed rate)", "rmn", 12, 0.0012), ("Rhamnose (rate 10)", "rmn", 10, 0.0012),
         ("GlcNAc-Core 2 (printed rate)", "gncore2", 2.7273, 0.14), ("GlcNAc-Core 2 (60/C = 2)", "gncore2", 2.0, 0.14),
         ("Amylopectin as strch1", "strch1", 0.91, 0.25), ("Amylopectin as starch1200", "starch1200", 0.00833, 0.25)]


def ex_of(m, base):
    """The exchange reaction of metabolite <base>_e (published exchange ids do not always follow the metabolite)."""
    met = m.metabolites.get_by_id(f"{base}_e")
    exs = [r for r in met.reactions if r.boundary]
    assert len(exs) == 1, (base, exs)
    return exs[0]


def setup(m, carbon=None, rate=None):
    for r in m.exchanges:
        r.lower_bound = 0.0
    for k, v in S8A.items():
        ex_of(m, k).lower_bound = -float(v)
    if carbon:
        ex_of(m, carbon).lower_bound = -float(rate)


def main():
    m = cobra.io.read_sbml_model(os.path.join(ROOT, "models/curated/iAH991/iAH991_rebuilt.xml"))
    m.solver = "glpk"
    rep = {"objective": str(m.objective.expression), "atpm_bounds": list(m.reactions.get_by_id("ATPM").bounds)}
    rows = []
    for name, cid, rate, pub in TESTS:
        with m:
            setup(m, cid, rate)
            g = m.slim_optimize(error_value=float("nan"))
        rows.append({"condition": name, "substrate": cid, "rate": rate, "published": pub, "rebuilt": g,
                     "rounds_to_published": round(g, 2 if pub >= 0.01 else 4) == pub})
    rep["table_S3"] = rows
    atpm = {}
    for lb in (0.0, 8.43, 84.3):
        with m:
            setup(m, "glc__D", 10)
            m.reactions.get_by_id("ATPM").lower_bound = lb
            atpm[str(lb)] = m.slim_optimize(error_value=float("nan"))
    rep["glucose_growth_by_ATPM_lower_bound"] = atpm
    with m:
        for r in m.exchanges:
            r.bounds = (-1000.0, 1000.0)
        wt = m.slim_optimize()
        ess = []
        for gene in m.genes:
            with m:
                gene.knock_out()
                v = m.slim_optimize(error_value=0.0)
                if not (v >= 1e-6):
                    ess.append(gene.id)
    rep["rich_medium"] = {"wt_growth": wt, "essential_genes": len(ess), "published": 61}
    json.dump(rep, open(os.path.join(HERE, "validate_rebuild.json"), "w"), indent=1, default=float)
    for r in rows:
        print(f"{r['condition']:34s} rate {r['rate']:<8} published {r['published']:<7} rebuilt {r['rebuilt']:.4f} {'OK' if r['rounds_to_published'] else 'differs'}")
    print("ATPM:", atpm)
    print("rich medium:", rep["rich_medium"])


if __name__ == "__main__":
    main()
