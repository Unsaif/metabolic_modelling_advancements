"""Minimally gap-fill the EMBL draft models for growth on the base medium of their Fitness Browser
carbon-source experiments plus one reference carbon source (glucose; L-lactate for S. oneidensis).

Writes models/gapfilled/<org>.xml.gz and models/gapfilled/<org>_gapfill.json (added reactions, growth
before/after, MILP status). Universe: CarveMe universe_bacteria (external/carveme/universe_bacteria.xml.gz).
"""
from __future__ import annotations

import glob
import gzip
import json
import logging
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
logging.getLogger("cobra").setLevel(logging.ERROR)

import cobra  # noqa: E402

from gembench.cards import now, sha256_of  # noqa: E402
from gembench.fitness_browser import base_medium  # noqa: E402
from gembench.gapfill import apply_gapfill, gapfill  # noqa: E402

EMBL = {"Btheta": "Bacteroides_thetaiotaomicron_VPI_5482", "Putida": "Pseudomonas_putida_KT2440",
        "MR1": "Shewanella_oneidensis_MR_1", "Smeli": "Sinorhizobium_meliloti_1021"}
REFERENCE = {"Btheta": ("Varel_Bryant_medium", "EX_glc__D_e"), "Putida": ("MOPS minimal media_noCarbon", "EX_glc__D_e"),
             "MR1": ("ShewMM_noCarbon", "EX_lac__L_e"), "Smeli": ("RCH2_defined_noCarbon", "EX_glc__D_e")}
OUT = os.path.join(ROOT, "models", "gapfilled")


def main() -> None:
    orgs = sys.argv[1].split(",") if len(sys.argv) > 1 else list(EMBL)
    os.makedirs(OUT, exist_ok=True)
    upath = os.path.join(ROOT, "external", "carveme", "universe_bacteria.xml.gz")
    with gzip.open(upath, "rt") as fh:
        universe = cobra.io.read_sbml_model(fh)
    for org in orgs:
        mpath = glob.glob(os.path.join(ROOT, "models", "embl", f"{EMBL[org]}*.xml.gz"))[0]
        with gzip.open(mpath, "rt") as fh:
            model = cobra.io.read_sbml_model(fh)
        model.solver = "glpk"
        medname, ref_ex = REFERENCE[org]
        med = base_medium(medname)
        res = gapfill(model, universe, med, extra_uptakes={ref_ex: -10.0}, min_growth=0.05, time_limit_s=900)
        print(f"{org}: {res.status}; added {res.added_reactions}; growth {res.growth_before:.4f} -> {res.growth_after:.4f} ({res.seconds:.0f}s)")
        gm = apply_gapfill(model, universe, res.added_reactions,
                           note=f"gembench minimal gap-fill for growth on {medname} + {ref_ex}, {now()}")
        gm.id = f"{model.id}_gapfilled"
        out_xml = os.path.join(OUT, f"{org}.xml.gz")
        with gzip.open(out_xml, "wt") as fh:
            cobra.io.write_sbml_model(gm, fh)
        info = {"org_id": org, "source_model": os.path.relpath(mpath, ROOT), "source_sha256": sha256_of(mpath),
                "universe": os.path.relpath(upath, ROOT), "universe_sha256": sha256_of(upath),
                "medium": medname, "reference_carbon_exchange": ref_ex, "min_growth": 0.05,
                "added_reactions": [{"id": a, "name": universe.reactions.get_by_id(a).name,
                                     "reaction": universe.reactions.get_by_id(a).reaction,
                                     "weight": res.weights_used.get(a)} for a in res.added_reactions],
                "growth_before": res.growth_before, "growth_after": res.growth_after, "milp_status": res.status,
                "n_candidates": res.n_candidates, "seconds": res.seconds, "created": now()}
        with open(os.path.join(OUT, f"{org}_gapfill.json"), "w") as fh:
            json.dump(info, fh, indent=2)


if __name__ == "__main__":
    main()
