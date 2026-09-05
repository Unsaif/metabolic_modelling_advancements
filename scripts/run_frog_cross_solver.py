"""Cross-solver FROG-style reproducibility pass over a set of models.

For each model: load as shipped, run FROG under each solver, write fixtures, and compare solvers pairwise.
Large models get FVA and reaction deletion on a seeded random subset (documented in the fixture).
"""
from __future__ import annotations

import argparse
import glob
import gzip
import json
import os
import shutil
import sys
import tempfile
import time
import traceback

import cobra

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gembench import frog as F  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "frog_cross_solver")


def load_any(path: str) -> cobra.Model:
    if path.endswith(".xml.gz"):
        with gzip.open(path, "rb") as fh, tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as tmp:
            shutil.copyfileobj(fh, tmp)
            tmpname = tmp.name
        try:
            return cobra.io.read_sbml_model(tmpname)
        finally:
            os.unlink(tmpname)
    if path.endswith(".xml"):
        return cobra.io.read_sbml_model(path)
    if path.endswith(".mat"):
        return cobra.io.load_matlab_model(path)
    if path.endswith(".json"):
        return cobra.io.load_json_model(path)
    raise ValueError(path)


def model_set() -> list:
    ext = os.path.join(ROOT, "external")
    paths = []
    paths += [os.path.join(ROOT, "models", m) for m in ["iJR904.xml", "iAF1260.xml", "iJO1366.xml", "iML1515.xml", "yeast-GEM.xml"]]
    paths += sorted(glob.glob(os.path.join(ext, "COBRA.models", "xml", "*.xml")))
    paths += sorted(glob.glob(os.path.join(ext, "COBRA.models", "mat", "*.mat")))
    paths += sorted(glob.glob(os.path.join(ROOT, "models", "embl", "*.xml.gz")))
    paths += [os.path.join(ROOT, "models", "Human-GEM.xml")]
    return paths


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--solvers", nargs="+", default=["glpk", "hybrid"])
    ap.add_argument("--max-fva", type=int, default=300, help="FVA on a seeded random subset if more reactions than this")
    ap.add_argument("--max-del", type=int, default=600, help="reaction deletion subset if more reactions than this")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--skip-existing", action="store_true")
    args = ap.parse_args()

    paths = model_set()
    if args.limit:
        paths = paths[: args.limit]
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for path in paths:
        name = os.path.basename(path).replace(".xml.gz", "").replace(".xml", "").replace(".mat", "")
        mdir = os.path.join(OUT, name)
        cmp_path = os.path.join(mdir, "comparison.json")
        if args.skip_existing and os.path.exists(cmp_path):
            rows.append(json.load(open(cmp_path)))
            continue
        t0 = time.time()
        try:
            model = load_any(path)
            model.id = name
            reports = {}
            for s in args.solvers:
                rep = F.frog(model, solver=s, max_fva_reactions=args.max_fva, max_deletion_reactions=args.max_del, processes=1)
                rep.write(os.path.join(mdir, f"cobrapy_{s}"))
                reports[s] = rep
            comp = F.compare(reports[args.solvers[0]], reports[args.solvers[1]])
            comp.update({"file": os.path.relpath(path, ROOT), "n_reactions": len(model.reactions), "n_genes": len(model.genes),
                         "objective_" + args.solvers[0]: reports[args.solvers[0]].objective_value,
                         "objective_" + args.solvers[1]: reports[args.solvers[1]].objective_value,
                         "seconds": time.time() - t0, "error": ""})
        except Exception as e:  # keep going; record the failure
            comp = {"model_id": name, "file": os.path.relpath(path, ROOT), "error": f"{type(e).__name__}: {e}",
                    "traceback": traceback.format_exc()[-800:], "seconds": time.time() - t0}
        os.makedirs(mdir, exist_ok=True)
        with open(cmp_path, "w") as fh:
            json.dump(comp, fh, indent=2, default=str)
        rows.append(comp)
        print(f"{name:45s} rxns={comp.get('n_reactions','?'):>6} obj_agree={comp.get('objective_agree','?')} "
              f"gene_dis={comp.get('gene_deletion_growth_call_disagreements','?')} rxn_dis={comp.get('reaction_deletion_growth_call_disagreements','?')} "
              f"fva>1e-6={comp.get('fva_n_exceeding_1e-6','?')} {comp.get('seconds',0):.0f}s {comp.get('error','')[:80]}", flush=True)
    with open(os.path.join(OUT, "divergence_matrix.json"), "w") as fh:
        json.dump(rows, fh, indent=2, default=str)


if __name__ == "__main__":
    main()
