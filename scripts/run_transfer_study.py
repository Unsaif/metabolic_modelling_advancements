"""Transfer study v1: fitness-blind correction rules evaluated on development organisms and, once frozen,
on an independent panel of new organisms (docs/studies/transfer-method-v1.md).

Subcommands
-----------
prepare  --org X           minimal gap-fill of the pinned EMBL draft for growth on the organism's reference condition
                           (models/transfer_v1/<org>_base.xml.gz + .json). Uses only the draft, the CarveMe universe,
                           the medium recipe and the reference carbon source; never a fitness value.
run      --org X --arms A,B score each arm (the base model plus the arm's transforms) against the organism's fitness
                           data with the fixed protocol and write results/transfer_v1/<role>/<org>/<arm>/.

Arms are defined in data/studies/transfer_v1/arms.json as ordered lists of transforms from gembench.transfer, plus
optional per-organism gene-rule files produced by the blind adjudication procedure.
"""
from __future__ import annotations

import argparse
import contextlib
import gzip
import hashlib
import json
import logging
import os
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
logging.getLogger("cobra").setLevel(logging.ERROR)

import cobra  # noqa: E402

from gembench import fitness_browser as FB  # noqa: E402
from gembench import transfer as T  # noqa: E402
from gembench.cards import BenchmarkCard, LeakageCard, ModelProvenance, now, sha256_of  # noqa: E402
from gembench.gapfill import apply_gapfill, gapfill  # noqa: E402
from gembench.gene_mapping import build_gene_map  # noqa: E402
from gembench.patches import apply_gpr_patches, load_patch_files  # noqa: E402
from gembench.protocols import carbon_fitness_generic as P  # noqa: E402
from run_carbon_fitness_generic import per_condition, score  # noqa: E402

STUDY = os.path.join(ROOT, "data", "studies", "transfer_v1")
UNIVERSE = os.path.join(ROOT, "external", "carveme", "universe_bacteria.xml.gz")
MODELS = os.path.join(ROOT, "models", "transfer_v1")
RESULTS = os.path.join(ROOT, "results", "transfer_v1")
FIXED = dict(carbon_uptake=-10.0, growth_threshold=1e-3, fitness_threshold=-2.0, drop_rich_medium_essentials=False,
             complete_medium_transport=True, solver="glpk")


def study_leakage(role):
    if role == "development":
        used = ("yes: development organism; its phenotypes guided the correction rules under test (development cycles 1-7 "
                "and this study's rule selection); scores are retrospective")
    else:
        used = ("no for the correction rules: frozen before this organism's fitness tables were downloaded (see the study "
                "freeze manifests); the base model's gap-fill uses observed wild-type growth on one reference condition")
    return LeakageCard(
        ground_truth_used_in_model_curation=used,
        ground_truth_public_since="Fitness Browser releases include Price et al. 2018; exact dates differ by experiment",
        frontier_model_training_exposure="unknown; public accessibility does not establish inclusion in a particular model's training data",
        held_out_recommendation="see docs/studies/transfer-method-v1.md and the panel selection record",
        notes=["Fitness < -2 is an operational importance threshold in a pooled mutant assay, not proof of lethality.",
               "Compare arms on matched genes, conditions and finite observations; report wild-type growth coverage separately."])


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


CONFIG_FILES = ("organisms.json", "organisms_panel_A.json")   # development (frozen with the method); panel A (inputs freeze)


def load_config(org):
    found = None
    for name in CONFIG_FILES:
        path = os.path.join(STUDY, name)
        if os.path.exists(path):
            cfg = json.load(open(path))
            if org in cfg["organisms"]:
                if found is not None:
                    raise SystemExit(f"{org} is configured twice")
                found = cfg["organisms"][org]
    if found is None:
        raise SystemExit(f"{org} is not configured in {STUDY}/{' or '.join(CONFIG_FILES)}")
    c = dict(found)
    c["org"] = org
    return c


@contextlib.contextmanager
def reference_tables(cfg):
    """Use the study's media and carbon-source tables (the development defaults unless the organism config names others)."""
    media = os.path.join(ROOT, cfg.get("media_table", "data/reference/fitness_browser_media_bigg.tsv"))
    carbon = os.path.join(ROOT, cfg.get("carbon_table", "data/reference/fitness_browser_carbon_sources_bigg.tsv"))
    orig_m, orig_c = FB.load_media_map, FB.load_carbon_source_map
    FB.load_media_map = lambda path=media: orig_m(path)
    FB.load_carbon_source_map = lambda path=carbon: orig_c(path)
    try:
        yield media, carbon
    finally:
        FB.load_media_map, FB.load_carbon_source_map = orig_m, orig_c


def read_sbml(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            return cobra.io.read_sbml_model(fh)
    return cobra.io.read_sbml_model(path)


def write_sbml(model, path):
    with gzip.open(path, "wt") as fh:
        cobra.io.write_sbml_model(model, fh)


def genes_table(cfg):
    return pd.read_table(os.path.join(ROOT, cfg["fb_dir"], "genes.tsv"), dtype=str, keep_default_na=False)


# ---------------------------------------------------------------------------------------------------------------------

def prepare(org):
    cfg = load_config(org)
    os.makedirs(MODELS, exist_ok=True)
    out_xml = os.path.join(MODELS, f"{org}_base.xml.gz")
    if os.path.exists(out_xml):
        raise SystemExit(f"{out_xml} exists; base models are never overwritten")
    draft_path = os.path.join(ROOT, cfg["embl_model"])
    draft = read_sbml(draft_path)
    draft.solver = "glpk"
    with gzip.open(UNIVERSE, "rt") as fh:
        universe = cobra.io.read_sbml_model(fh)
    ref = cfg["reference"]
    with reference_tables(cfg):
        med = FB.base_medium(ref["medium"])
        t0 = time.time()
        res = gapfill(draft, universe, med, extra_uptakes={ref["carbon_exchange"]: -10.0}, min_growth=0.05, time_limit_s=900)
    print(f"{org}: gap-fill {res.status}; added {res.added_reactions}; growth {res.growth_before:.4f} -> {res.growth_after:.4f} "
          f"({time.time()-t0:.0f}s); rejected energy-cycle sets {res.rejected_energy_cycles}", flush=True)
    if not res.status.lower().startswith("optimal") and not res.status == "not_needed":
        raise SystemExit(f"{org}: gap-fill failed: {res.status}")
    base = apply_gapfill(draft, universe, res.added_reactions, note=f"transfer_v1 base: minimal gap-fill on {ref['medium']} + {ref['carbon_exchange']}")
    base.id = f"{org}_transfer_v1_base"
    write_sbml(base, out_xml)
    info = {"org": org, "role": cfg["role"], "draft": cfg["embl_model"], "draft_sha256": sha256_file(draft_path),
            "universe": os.path.relpath(UNIVERSE, ROOT), "universe_sha256": sha256_file(UNIVERSE),
            "reference": ref, "min_growth": 0.05, "status": res.status,
            "added_reactions": [{"id": a, "reaction": universe.reactions.get_by_id(a).reaction,
                                 "name": universe.reactions.get_by_id(a).name} for a in res.added_reactions],
            "rejected_energy_cycle_sets": res.rejected_energy_cycles,
            "growth_before": res.growth_before, "growth_after": res.growth_after,
            "base_sha256": sha256_file(out_xml), "created": now()}
    json.dump(info, open(os.path.join(MODELS, f"{org}_base.json"), "w"), indent=1)


# ---------------------------------------------------------------------------------------------------------------------

def apply_arm(model, cfg, arm, gm_builder):
    """Apply an arm's transforms in order; returns (model, applied records, gene map)."""
    genes = genes_table(cfg)
    universe_patches = json.load(open(os.path.join(ROOT, "data/reference/universe_patches_v0.1.json")))
    model_patches = json.load(open(os.path.join(ROOT, "data/reference/model_patches_v0.1.json")))
    applied = []
    gm = gm_builder(model)
    for step in arm["transforms"]:
        if step == "universal_reaction_patches":
            rec = T.universal_reaction_patches(model, universe_patches)
        elif step == "universal_model_additions":
            rec = T.universal_model_additions(model, model_patches, gm)
        elif step == "add_atp_synthase_if_annotated":
            rec = T.add_atp_synthase_if_annotated(model, genes, gm)
        elif step == "normalize_conjunction_rules":
            rec = T.normalize_conjunction_rules(model)
        elif step == "remove_menaquinol_if_no_pathway":
            rec = T.remove_menaquinol_if_no_pathway(model, genes)
        elif step.startswith("decisions:"):
            # "decisions:<path>" applies every decision; "decisions:<path>#R6" or "#R1,R2,R2+R1" only those kinds
            spec, _, kinds = step.split(":", 1)[1].partition("#")
            path = os.path.join(ROOT, spec.format(org=cfg["org"]))
            decisions = json.load(open(path))
            if kinds:
                keep = set(kinds.split(","))
                decisions = dict(decisions, decisions=[d for d in decisions["decisions"] if d.get("decision") in keep])
            rec = T.apply_decisions(model, decisions, gm, set(genes["sysName"]))
        elif step.startswith("gene_rules:"):
            path = step.split(":", 1)[1].format(org=cfg["org"])
            gpr = load_patch_files([os.path.join(ROOT, path)])
            rec = apply_gpr_patches(model, cfg["org"], gm, gpr, statuses=("accepted",), verbose=True)
        else:
            raise ValueError(f"unknown transform {step}")
        applied.append({"transform": step, "records": rec})
        gm = gm_builder(model)          # genes added under Browser locus tags map to themselves
    return model, applied, gm


def run(org, arm_names, out_root, processes):
    cfg = load_config(org)
    arms = json.load(open(os.path.join(STUDY, "arms.json")))["arms"]
    base_path = os.path.join(MODELS, f"{org}_base.xml.gz")
    base_info = json.load(open(os.path.join(MODELS, f"{org}_base.json")))
    fb = FB.load_organism(org, data_dir=os.path.dirname(os.path.join(ROOT, cfg["fb_dir"])))
    sysnames = set(fb.genes["sysName"])
    genpept = os.path.join(ROOT, cfg["genpept_map"])

    def gm_builder(model):
        return build_gene_map(org, [g.id for g in model.genes], sysnames, genpept_path=genpept)

    for name in arm_names:
        arm = arms[name]
        outdir = os.path.join(out_root, cfg["role"], org, name)
        if os.path.exists(os.path.join(outdir, "card.json")):
            print(f"{org}/{name}: exists, skipped", flush=True)
            continue
        os.makedirs(outdir, exist_ok=True)
        t0 = time.time()
        model = read_sbml(base_path)
        model.solver = "glpk"
        with reference_tables(cfg) as (media_table, carbon_table):
            model, applied, gm = apply_arm(model, cfg, arm, gm_builder)
            gm.stats["mapped_with_fitness_data"] = sum(1 for v in gm.model_to_browser.values() if v in set(fb.fitness.index))
            conds = FB.carbon_source_conditions(fb)
            params = P.GenericParams(processes=processes, **FIXED)
            res = P.run(model, fb, conds, gm, params, verbose=False)
        grows = res.wt_growth >= params.growth_threshold
        results = {"condition_level": {"n_conditions_mapped": int(len(res.conditions)), "n_conditions_wt_grows": int(grows.sum()),
                                       "conditions_with_absent_exchange": int(sum(1 for c in res.conditions if res.missing_carbon_exchanges.get(c.key)))},
                   "gene_level_conditions_where_wt_grows": score(res, grows),
                   "gene_map": gm.stats, "counts": res.counts, "timings_s": res.timings_s}
        card = BenchmarkCard(
            benchmark=f"carbon-source fitness benchmark (Fitness Browser RB-TnSeq), organism {org}",
            created=now(), dataset_provenance=fb.provenance,
            model=ModelProvenance(model_id=model.id, file=os.path.relpath(base_path, ROOT), source=f"transfer_v1 base ({base_info['draft']})",
                                  n_reactions=len(model.reactions), n_metabolites=len(model.metabolites), n_genes=len(model.genes),
                                  sha256=sha256_of(base_path)),
            protocol={"study": "transfer_v1", "arm": name, "arm_definition": arm, "applied": applied, "params": res.params.__dict__,
                      "media_mapping": os.path.relpath(media_table, ROOT), "carbon_source_mapping": os.path.relpath(carbon_table, ROOT),
                      "gene_mapping": gm.provenance, "role": cfg["role"], "base_model": base_info},
            leakage=study_leakage(cfg["role"]), results=results, warnings=[])
        card.write(os.path.join(outdir, "card.json"), os.path.join(outdir, "card.md"))
        np.savez_compressed(os.path.join(outdir, "matrices.npz"), sim_growth=res.sim_growth, wt_growth=res.wt_growth,
                            fitness=res.fitness, model_genes=np.array(res.model_genes), browser_genes=np.array(res.browser_genes),
                            conditions=np.array([c.key for c in res.conditions]))
        res.condition_table().to_csv(os.path.join(outdir, "conditions.tsv"), sep="\t", index=False)
        per_condition(res).to_csv(os.path.join(outdir, "per_condition_metrics.tsv"), sep="\t", index=False)
        write_sbml(model, os.path.join(outdir, "model.xml.gz"))
        g = results["gene_level_conditions_where_wt_grows"]
        print(f"{org}/{name}: WT grows {int(grows.sum())}/{len(grows)}; MCC={g.get('mcc', {}).get('point', float('nan')):.4f} "
              f"pairs={g.get('n_gene_condition_pairs')} ({time.time()-t0:.0f}s)", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("prepare")
    p1.add_argument("--org", required=True)
    p2 = sub.add_parser("run")
    p2.add_argument("--org", required=True)
    p2.add_argument("--arms", required=True)
    p2.add_argument("--output-dir", default=RESULTS)
    p2.add_argument("--processes", type=int, default=2)
    a = ap.parse_args()
    if a.cmd == "prepare":
        prepare(a.org)
    else:
        run(a.org, a.arms.split(","), a.output_dir, a.processes)


if __name__ == "__main__":
    main()
