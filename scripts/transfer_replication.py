"""Transfer study v1, replication on panel B (docs/studies/transfer-v1-replication-plan.md).

Runs the frozen v1 method unchanged on the reserved panel B organisms. This wrapper only points the frozen code at the
replication's own configuration files; every rule (media and carbon-source mapping, reference condition, gap-fill,
transforms, adjudication packets, scoring) is the frozen v1 implementation.

Subcommands
-----------
media       --orgs ...   FEBA recipes -> data/studies/transfer_v1_replication/media_bigg_panel_B.tsv (frozen rule; dictionary
                         sections 1-3; section 3 = feba_component_bigg_panel_B.tsv, added after the replication freeze)
conditions  --orgs ...   carbon-source names not yet mapped by carbon_sources_bigg_panel_B.tsv
configure   --orgs ...   reference condition (frozen rule) -> organisms_panel_B.json (role evaluation_panel_B)
prepare     --org X      minimal gap-fill of the pinned draft (frozen run_transfer_study.prepare)
candidates  --org X --out DIR    blind adjudication packet from the U' (UNQ) model (frozen transfer_candidates)
filter      --orgs ...   decisions without R6 assignments to inorganic-ion transport reactions (secondary arm M_noIonR6)
run         --org X --arms A,B   score arms (frozen run_transfer_study.run; arms from the replication's arms.json)
pool        --inputs a.json,b.json --out pooled.json   concatenate per-organism comparison rows of two panels
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
logging.getLogger("cobra").setLevel(logging.ERROR)

import run_transfer_study as RTS  # noqa: E402
from gembench import feba_media as FM  # noqa: E402
from gembench import transfer as T  # noqa: E402

REP = os.path.join(ROOT, "data", "studies", "transfer_v1_replication")
V1 = os.path.join(ROOT, "data", "studies", "transfer_v1")
PANEL_FB = os.path.join(ROOT, "data", "fitness_browser_panel")
MEDIA_TABLE = os.path.join(REP, "media_bigg_panel_B.tsv")
CARBON_TABLE = os.path.join(REP, "carbon_sources_bigg_panel_B.tsv")
CONFIG = "organisms_panel_B.json"
COMPONENTS = [os.path.join(V1, "feba_component_bigg.tsv"),                 # section 1 (frozen with the v1 method)
              os.path.join(V1, "feba_component_bigg_panel_A.tsv"),         # section 2 (panel A inputs freeze)
              os.path.join(REP, "feba_component_bigg_panel_B.tsv")]        # section 3 (after the replication freeze)
ALLOWED_C2 = ("", "Dimethyl Sulfoxide")

# Secondary hypothesis H3 (declared in the replication plan): R6 assignments to transport reactions whose transported
# species are all inorganic ions are dropped. A species is transported if its metabolite occurs in two or more
# compartments of the reaction; protons and water are ignored.
INORGANIC_IONS = frozenset({"zn2", "cu2", "cu", "cobalt2", "mn2", "fe2", "fe3", "ni2", "mobd", "mg2", "ca2", "k", "na1",
                            "cl", "so4", "so3", "pi", "ppi", "nh4", "no3", "no2", "tsul", "slnt", "sel", "tungs", "cd2",
                            "hg2", "pb", "cro4", "aso3", "aso4", "hco3", "co2", "h2s", "n2", "o2"})
IGNORED_CARRIED = frozenset({"h", "h2o"})


def point_to_replication():
    """Make the frozen runner read the replication's arms.json and panel B configuration."""
    RTS.STUDY = REP
    RTS.CONFIG_FILES = (CONFIG,)


def experiments(org):
    return pd.read_table(os.path.join(PANEL_FB, org, "experiments.tsv"), dtype=str, keep_default_na=False)


def cs_experiments(e):
    return e[(e["expGroup"] == "carbon source") & (e["condition_2"].isin(ALLOWED_C2))]


# --------------------------------------------------------------------------------------------------------------------

def cmd_media(orgs):
    exps = {o: experiments(o) for o in orgs}
    names = sorted({m for e in exps.values() for m in cs_experiments(e)["media"]})
    comps = [p for p in COMPONENTS if os.path.exists(p)]
    rows = FM.media_table_rows(names, os.path.join(V1, "feba", "media"), os.path.join(V1, "feba", "mixes"), comps, exps)
    ok = [r for r in rows if r["status"] == "mapped"]
    pd.DataFrame(ok, columns=["media", "aerobic", "bigg_ids", "note", "trace_components"]).to_csv(MEDIA_TABLE, sep="\t", index=False)
    json.dump({r["media"]: {k: r[k] for k in ("status", "aerobic", "unmapped", "undefined", "ignored", "bigg_ids", "trace_components")}
               for r in rows}, open(os.path.join(REP, "media_mapping_report_panel_B.json"), "w"), indent=1)
    for r in rows:
        print(f"{r['media']}: {r['status']} aerobic={r['aerobic']} unmapped={r['unmapped']} undefined={r['undefined']}")


def carbon_table():
    if os.path.exists(CARBON_TABLE):
        return pd.read_table(CARBON_TABLE, dtype=str, keep_default_na=False)
    return pd.read_table(os.path.join(V1, "carbon_sources_bigg.tsv"), dtype=str, keep_default_na=False)


def cmd_conditions(orgs):
    known = set(carbon_table()["condition"])
    for o in orgs:
        cs = cs_experiments(experiments(o))
        new = sorted(set(cs["condition_1"]) - known)
        print(f"{o}: {cs['condition_1'].nunique()} carbon-source conditions; {len(new)} not yet mapped")
        for n in new:
            print(f"  {n}")


def cmd_configure(orgs):
    import gzip
    ux = set()
    with gzip.open(RTS.UNIVERSE, "rt") as fh:
        for line in fh:
            if "<reaction " in line and 'id="R_EX_' in line:
                i = line.index('id="R_') + 6
                ux.add(line[i:line.index('"', i)])
    media_map = pd.read_table(MEDIA_TABLE, dtype=str, keep_default_na=False)
    cs_map = pd.read_table(CARBON_TABLE, dtype=str, keep_default_na=False)
    path = os.path.join(REP, CONFIG)
    cfg = json.load(open(path)) if os.path.exists(path) else {
        "note": "Panel B organisms of the transfer v1 replication, configured after the replication freeze", "organisms": {}}
    for o in orgs:
        conds = T.conditions_from_metadata(experiments(o), cs_map, media_map)
        ref = T.choose_reference_condition(conds, ux)
        cfg["organisms"][o] = {
            "role": "evaluation_panel_B", "fb_dir": f"data/fitness_browser_panel/{o}", "embl_model": f"models/embl_pinned/{o}.xml.gz",
            "genpept_map": f"data/genpept_panel/{o}_genpept_map.tsv",
            "media_table": os.path.relpath(MEDIA_TABLE, ROOT), "carbon_table": os.path.relpath(CARBON_TABLE, ROOT),
            "reference": ({"medium": ref["media"], "carbon_exchange": f"EX_{ref['bigg_ids'][0]}_e", "condition": ref["name"],
                           "rule": "gembench.transfer.choose_reference_condition"} if ref else None),
            "n_mapped_conditions_from_metadata": sum(1 for c in conds if c["bigg_ids"]), "n_conditions_from_metadata": len(conds)}
        print(o, json.dumps(cfg["organisms"][o]["reference"]), cfg["organisms"][o]["n_mapped_conditions_from_metadata"], "/", len(conds))
    json.dump(cfg, open(path, "w"), indent=1)


# --------------------------------------------------------------------------------------------------------------------

def transported_species(reaction) -> set:
    comps = {}
    for met in reaction.metabolites:
        base = met.id.rsplit("_", 1)[0] if "_" in met.id else met.id
        comps.setdefault(base, set()).add(met.compartment or met.id.rsplit("_", 1)[-1])
    return {b for b, cs in comps.items() if len(cs) >= 2} - IGNORED_CARRIED


def is_inorganic_ion_transport(reaction) -> bool:
    carried = transported_species(reaction)
    return bool(carried) and carried <= INORGANIC_IONS


def any_config(org):
    """Configuration of a panel B, panel A or development organism (the filter is also described retrospectively)."""
    for d, name in ((REP, CONFIG), (V1, "organisms.json"), (V1, "organisms_panel_A.json")):
        path = os.path.join(d, name)
        if os.path.exists(path):
            cfg = json.load(open(path))["organisms"]
            if org in cfg:
                return dict(cfg[org], org=org)
    raise SystemExit(f"{org} is not configured")


def cmd_filter(orgs):
    """Write decisions_noIonR6/<org>.json: R6 decisions on inorganic-ion transport reactions become abstentions."""
    from gembench.gene_mapping import build_gene_map
    os.makedirs(os.path.join(REP, "decisions_noIonR6"), exist_ok=True)
    arms = json.load(open(os.path.join(V1, "arms.json")))["arms"]
    report = {}
    for org in orgs:
        cfg = any_config(org)
        genes = RTS.genes_table(cfg)
        model = RTS.read_sbml(os.path.join(RTS.MODELS, f"{org}_base.xml.gz"))
        gpath = os.path.join(ROOT, cfg["genpept_map"])
        with RTS.reference_tables(cfg):
            model, _, _ = RTS.apply_arm(model, cfg, arms["UNQ"],
                                        lambda m: build_gene_map(org, [g.id for g in m.genes], set(genes["sysName"]), genpept_path=gpath))
        dec = json.load(open(os.path.join(V1, "decisions", f"{org}.json")))
        out, dropped = [], []
        for d in dec["decisions"]:
            d = dict(d)
            if d.get("decision") == "R6" and d["reaction"] in model.reactions and \
                    is_inorganic_ion_transport(model.reactions.get_by_id(d["reaction"])):
                dropped.append({"reaction": d["reaction"], "new_rule": d.get("new_rule"),
                                "carried": sorted(transported_species(model.reactions.get_by_id(d["reaction"])))})
                d.update({"decision": "abstain", "new_rule": None,
                          "evidence": "removed by the noIonR6 filter (R6 on an inorganic-ion transport reaction); original: " + str(d.get("evidence", ""))})
            out.append(d)
        json.dump(dict(dec, decisions=out, filter="noIonR6 (docs/studies/transfer-v1-replication-plan.md)"),
                  open(os.path.join(REP, "decisions_noIonR6", f"{org}.json"), "w"), indent=1)
        report[org] = dropped
        print(org, "dropped", [x["reaction"] for x in dropped])
    path = os.path.join(REP, "decisions_noIonR6", "filter_report.json")
    old = json.load(open(path)) if os.path.exists(path) else {}
    old.update(report)
    json.dump(old, open(path, "w"), indent=1)


def cmd_pool(inputs, out):
    rows = []
    for p in inputs.split(","):
        rows += json.load(open(p))
    orgs = [r["org"] for r in rows]
    pairs = {}
    for r in rows:
        pairs.setdefault((r["A"], r["B"]), []).append(r["org"])
    for k, v in pairs.items():
        if len(v) != len(set(v)):
            raise SystemExit(f"duplicate organisms for pair {k}")
    json.dump(rows, open(out, "w"), indent=1)
    print(f"pooled {len(rows)} rows from {len(set(orgs))} organisms -> {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("cmd", choices=["media", "conditions", "configure", "prepare", "candidates", "filter", "run", "pool"])
    ap.add_argument("--orgs")
    ap.add_argument("--org")
    ap.add_argument("--arms")
    ap.add_argument("--out")
    ap.add_argument("--inputs")
    ap.add_argument("--processes", type=int, default=2)
    a = ap.parse_args()
    point_to_replication()
    if a.cmd == "media":
        cmd_media(a.orgs.split(","))
    elif a.cmd == "conditions":
        cmd_conditions(a.orgs.split(","))
    elif a.cmd == "configure":
        cmd_configure(a.orgs.split(","))
    elif a.cmd == "prepare":
        RTS.prepare(a.org)
    elif a.cmd == "candidates":
        import transfer_candidates as TC
        sys.argv = ["transfer_candidates.py", "--org", a.org, "--arm", "UNQ", "--out", a.out, "--processes", str(a.processes)]
        TC.main()
    elif a.cmd == "filter":
        cmd_filter(a.orgs.split(","))
    elif a.cmd == "run":
        RTS.run(a.org, a.arms.split(","), RTS.RESULTS, a.processes)
    elif a.cmd == "pool":
        cmd_pool(a.inputs, a.out)


if __name__ == "__main__":
    main()
