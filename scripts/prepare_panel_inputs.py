"""Prepare an evaluation organism's non-outcome inputs for transfer study v1 (docs/studies/transfer-method-v1.md, 3.2-3.4).

Run only after the method freeze. No fitness table is read.

Subcommands
-----------
media      --orgs A,B   map the base media of the organisms' carbon-source experiments with gembench.feba_media and write
                        data/studies/transfer_v1/media_bigg.tsv (+ media_mapping_report.json); components missing from
                        the dictionary are listed so that they can be added by chemical identity (dictionary section 2)
conditions --orgs A,B   list carbon-source condition names (condition_2 empty or DMSO) that neither the development table
                        nor data/studies/transfer_v1/carbon_sources_bigg.tsv maps yet
configure  --orgs A,B   choose each organism's gap-fill reference condition with the frozen rule and add the organism to
                        data/studies/transfer_v1/organisms_panel_A.json (study media and carbon tables, panel gene map)
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from gembench import feba_media as FM  # noqa: E402
from gembench import fitness_browser as FB  # noqa: E402
from gembench import transfer as T  # noqa: E402

STUDY = os.path.join(ROOT, "data", "studies", "transfer_v1")
PANEL_FB = os.path.join(ROOT, "data", "fitness_browser_panel")
DEV_CARBON = os.path.join(ROOT, "data", "reference", "fitness_browser_carbon_sources_bigg.tsv")
STUDY_CARBON = os.path.join(STUDY, "carbon_sources_bigg.tsv")
STUDY_MEDIA = os.path.join(STUDY, "media_bigg.tsv")
COMPONENTS = [os.path.join(STUDY, "feba_component_bigg.tsv"),            # frozen with the method
              os.path.join(STUDY, "feba_component_bigg_panel_A.tsv")]    # additions after the method freeze
UNIVERSE = os.path.join(ROOT, "external", "carveme", "universe_bacteria.xml.gz")
ALLOWED_C2 = ("", "Dimethyl Sulfoxide")


def experiments(org):
    return pd.read_table(os.path.join(PANEL_FB, org, "experiments.tsv"), dtype=str, keep_default_na=False)


def carbon_source_experiments(e):
    return e[(e["expGroup"] == "carbon source") & (e["condition_2"].isin(ALLOWED_C2))]


def universe_exchanges():
    ux = set()
    with gzip.open(UNIVERSE, "rt") as fh:
        for line in fh:
            if "<reaction " in line and 'id="R_EX_' in line:
                i = line.index('id="R_') + 6
                ux.add(line[i:line.index('"', i)])
    return ux


def cmd_media(orgs):
    exps = {o: experiments(o) for o in orgs}
    names = sorted({m for e in exps.values() for m in carbon_source_experiments(e)["media"]})
    comps = [p for p in COMPONENTS if os.path.exists(p)]
    rows = FM.media_table_rows(names, os.path.join(STUDY, "feba", "media"), os.path.join(STUDY, "feba", "mixes"), comps, exps)
    ok = [r for r in rows if r["status"] == "mapped"]
    pd.DataFrame(ok, columns=["media", "aerobic", "bigg_ids", "note", "trace_components"]).to_csv(STUDY_MEDIA, sep="\t", index=False)
    report = {r["media"]: {k: r[k] for k in ("status", "aerobic", "unmapped", "undefined", "ignored", "bigg_ids", "trace_components")}
              for r in rows}
    json.dump(report, open(os.path.join(STUDY, "media_mapping_report.json"), "w"), indent=1)
    for r in rows:
        print(f"{r['media']}: {r['status']} aerobic={r['aerobic']} unmapped={r['unmapped']} undefined={r['undefined']}")


def carbon_tables():
    dev = pd.read_table(DEV_CARBON, dtype=str, keep_default_na=False)
    if os.path.exists(STUDY_CARBON):
        study = pd.read_table(STUDY_CARBON, dtype=str, keep_default_na=False, comment="#")
    else:
        study = dev.iloc[0:0]
    return dev, study


def cmd_conditions(orgs):
    dev, study = carbon_tables()
    known = set(dev["condition"]) | set(study["condition"])
    for o in orgs:
        cs = carbon_source_experiments(experiments(o))
        new = sorted(set(cs["condition_1"]) - known)
        print(f"{o}: {cs['condition_1'].nunique()} carbon-source conditions; {len(new)} not yet mapped")
        for n in new:
            print(f"  {n}")


def cmd_configure(orgs):
    cfg_path = os.path.join(STUDY, "organisms_panel_A.json")      # organisms.json is frozen with the method
    cfg = json.load(open(cfg_path)) if os.path.exists(cfg_path) else {
        "note": "Panel A organisms of transfer study v1, configured after the method freeze by scripts/prepare_panel_inputs.py",
        "organisms": {}}
    ux = universe_exchanges()
    media_map = pd.read_table(STUDY_MEDIA, dtype=str, keep_default_na=False)
    cs_map = pd.read_table(STUDY_CARBON, dtype=str, keep_default_na=False, comment="#")
    pin = json.load(open(os.path.join(ROOT, "models", "embl_pinned", "MANIFEST.json")))
    for o in orgs:
        conds = T.conditions_from_metadata(experiments(o), cs_map, media_map)
        ref = T.choose_reference_condition(conds, ux)
        entry = {"role": "evaluation_panel_A", "fb_dir": f"data/fitness_browser_panel/{o}",
                 "embl_model": f"models/embl_pinned/{o}.xml.gz", "genpept_map": f"data/genpept_panel/{o}_genpept_map.tsv",
                 "media_table": os.path.relpath(STUDY_MEDIA, ROOT), "carbon_table": os.path.relpath(STUDY_CARBON, ROOT),
                 "reference": ({"medium": ref["media"], "carbon_exchange": f"EX_{ref['bigg_ids'][0]}_e",
                                "condition": ref["name"], "rule": "gembench.transfer.choose_reference_condition"} if ref else None),
                 "n_mapped_conditions_from_metadata": sum(1 for c in conds if c["bigg_ids"]),
                 "n_conditions_from_metadata": len(conds)}
        if o not in pin.get("models", pin):
            entry["warning"] = "pinned model not listed in MANIFEST.json"
        cfg["organisms"][o] = entry
        print(o, json.dumps(entry["reference"]), entry["n_mapped_conditions_from_metadata"], "/", entry["n_conditions_from_metadata"])
    json.dump(cfg, open(cfg_path, "w"), indent=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("cmd", choices=["media", "conditions", "configure"])
    ap.add_argument("--orgs", required=True)
    a = ap.parse_args()
    orgs = a.orgs.split(",")
    {"media": cmd_media, "conditions": cmd_conditions, "configure": cmd_configure}[a.cmd](orgs)


if __name__ == "__main__":
    main()
