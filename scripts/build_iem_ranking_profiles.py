"""Build the biomarker profiles and the readout panel for the IEM disease-ranking study.

Usage: python scripts/build_iem_ranking_profiles.py [--model-file external/COBRA.models/mat/Harvey_1_03d.mat]

Profiles (one per simulated IEM; see docs/studies/wbm-iem-ranking-plan.md):
  lab                   the protocol's biomarker tuples (iem_protocol_v0.2.json) with their expected directions
  lab_hpo_corroborated  the lab tuples whose linked lab-table row is corroborated by an HPO annotation
  hpo                   HPO metabolite-concentration annotations of the disorder (via its Orphanet code, as linked
                        in v0.1), mapped to readouts by data/iem/hpo_term_readout_map_v0.1.tsv; annotations with
                        frequency "Excluded (0%)" are dropped
  hpo_frequent          as hpo, restricted to Obligate, Very frequent and Frequent annotations
A readout with both directions in one profile is dropped from that profile. Readouts that HPO profiles need and
the protocol panel lacks are written to data/iem/iem_ranking_extra_readouts_v0.1.txt.

--version v0.2 (plan section "HPO profiles v0.2") builds the HPO sets directly from Orphanet instead:
  links  data/iem/iem_orphanet_links_v0.2.tsv (IEM -> Orphanet disorders, every change with its reason)
  terms  every annotation of a linked disorder lying under HP:0001939 (Abnormality of metabolism/homeostasis),
         HP:0003117 (Abnormal circulating hormone concentration) or HP:0040085 (Abnormal circulating aldosterone
         concentration), mapped or excluded by data/iem/hpo_term_readout_map_v0.2.tsv; an unmapped term is an error
  lab_hpo_corroborated  lab tuples whose readout and direction appear in the disorder's hpo profile
It needs --orphanet en_product4.xml and --hpo-obo hp.obo (kept outside the repository; checksums are recorded) and
writes data/iem/iem_ranking_profiles_v0.2.json, the HPO readouts outside the protocol panel
(iem_ranking_extra_readouts_v0.2.txt; the main run's panel is the protocol's readouts plus this list).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "iem")
PROTOCOL = os.path.join(DATA, "iem_protocol_v0.2.json")
LAB = os.path.join(DATA, "iem_biomarkers_v0.1_linked.tsv")
HPO_ONLY = os.path.join(DATA, "hpo_metabolite_tuples_not_in_lab_table.tsv")
MAP = os.path.join(DATA, "hpo_term_readout_map_v0.1.tsv")
OUT_JSON = os.path.join(DATA, "iem_ranking_profiles_v0.1.json")
OUT_EXTRA = os.path.join(DATA, "iem_ranking_extra_readouts_v0.1.txt")
FREQUENT = {"Obligate (100%)", "Very frequent (99-80%)", "Frequent (79-30%)"}
EXCLUDED = "Excluded (0%)"
TERM = re.compile(r"(HP:\d+)\s+(.*?)\s*\[(.*)\]$")


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def direction(label):
    l = label.lower()
    return "Increased" if "incre" in l else "Decreased" if "decre" in l else "Unchanged"


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def merge(tuples):
    """(readout, direction) list -> sorted unique list without conflicting readouts; also the dropped readouts."""
    dirs = {}
    for r, d in tuples:
        dirs.setdefault(r, set()).add(d)
    kept = sorted((r, next(iter(d))) for r, d in dirs.items() if len(d) == 1)
    dropped = sorted(r for r, d in dirs.items() if len(d) > 1)
    return [list(t) for t in kept], dropped


ROOTS_V02 = ("HP:0001939", "HP:0003117", "HP:0040085")


def load_obo(path):
    """is_a parents, names, data-version, and replaced_by for obsolete terms."""
    parents, names, data_version = {}, {}, None
    replaced, obsolete = {}, set()
    cur = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("data-version:"):
                data_version = line.split(":", 1)[1].strip()
            if line.startswith("["):
                cur = "" if line == "[Term]" else None
                continue
            if cur is None:
                continue
            if line.startswith("id: "):
                cur = line[4:]
                parents.setdefault(cur, set())
            elif cur and line.startswith("name: "):
                names[cur] = line[6:]
            elif cur and line.startswith("is_a: "):
                parents[cur].add(line[6:].split(" ")[0])
            elif cur and line.startswith("is_obsolete: true"):
                obsolete.add(cur)
            elif cur and line.startswith("replaced_by: "):
                replaced[cur] = line[13:].split(" ")[0]
    return parents, names, data_version, replaced, obsolete


def under_roots(term, parents, roots=ROOTS_V02):
    seen, stack = set(), [term]
    while stack:
        x = stack.pop()
        if x in roots:
            return True
        for q in parents.get(x, ()):
            if q not in seen:
                seen.add(q)
                stack.append(q)
    return False


def load_orphanet(path):
    import xml.etree.ElementTree as ET
    root = ET.parse(path).getroot()
    out = {}
    for d in root.iter("Disorder"):
        code = d.findtext("OrphaCode")
        if not code:
            continue
        anns = []
        for a in d.iter("HPODisorderAssociation"):
            hid, term, freq = a.findtext("HPO/HPOId"), a.findtext("HPO/HPOTerm"), a.findtext("HPOFrequency/Name")
            if not hid or not term:
                raise ValueError(f"Orphanet {code}: annotation without an HPO id or term")
            anns.append((hid, term, freq or ""))
        out[code] = {"name": d.findtext("Name"), "annotations": anns}
    return out, root.attrib.get("date")


def build_v02(args, protocol, panel, profiles_v01):
    orpha, orpha_date = load_orphanet(args.orphanet)
    parents, names, hpo_version, replaced, obsolete = load_obo(args.hpo_obo)
    term_map = {r["hpo_id"]: r for r in tsv(os.path.join(DATA, "hpo_term_readout_map_v0.2.tsv"))}
    links = {r["iem"]: r for r in tsv(os.path.join(DATA, "iem_orphanet_links_v0.2.tsv"))}
    if set(links) != {p["iem"] for p in protocol}:
        raise ValueError("the links table must list every protocol IEM")
    annotations, tuples = [], {"hpo": [], "hpo_frequent": []}
    for p in protocol:
        codes = [c for c in links[p["iem"]]["v02_orpha_codes"].split(";") if c]
        seen = set()
        for code in codes:
            if code not in orpha:
                raise ValueError(f"{p['iem']}: Orphanet {code} is not in en_product4")
            for hid, term, freq in orpha[code]["annotations"]:
                if hid not in parents:
                    raise ValueError(f"{p['iem']}: {hid} {term} is not in hp.obo")
                original = None
                if hid in obsolete:
                    if hid not in replaced:
                        raise ValueError(f"{p['iem']}: {hid} {term} is obsolete without a replacement")
                    original, hid, term = hid, replaced[hid], names.get(replaced[hid], term)
                if not under_roots(hid, parents) or (hid, freq) in seen:
                    continue
                seen.add((hid, freq))
                m = term_map.get(hid)
                if m is None:
                    raise ValueError(f"{p['iem']}: {hid} {term} has no row in hpo_term_readout_map_v0.2.tsv")
                e = {"iem": p["iem"], "orpha_code": code, "hpo_id": hid, "hpo_term": term, "frequency": freq}
                if original:
                    e["obsolete_hpo_id_replaced"] = original
                if m["decision"] == "exclude":
                    e.update(used=False, reason=m["reason"])
                elif freq == EXCLUDED:
                    e.update(used=False, reason="frequency Excluded (0%)")
                else:
                    e.update(used=True, readout=m["readout"], direction=m["direction"], reason="")
                    tuples["hpo"].append((p["iem"], m["readout"], m["direction"]))
                    if freq in FREQUENT:
                        tuples["hpo_frequent"].append((p["iem"], m["readout"], m["direction"]))
                annotations.append(e)
    profiles = {"lab": profiles_v01["lab"], "hpo": {}, "hpo_frequent": {}, "lab_hpo_corroborated": {}}
    conflicts = {k: {} for k in profiles}
    for key in ("hpo", "hpo_frequent"):
        for p in protocol:
            prof, dropped = merge([(r, d) for i, r, d in tuples[key] if i == p["iem"]])
            if prof or dropped:
                profiles[key][p["iem"]], conflicts[key][p["iem"]] = prof, dropped
    for p in protocol:
        hpo = {tuple(t) for t in profiles["hpo"].get(p["iem"], [])}
        corr = [t for t in profiles["lab"][p["iem"]] if tuple(t) in hpo]
        if corr:
            profiles["lab_hpo_corroborated"][p["iem"]] = corr
    extra = sorted({r for key in ("hpo", "hpo_frequent") for prof in profiles[key].values() for r, _ in prof} - set(panel))
    with open(OUT_EXTRA) as fh:
        extra_v01 = [line.strip() for line in fh if line.strip()]
    supplement = sorted(set(extra) - set(extra_v01))
    provenance = {"protocol_sha256": sha256(PROTOCOL), "lab_table_sha256": sha256(LAB),
                  "hpo_map_sha256": sha256(os.path.join(DATA, "hpo_term_readout_map_v0.2.tsv")),
                  "links_sha256": sha256(os.path.join(DATA, "iem_orphanet_links_v0.2.tsv")),
                  "orphanet_en_product4_sha256": sha256(args.orphanet), "orphanet_date": orpha_date,
                  "hpo_obo_sha256": sha256(args.hpo_obo), "hpo_data_version": hpo_version, "hpo_roots": list(ROOTS_V02)}
    return profiles, conflicts, annotations, extra, supplement, provenance


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model-file", default=os.path.join(ROOT, "external", "COBRA.models", "mat", "Harvey_1_03d.mat"))
    ap.add_argument("--version", choices=["v0.1", "v0.2"], default="v0.1")
    ap.add_argument("--orphanet", help="v0.2: Orphanet en_product4.xml")
    ap.add_argument("--hpo-obo", help="v0.2: HPO hp.obo")
    args = ap.parse_args()
    if args.version == "v0.2":
        return main_v02(args)
    with open(PROTOCOL) as fh:
        protocol = json.load(fh)
    candidates = [{"iem": p["iem"], "call_index": p["call_index"]} for p in protocol]
    names = {p["iem"] for p in protocol}
    panel = list(dict.fromkeys(r for p in protocol for r, _ in p["biomarkers"]))
    lab_rows = tsv(LAB)
    hpo_only = tsv(HPO_ONLY)
    term_map = {r["hpo_id"]: r for r in tsv(MAP)}

    profiles = {"lab": {}, "lab_hpo_corroborated": {}, "hpo": {}, "hpo_frequent": {}}
    conflicts = {k: {} for k in profiles}
    corroborated = {(r["iem_abbr"], r["biomarker_reaction"]) for r in lab_rows if r["hpo_check"] == "corroborated"}
    for p in protocol:
        lab = [(r, direction(label)) for r, label in p["biomarkers"] if direction(label) != "Unchanged"]
        profiles["lab"][p["iem"]], conflicts["lab"][p["iem"]] = merge(lab)
        corr = [(r, d) for r, d in lab if (p["iem"], r) in corroborated]
        profiles["lab_hpo_corroborated"][p["iem"]], conflicts["lab_hpo_corroborated"][p["iem"]] = merge(corr)

    # HPO annotations: evidence terms of lab rows plus the HPO-only tuples, one entry per (IEM, HPO term).
    pairs = {}
    for r in lab_rows:
        if r["iem_abbr"] not in names or r["hpo_check"] not in ("corroborated", "corroborated_other_biofluid", "direction_conflict"):
            continue
        for term in re.split(r";\s*|\|\s*", r["hpo_evidence"]):
            if not term.strip():
                continue
            m = TERM.match(term.strip())
            if not m:
                raise ValueError(f"cannot parse HPO evidence {term!r} ({r['iem_abbr']} {r['biomarker_reaction']})")
            hid, label, freq = m.groups()
            e = pairs.setdefault((r["iem_abbr"], hid), {"iem": r["iem_abbr"], "hpo_id": hid, "hpo_term": label,
                                                         "frequency": freq, "sources": []})
            e["sources"].append(f"lab:{r['biomarker_reaction']}")
    for r in hpo_only:
        if r["iem_abbr"] not in names:
            continue
        e = pairs.setdefault((r["iem_abbr"], r["hpo_id"]), {"iem": r["iem_abbr"], "hpo_id": r["hpo_id"], "hpo_term": r["hpo_term"],
                                                             "frequency": r["frequency"], "sources": []})
        if e["frequency"] != r["frequency"]:
            raise ValueError(f"frequency differs between sources for {r['iem_abbr']} {r['hpo_id']}")
        e["sources"].append("hpo_only")
    hpo_tuples = {k: [] for k in ("hpo", "hpo_frequent")}
    annotations = []
    for (iem, hid), e in sorted(pairs.items()):
        e["sources"] = sorted(set(e["sources"]))
        m = term_map.get(hid)
        if m is None:
            raise ValueError(f"HPO term {hid} {e['hpo_term']} has no row in {os.path.basename(MAP)}")
        if m["decision"] == "exclude":
            e["used"], e["reason"] = False, m["reason"]
        elif e["frequency"] == EXCLUDED:
            e["used"], e["reason"] = False, "frequency Excluded (0%)"
        else:
            e["used"], e["reason"] = True, ""
            e["readout"], e["direction"] = m["readout"], m["direction"]
            hpo_tuples["hpo"].append((iem, m["readout"], m["direction"]))
            if e["frequency"] in FREQUENT:
                hpo_tuples["hpo_frequent"].append((iem, m["readout"], m["direction"]))
        annotations.append(e)
    for key, tuples in hpo_tuples.items():
        for iem in sorted(names):
            prof, dropped = merge([(r, d) for i, r, d in tuples if i == iem])
            if prof or dropped:
                profiles[key][iem], conflicts[key][iem] = prof, dropped

    extra = sorted({r for key in ("hpo", "hpo_frequent") for prof in profiles[key].values() for r, _ in prof} - set(panel))

    missing = None
    if os.path.exists(args.model_file):
        import scipy.io as sio
        d = sio.loadmat(args.model_file, squeeze_me=False, struct_as_record=False)
        m = d[[k for k in d if not k.startswith("__")][0]][0, 0]
        rxns = {str(x[0]) for x in m.rxns[:, 0] if len(x)}
        mets = {str(x[0]) for x in m.mets[:, 0] if len(x)}

        def exists(r):
            return r in rxns or (r.startswith("DM_") and r[3:] in mets)
        missing = {"panel": [r for r in panel if not exists(r)], "extra": [r for r in extra if not exists(r)]}

    def size(k):
        return {"n_profiles": sum(1 for v in profiles[k].values() if v), "n_tuples": sum(len(v) for v in profiles[k].values()),
                "n_conflicts_dropped": sum(len(v) for v in conflicts[k].values())}

    out = {"provenance": {"protocol_sha256": sha256(PROTOCOL), "lab_table_sha256": sha256(LAB),
                          "hpo_only_sha256": sha256(HPO_ONLY), "hpo_map_sha256": sha256(MAP),
                          "model_file": os.path.relpath(args.model_file, ROOT) if missing is not None else None},
           "candidates": candidates, "panel_protocol": panel, "extra_readouts": extra,
           "readouts_missing_from_model": missing,
           "summary": {k: size(k) for k in profiles},
           "profiles": profiles, "conflicting_readouts_dropped": {k: {i: v for i, v in c.items() if v} for k, c in conflicts.items()},
           "hpo_annotations": annotations}
    with open(OUT_JSON, "w") as fh:
        json.dump(out, fh, indent=1)
    with open(OUT_EXTRA, "w") as fh:
        fh.write("".join(r + "\n" for r in extra))
    print(json.dumps({"summary": out["summary"], "n_panel": len(panel), "extra_readouts": extra,
                      "readouts_missing_from_model": missing,
                      "n_hpo_annotations": len(annotations), "n_used": sum(a["used"] for a in annotations)}, indent=1))



def model_check(model_file, readouts):
    import scipy.io as sio
    d = sio.loadmat(model_file, squeeze_me=False, struct_as_record=False)
    m = d[[k for k in d if not k.startswith("__")][0]][0, 0]
    rxns = {str(x[0]) for x in m.rxns[:, 0] if len(x)}
    mets = {str(x[0]) for x in m.mets[:, 0] if len(x)}
    return [r for r in readouts if not (r in rxns or (r.startswith("DM_") and r[3:] in mets))]


def main_v02(args):
    if not (args.orphanet and args.hpo_obo):
        raise SystemExit("--version v0.2 needs --orphanet and --hpo-obo")
    with open(PROTOCOL) as fh:
        protocol = json.load(fh)
    with open(OUT_JSON) as fh:
        v01 = json.load(fh)
    panel = v01["panel_protocol"]
    profiles, conflicts, annotations, extra, supplement, provenance = build_v02(args, protocol, panel, v01["profiles"])
    missing = model_check(args.model_file, extra) if os.path.exists(args.model_file) else None
    if missing:
        raise SystemExit(f"HPO readouts missing from the model: {missing}")
    summary = {k: {"n_profiles": sum(1 for v in profiles[k].values() if v), "n_tuples": sum(len(v) for v in profiles[k].values()),
                   "n_conflicts_dropped": sum(len(v) for v in conflicts[k].values())} for k in profiles}
    out = {"provenance": provenance, "candidates": v01["candidates"], "panel_protocol": panel, "extra_readouts": extra,
           "supplement_readouts": supplement, "summary": summary, "profiles": profiles,
           "conflicting_readouts_dropped": {k: {i: v for i, v in c.items() if v} for k, c in conflicts.items()},
           "hpo_annotations": annotations}
    path = os.path.join(DATA, "iem_ranking_profiles_v0.2.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    with open(os.path.join(DATA, "iem_ranking_extra_readouts_v0.2.txt"), "w") as fh:
        fh.write("".join(r + "\n" for r in extra))
    print(json.dumps({"summary": summary, "n_extra": len(extra), "supplement": supplement,
                      "n_annotations": len(annotations), "n_used": sum(a["used"] for a in annotations)}, indent=1))


if __name__ == "__main__":
    main()
