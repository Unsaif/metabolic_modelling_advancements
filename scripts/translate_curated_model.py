"""Give a published curated model the BiGG identifiers the transfer study's tables use (Paper 1, post hoc).

The transfer protocol applies media and carbon sources by BiGG exchange identifier (EX_glc__D_e) and checks for
energy-generating cycles with BiGG currency metabolites (atp_c, nadh_c ...). Curated models outside BiGG name
things differently: iSO783 (BioModels MODEL1507180036) uses older BiGG-style identifiers with escaped exchange
names (EX_ac_LPAREN_e_RPAREN_), iGD1575 (diCenzo et al. 2016) uses ModelSEED compounds (cpd00027_e0) and
KBase-style two-step exchanges (cpd00027_e0 <=> cpd00027_b, cpd00027_b <=>).

This script only relabels. Metabolite, compartment and exchange-reaction identifiers change; stoichiometry, bounds,
gene rules and the objective do not. A two-step exchange is collapsed to one (X_e <=>) with the published bounds,
which is the same constraint because the boundary species is unconstrained. The script checks the result: the
optimum must be unchanged with the published bounds and with every exchange open.

Identifier sources, in order: the identifier itself when it is a BiGG identifier (the CarveMe BiGG universe, the
BiGG models in models/, or a BiGG alias in the ModelSEED database); for older BiGG identifiers, the ModelSEED
database's BiGG1 -> compound -> BiGG aliases; for ModelSEED compounds, the compound's BiGG aliases (one candidate,
or one candidate that the BiGG models in models/ use); then explicit decisions in
data/reference/namespace/curated_model_id_overrides.tsv (each with its evidence). Two metabolites that would get
the same identifier keep their own identifiers and are reported.

Usage:
  python -I scripts/translate_curated_model.py --model models/curated/iSO783/MODEL1507180036_url.xml \
      --label iSO783 --scheme bigg --org MR1
  python -I scripts/translate_curated_model.py --model models/curated/iGD1575/ncomms12219-s7.xml \
      --label iGD1575 --scheme modelseed --org Smeli
Writes models/curated/<label>/<label>_bigg_view.xml.gz and <label>_translation.json.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import logging
import os
import re
import sys

import cobra
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
logging.getLogger("cobra").setLevel(logging.ERROR)

NS = os.path.join(ROOT, "data", "reference", "namespace")
UNIVERSE = os.path.join(ROOT, "results", "ppnp_repair_2026_09_06", "evidence", "source", "universe_bacteria.xml.gz")
LOCAL_BIGG = [os.path.join(ROOT, "models", "bigg", "iJN1463.xml"), os.path.join(ROOT, "models", "iML1515", "iML1515.xml")]
OVERRIDES = os.path.join(NS, "curated_model_id_overrides.tsv")
GPR_REPAIRS = os.path.join(NS, "curated_model_gpr_repairs.tsv")
COMPARTMENT = {"c": "c", "c0": "c", "e": "e", "e0": "e", "p": "p", "p0": "p"}
CURRENCY_IDS = ["atp_c", "adp_c", "pi_c", "h_c", "h2o_c", "nadh_c", "nad_c", "nadph_c", "nadp_c", "q8h2_c", "q8_c"]


def species_ids(path):
    """Metabolite base identifiers (no compartment) of a BiGG-namespace SBML file, read without building a model."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        text = fh.read()
    return {m.group(1) for m in re.finditer(r'<species[^>]*\sid="M_([^"]+)_[a-z]{1,2}"', text)}


def load_namespace():
    universe = species_ids(UNIVERSE)
    local = set().union(*[species_ids(p) for p in LOCAL_BIGG if os.path.exists(p)])
    al = pd.read_table(os.path.join(NS, "Unique_ModelSEED_Compound_Aliases.txt"), dtype=str, keep_default_na=False)
    al.columns = ["cpd", "ext", "source"]
    bigg = al[al["source"] == "BiGG"]
    bigg1 = al[al["source"] == "BiGG1"]
    cpd_to_bigg = collections.defaultdict(set)
    for c, e in zip(bigg["cpd"], bigg["ext"]):
        cpd_to_bigg[c].add(e)
    bigg1_to_cpd = collections.defaultdict(set)
    for c, e in zip(bigg1["cpd"], bigg1["ext"]):
        bigg1_to_cpd[e].add(c)
    known = universe | local | set(bigg["ext"])
    return {"universe": universe, "local": local, "known": known, "cpd_to_bigg": cpd_to_bigg, "bigg1_to_cpd": bigg1_to_cpd}


def load_overrides(label):
    if not os.path.exists(OVERRIDES):
        return {}
    t = pd.read_table(OVERRIDES, dtype=str, keep_default_na=False, comment="#")
    t = t[t["model"] == label]
    return dict(zip(t["original_base_id"], t["bigg_id"]))


def repair_gprs(model, label):
    """Gene rules the SBML reader could not parse, restored to the rule the authors wrote (each repair listed with
    its evidence in data/reference/namespace/curated_model_gpr_repairs.tsv)."""
    if not os.path.exists(GPR_REPAIRS):
        return []
    t = pd.read_table(GPR_REPAIRS, dtype=str, keep_default_na=False, comment="#")
    done = []
    for _, row in t[t["model"] == label].iterrows():
        r = model.reactions.get_by_id(row["reaction"])
        before = r.gene_reaction_rule
        if before.strip():
            raise SystemExit(f"{row['reaction']}: expected an unparsed (empty) rule, found {before!r}")
        r.gene_reaction_rule = row["repaired_rule"]
        done.append({"reaction": r.id, "published_text": row["published_text"], "rule_read_by_cobra": before,
                     "repaired_rule": r.gene_reaction_rule, "evidence": row["evidence"]})
    return done


def split_id(met):
    """(base, compartment) of a metabolite identifier such as glc__D_e, cpd00027_e0 or cpd00027_b."""
    m = re.match(r"^(.*)_([a-z]\d?)$", met.id)
    if not m:
        return met.id, met.compartment
    return m.group(1), m.group(2)


def candidate(base, scheme, ns, overrides):
    """(BiGG id or None, how)."""
    if base in overrides:
        return overrides[base], "override"
    if scheme == "bigg":
        if base in ns["known"]:
            return base, "identical"
        cpds = ns["bigg1_to_cpd"].get(base, set()) or ns["bigg1_to_cpd"].get(base.replace("__", "-"), set())
        cands = set().union(*[ns["cpd_to_bigg"][c] for c in cpds]) if cpds else set()
        cands &= ns["known"]
        if len(cands) == 1:
            return next(iter(cands)), "BiGG1 alias"
        return None, ("ambiguous BiGG1 alias: " + ",".join(sorted(cands))) if cands else "no alias"
    if scheme == "modelseed":
        cands = ns["cpd_to_bigg"].get(base, set()) & ns["known"]
        if len(cands) == 1:
            return next(iter(cands)), "ModelSEED BiGG alias"
        if len(cands) > 1:
            in_local = cands & ns["local"]
            if len(in_local) == 1:
                return next(iter(in_local)), "ModelSEED BiGG alias (the one used by the BiGG models in models/)"
            return None, "ambiguous ModelSEED BiGG alias: " + ",".join(sorted(cands))
        return None, "no alias"
    raise ValueError(scheme)


def collapse_two_step_exchanges(model):
    """X_e <=> X_b plus X_b <=> (boundary species exported by KBase/ModelSEED) -> X_e <=> with the published bounds."""
    collapsed = []
    for met in list(model.metabolites):
        if not met.id.endswith("_b"):       # only the boundary species of a KBase/ModelSEED export (cpd00027_b)
            continue
        rxns = list(met.reactions)
        if len(rxns) != 2:
            continue
        bnd = [r for r in rxns if r.boundary and len(r.metabolites) == 1]
        other = [r for r in rxns if r not in bnd]
        if len(bnd) != 1 or len(other) != 1 or len(other[0].metabolites) != 2:
            continue
        ex, tr = bnd[0], other[0]
        partner = [m for m in tr.metabolites if m is not met][0]
        if abs(tr.metabolites[met]) != 1 or abs(tr.metabolites[partner]) != 1 or tr.metabolites[met] * tr.metabolites[partner] > 0:
            continue
        if ex.bounds != (-1000.0, 1000.0) and ex.bounds != (-float("inf"), float("inf")):
            continue                        # only an unconstrained boundary species can be dropped without changing the model
        # orient as partner <=> (partner consumed when flux is positive), keeping the published bounds of `tr`
        sign = -1.0 if tr.metabolites[partner] < 0 else 1.0
        lb, ub = tr.lower_bound, tr.upper_bound
        tr.subtract_metabolites({met: tr.metabolites[met]})
        if sign > 0:                        # tr was "X_b -> partner": flip so positive flux exports partner
            tr.add_metabolites({partner: -2.0 * tr.metabolites[partner]})
            lb, ub = -ub, -lb
        tr.bounds = (lb, ub)
        model.remove_reactions([ex])
        model.remove_metabolites([met])
        collapsed.append({"exchange": tr.id, "dropped_boundary_species": met.id, "dropped_reaction": ex.id})
    return collapsed


def tables_for(org):
    import run_transfer_study as RT
    from gembench import fitness_browser as FB
    cfg = RT.load_config(org)
    with RT.reference_tables(cfg):
        fb = FB.load_organism(org, data_dir=os.path.dirname(os.path.join(ROOT, cfg["fb_dir"])))
        conds = [c for c in FB.carbon_source_conditions(fb) if c.bigg_ids]
        media = sorted({c.media for c in conds})
        med = {m: sorted(FB.base_medium(m).uptakes) for m in media}
    carbon = sorted({ex for c in conds for ex in c.exchanges})
    return {"media": med, "carbon_exchanges": carbon, "n_conditions": len(conds)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--scheme", choices=("bigg", "modelseed"), required=True)
    ap.add_argument("--org", required=True, help="Fitness Browser organism whose tables the result must serve")
    args = ap.parse_args()
    outdir = os.path.join(ROOT, "models", "curated", args.label)
    out_model = os.path.join(outdir, f"{args.label}_bigg_view.xml.gz")
    model = cobra.io.read_sbml_model(os.path.join(ROOT, args.model))
    original = model.copy()
    ns = load_namespace()
    overrides = load_overrides(args.label)
    report = {"model": args.model, "label": args.label, "scheme": args.scheme,
              "counts_before": {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)}}

    report["gpr_repairs"] = repair_gprs(model, args.label)
    original = model.copy()                 # equivalence below is checked against the published model as repaired
    report["collapsed_two_step_exchanges"] = collapse_two_step_exchanges(model)

    # metabolites
    proposals, how = {}, {}
    for met in model.metabolites:
        base, comp = split_id(met)
        c = COMPARTMENT.get(met.compartment, COMPARTMENT.get(comp))
        if c is None:
            how[met.id] = f"unknown compartment {met.compartment}"
            continue
        new_base, why = candidate(base, args.scheme, ns, overrides)
        how[met.id] = why
        if new_base:
            proposals[met.id] = (f"{new_base}_{c}", c)
    targets = collections.Counter(v[0] for v in proposals.values())
    existing = {m.id for m in model.metabolites} - set(proposals)
    clashes = sorted({k for k, v in proposals.items() if targets[v[0]] > 1 or v[0] in existing})
    renamed = {}
    for met in model.metabolites:
        if met.id in proposals and met.id not in clashes:
            new_id, c = proposals[met.id]
            if new_id != met.id:
                renamed[met.id] = new_id
            met.id = new_id
            met.compartment = c
        elif COMPARTMENT.get(met.compartment):
            met.compartment = COMPARTMENT[met.compartment]
    model.compartments = {c: {"c": "cytosol", "e": "extracellular", "p": "periplasm"}[c] for c in sorted({m.compartment for m in model.metabolites})}
    model.repair()

    # exchanges
    ex_renamed, ex_clash = {}, []
    rxn_ids = {r.id for r in model.reactions}
    for r in model.reactions:
        if r.boundary and len(r.metabolites) == 1:
            met = next(iter(r.metabolites))
            if met.compartment != "e":
                continue
            new = f"EX_{met.id}"
            if new == r.id:
                continue
            if new in rxn_ids:
                ex_clash.append([r.id, new])
                continue
            ex_renamed[r.id] = new
            rxn_ids.discard(r.id); rxn_ids.add(new)
            r.id = new
    model.repair()

    # equivalence check: same optimum with published bounds and with every exchange open
    def opt(m, open_all):
        with m:
            m.solver = "glpk"
            if open_all:
                for r in m.reactions:
                    if r.boundary:
                        r.bounds = (-1000.0, 1000.0)
            return m.slim_optimize()
    eq = {"published_bounds": [opt(original, False), opt(model, False)], "all_exchanges_open": [opt(original, True), opt(model, True)]}
    for k, (a, b) in eq.items():
        if not (abs(a - b) <= 1e-6 * max(1.0, abs(a))):
            raise SystemExit(f"translation changed the optimum ({k}): {a} vs {b}")
    report["equivalence_check"] = eq

    # what the transfer protocol will look for
    need = tables_for(args.org)
    ids = {r.id for r in model.reactions}
    report["protocol_needs"] = {
        "n_conditions_mapped": need["n_conditions"],
        "medium_exchanges_present": {m: [x for x in v if x in ids] for m, v in need["media"].items()},
        "medium_exchanges_absent": {m: [x for x in v if x not in ids] for m, v in need["media"].items()},
        "carbon_exchanges_present": [x for x in need["carbon_exchanges"] if x in ids],
        "carbon_exchanges_absent": [x for x in need["carbon_exchanges"] if x not in ids],
        "energy_check_metabolites_present": [m for m in CURRENCY_IDS if m in model.metabolites],
        "energy_check_metabolites_absent": [m for m in CURRENCY_IDS if m not in model.metabolites],
    }
    # unmapped extracellular metabolites with names, for review of the absent exchanges
    report["unmapped_extracellular"] = sorted(
        [{"id": m.id, "name": m.name, "formula": m.formula, "why": how.get(m.id, "")} for m in model.metabolites
         if m.compartment == "e" and m.id not in set(renamed.values()) and how.get(m.id) not in ("identical",)],
        key=lambda d: d["id"])
    report["metabolites_renamed"] = len(renamed)
    report["metabolite_renames"] = renamed
    report["metabolite_clashes_kept_original"] = clashes
    report["exchange_renames"] = ex_renamed
    report["exchange_rename_clashes"] = ex_clash
    report["mapping_method_counts"] = dict(collections.Counter(how.values()).most_common())
    report["overrides_used"] = overrides
    report["counts_after"] = {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)}
    model.id = f"{args.label}_bigg_view"
    os.makedirs(outdir, exist_ok=True)
    with gzip.open(out_model, "wt") as fh:
        cobra.io.write_sbml_model(model, fh)
    with open(os.path.join(outdir, f"{args.label}_translation.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=float)
        fh.write("\n")
    pn = report["protocol_needs"]
    print(json.dumps({"label": args.label, "renamed_metabolites": len(renamed), "clashes": len(clashes),
                      "exchanges_renamed": len(ex_renamed), "collapsed": len(report["collapsed_two_step_exchanges"]),
                      "equivalence": eq, "carbon_absent": pn["carbon_exchanges_absent"],
                      "medium_absent": pn["medium_exchanges_absent"],
                      "energy_absent": pn["energy_check_metabolites_absent"]}, indent=1, default=float))


if __name__ == "__main__":
    main()
