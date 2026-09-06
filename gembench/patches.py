"""Apply universe-level reaction patches and per-organism gene-rule patches to a model.

Patch files live in data/reference/ (universe_patches_v*.json, gpr_patches_v*.json); every patch carries a
biochemical rationale written before re-scoring and, for gene-rule patches, a status ('accepted' or 'held')
set by the fitness-benchmark verifier. Reaction patches may be restricted to organisms with `scope_orgs`.
"""
from __future__ import annotations

import json
import re
from typing import Dict, Iterable, List, Optional

import cobra

from .gene_mapping import GeneMap

_TOKEN = re.compile(r"[A-Za-z][A-Za-z0-9_]*")


def load_patch_files(paths: Iterable[str]) -> Dict:
    """Merge one or more gene-rule patch files (applied in order); version strings are joined with '+'."""
    parts = [json.load(open(f)) for f in paths]
    return {"version": "+".join(p["version"] for p in parts), "patches": [x for p in parts for x in p["patches"]]}


def apply_universe_patches(model: cobra.Model, org: str, patches: Dict) -> List[dict]:
    """Reaction-level patches: bounds or gene-rule shape changes on reactions the model contains."""
    applied: List[dict] = []
    for pt in patches["patches"]:
        if pt["reaction"] not in model.reactions:
            continue
        if pt.get("scope_orgs") and org not in pt["scope_orgs"]:
            continue        # taxon-scoped patch (e.g. URIC closed for E. coli, PPCK reversible for Bacteroides)
        r = model.reactions.get_by_id(pt["reaction"])
        if "gpr" in pt["change"]:
            rule = r.gene_reaction_rule
            alts = [a.strip() for a in re.split(r"\s+or\s+(?![^()]*\))", rule)]
            conj = [a for a in alts if " and " in a]
            members = set()
            for a in conj:
                members |= set(a.strip("()").split(" and "))
            if pt["change"]["gpr"] == "keep_only_conjunctions":
                keep = [a for a in alts if " and " in a]
            else:   # drop_standalone_members_of_conjunctions
                keep = [a for a in alts if " and " in a or a.strip("()") not in members]
            if conj and len(keep) < len(alts):
                r.gene_reaction_rule = " or ".join(keep)
                applied.append({"reaction": pt["reaction"], "gpr_before": rule, "gpr_after": r.gene_reaction_rule})
        elif pt["reaction"] != "CBMKr" or "CBPS" in model.reactions:
            # the carbamate-kinase restriction presumes a carbamoyl-phosphate synthetase (see the patch caveat)
            before = r.bounds
            r.bounds = (pt["change"]["lower_bound"], pt["change"]["upper_bound"])
            applied.append({"reaction": pt["reaction"], "bounds_before": list(before), "bounds_after": list(r.bounds)})
    return applied


def translate_rule(rule: str, gm: GeneMap, model_ids: set) -> str:
    """Browser locus tags -> model gene ids; tokens that already are model ids (or unknown, R4 additions) stay."""
    inv = {b: m_ for m_, b in gm.model_to_browser.items()}
    return _TOKEN.sub(lambda mm: mm.group(0) if (mm.group(0) in ("and", "or") or mm.group(0) in model_ids)
                      else inv.get(mm.group(0), mm.group(0)), rule)


def apply_gpr_patches(model: cobra.Model, org: str, gm: GeneMap, gpr_patches: Dict, statuses: Iterable[str] = ("accepted",),
                      only: Optional[Iterable[tuple]] = None, verbose: bool = True) -> List[dict]:
    """Per-organism gene-rule patches. `statuses` selects by verifier status ('accepted' is the default; a patch
    without a status counts as accepted); `only` restricts to (org, reaction) pairs."""
    inv = {b: m_ for m_, b in gm.model_to_browser.items()}
    model_ids = {g.id for g in model.genes}
    sel = set(only) if only is not None else None
    applied: List[dict] = []
    for pt in gpr_patches["patches"]:
        if pt["org"] != org or pt.get("status", "accepted") not in statuses:
            continue
        if sel is not None and (org, pt["reaction"]) not in sel:
            continue
        if pt["reaction"] not in model.reactions:
            continue
        r = model.reactions.get_by_id(pt["reaction"])
        rule = pt["new_rule"]
        toks = [g for g in _TOKEN.findall(rule) if g not in ("and", "or")]
        new_genes = [g for g in toks if g not in inv and g not in model_ids]
        if new_genes and pt.get("rule") != "R4":
            if verbose:
                print(f"   gpr patch {pt['reaction']}: genes not in model map {new_genes}; skipped", flush=True)
            continue
        new_rule = translate_rule(rule, gm, model_ids)
        applied.append({"reaction": pt["reaction"], "gpr_before": r.gene_reaction_rule, "gpr_after": new_rule,
                        "rule": pt.get("rule"), "status": pt.get("status", "accepted"), "genes_added": new_genes})
        r.gene_reaction_rule = new_rule
    return applied


def apply_model_patches(model: cobra.Model, org: str, model_patches: Dict, gm: GeneMap, verbose: bool = True) -> List[dict]:
    """Per-organism reaction additions with gene evidence (data/reference/model_patches_v*.json).

    Metabolites missing from the model are created from the file's `new_metabolites` table; gene ids in the
    patch are Fitness Browser locus tags, translated to model ids where the draft carries the gene and otherwise
    added under the locus tag. `org` may be one organism, a list, or "*" (every model); a patch may be
    conditional on `requires_metabolites` (all present) and `skip_if_reactions` (none present)."""
    applied: List[dict] = []
    new_mets = model_patches.get("new_metabolites", {})
    model_ids = {g.id for g in model.genes}
    for pt in model_patches["patches"]:
        orgs = pt["org"] if isinstance(pt["org"], list) else [pt["org"]]
        if "*" not in orgs and org not in orgs:
            continue
        if pt["reaction"] in model.reactions:
            if verbose:
                print(f"   model patch {pt['reaction']}: already in the model; skipped", flush=True)
            continue
        if any(mid not in model.metabolites for mid in pt.get("requires_metabolites", [])):
            continue        # the model has no use for it (e.g. no glycolaldehyde at all)
        if any(rid in model.reactions for rid in pt.get("skip_if_reactions", [])):
            continue        # the model already has an equivalent route
        rxn = cobra.Reaction(pt["reaction"], name=pt.get("name", pt["reaction"]),
                             lower_bound=pt.get("lower_bound", 0.0), upper_bound=pt.get("upper_bound", 1000.0))
        stoich = {}
        for mid, coef in pt["metabolites"].items():
            if mid not in model.metabolites:
                spec = new_mets.get(mid)
                if spec is None:
                    raise KeyError(f"model patch {pt['reaction']}: metabolite {mid} is neither in the model nor in new_metabolites")
                model.add_metabolites([cobra.Metabolite(mid, name=spec["name"], formula=spec.get("formula"),
                                                        charge=spec.get("charge"), compartment=spec["compartment"])])
            stoich[model.metabolites.get_by_id(mid)] = coef
        rxn.add_metabolites(stoich)
        model.add_reactions([rxn])
        gpr = translate_rule(pt.get("gpr", ""), gm, model_ids) if pt.get("gpr") else ""
        rxn.gene_reaction_rule = gpr
        rxn.annotation["gembench_patch"] = f"model_patches v{model_patches['version']}"
        applied.append({"reaction": pt["reaction"], "added": rxn.reaction, "gpr": gpr, "bounds": list(rxn.bounds)})
    return applied
