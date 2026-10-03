"""Fitness-blind correction rules for transfer to new organisms (transfer study v1).

Every rule here uses only a draft model, the organism's genome annotation (Fitness Browser gene table:
sysName, scaffold, coordinates, description) and the experimental design (media and carbon sources tested).
None reads a fitness value. The rules were developed on the development organisms (Btheta, Putida, MR1,
Smeli); docs/studies/transfer-method-v1.md states which were kept for the frozen method and why.

Rules
-----
universal_reaction_patches   universe-level reaction corrections without a taxon scope (CBMKr catabolic direction,
                             CBPS subunit rule, MECDPDH4E irreversible) from data/reference/universe_patches_v0.1.json
universal_model_additions    organism-independent additions with org '*' from data/reference/model_patches_v0.1.json
                             (glycolaldehyde diffusion sink)
add_atp_synthase_if_annotated  add the F-type ATP synthase when the draft has no ATPS* reaction and the annotation
                             names at least six of the eight canonical F0F1 subunits
normalize_conjunction_rules  a gene that is both a stand-alone alternative and a member of an AND-complex in the same
                             rule is a merge artefact of multi-source homology mapping; drop the stand-alone alternative
remove_menaquinol_if_no_pathway  drop menaquinol from the biomass when the annotation shows no menaquinone pathway
                             (fewer than two classical men or futalosine mqn steps)
"""
from __future__ import annotations

import json
import re
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import cobra
import pandas as pd

from .gene_mapping import GeneMap
from .patches import apply_model_patches, apply_universe_patches, translate_rule

# ---------------------------------------------------------------------------------------------------------------------
# annotation helpers

_ANNOT_SUFFIX = re.compile(r"\s*\((?:NCBI(?: ptt file)?|VIMSS|RefSeq|TIGR|KEGG|SEED)[^)]*\)\s*$", re.I)


def clean_desc(desc: str) -> str:
    return _ANNOT_SUFFIX.sub("", str(desc or "")).strip()


# ---------------------------------------------------------------------------------------------------------------------
# F-type ATP synthase

_ATP_INCLUDE = re.compile(r"atp synthase|h\+-transporting two-sector atpase|f0f1|fof1|f1f0", re.I)
_ATP_EXCLUDE = re.compile(r"v-type|v/a-type|a/v-type|vacuolar|flagell|\bflii\b|protein i\b|subunit i\b|dispensable|"
                          r"assembly|atp ?synthase.*(?:inhibitor|associated)|kinase|potassium|copper|cation|"
                          r"type iii|secretion|ribosom|polymerase|helicase|protease|chaperone|regulator", re.I)
_GREEK = ["alpha", "beta", "gamma", "delta", "epsilon"]
_LATIN = re.compile(r"(?:subunit|chain|f0,?|fo,?)\s+(a|b|b'|b2|c)\b|\b(a|b|b'|c)\s+(?:subunit|chain)\b", re.I)


def atp_synthase_subunits(desc: str) -> List[str]:
    """Canonical F0F1 subunit labels named in one gene description ([] if it is not an F-type subunit)."""
    d = clean_desc(desc)
    if not _ATP_INCLUDE.search(d) or _ATP_EXCLUDE.search(d):
        return []
    low = d.lower()
    found = [g for g in _GREEK if re.search(rf"\b{g}\b", low)]
    if "b-delta" in low or "b/delta" in low:
        found += ["b", "delta"]
    for m in _LATIN.finditer(d):
        tok = (m.group(1) or m.group(2) or "").lower()
        tok = "b" if tok in ("b'", "b2") else tok
        if tok in ("a", "b", "c"):
            found.append(tok)
    return sorted(set(found))


def find_atp_synthase_genes(genes: pd.DataFrame, max_gap_genes: int = 3) -> Dict[str, object]:
    """Locate annotated F0F1 subunit genes and group them into operon-like clusters.

    `genes` is a Fitness Browser gene table (sysName, scaffoldId, begin, end, desc). Returns
    {"genes": [(sysName, [labels])...], "clusters": [[sysName...]...], "complete_clusters": [...], "subunit_types": [...]}."""
    hits = []
    for _, r in genes.iterrows():
        labs = atp_synthase_subunits(r.get("desc", ""))
        if labs:
            hits.append((r["sysName"] or r["locusId"], r["scaffoldId"], int(r["begin"]), labs))
    # order genes along each scaffold to define adjacency by rank, not by base pairs
    order = {}
    for scaf, sub in genes.assign(_b=genes["begin"].astype(int)).sort_values(["scaffoldId", "_b"]).groupby("scaffoldId"):
        for i, (_, r) in enumerate(sub.iterrows()):
            order[(r["sysName"] or r["locusId"])] = (scaf, i)
    hits.sort(key=lambda h: order.get(h[0], (h[1], h[2])))
    clusters: List[List[str]] = []
    for name, scaf, _, _ in hits:
        pos = order.get(name)
        if clusters:
            last = order.get(clusters[-1][-1])
            if pos and last and pos[0] == last[0] and pos[1] - last[1] <= max_gap_genes:
                clusters[-1].append(name)
                continue
        clusters.append([name])
    labels = {name: labs for name, _, _, labs in hits}

    def types(cl):
        return sorted({t for g in cl for t in labels[g]})
    complete = [cl for cl in clusters if len(types(cl)) >= 6]
    return {"genes": [(n, labels[n]) for n, _, _, _ in hits], "clusters": clusters, "complete_clusters": complete,
            "subunit_types": sorted({t for _, _, _, labs in hits for t in labs})}


def add_atp_synthase_if_annotated(model: cobra.Model, genes: pd.DataFrame, gm: GeneMap) -> List[dict]:
    """Add ATP synthase (ATPS4rpp with a periplasm, else ATPS4r) when the draft has none and the annotation
    names >= 6 of the 8 canonical F0F1 subunit types. Rule: each operon-like cluster that is complete on its own is
    one alternative (OR of ANDs); otherwise all annotated subunit genes form one AND."""
    if any(r.id.startswith("ATPS") for r in model.reactions):
        return []
    found = find_atp_synthase_genes(genes)
    if len(found["subunit_types"]) < 6:
        return [{"reaction": "ATPS", "skipped": f"only {len(found['subunit_types'])} F0F1 subunit types annotated"}]
    if found["complete_clusters"]:
        alts = [" and ".join(cl) for cl in found["complete_clusters"]]
    else:
        alts = [" and ".join(n for n, _ in found["genes"])]
    rule = " or ".join(f"({a})" if len(alts) > 1 else a for a in alts)
    proton_out = "h_p" if "h_p" in model.metabolites else "h_e"
    rid = "ATPS4rpp" if proton_out == "h_p" else "ATPS4r"
    need = ["adp_c", "pi_c", "atp_c", "h2o_c", "h_c", proton_out]
    missing = [m for m in need if m not in model.metabolites]
    if missing:
        return [{"reaction": rid, "skipped": f"metabolites missing from the draft: {missing}"}]
    rxn = cobra.Reaction(rid, name="ATP synthase (four protons for one ATP)", lower_bound=-1000.0, upper_bound=1000.0)
    M = model.metabolites
    rxn.add_metabolites({M.adp_c: -1.0, M.get_by_id(proton_out): -4.0, M.pi_c: -1.0, M.atp_c: 1.0, M.h2o_c: 1.0, M.h_c: 3.0})
    model.add_reactions([rxn])
    rxn.gene_reaction_rule = translate_rule(rule, gm, {g.id for g in model.genes})
    rxn.annotation["gembench_rule"] = "transfer-v1 add_atp_synthase_if_annotated"
    return [{"reaction": rid, "added": rxn.reaction, "gpr": rxn.gene_reaction_rule,
             "subunit_types": found["subunit_types"], "clusters": found["clusters"]}]


# ---------------------------------------------------------------------------------------------------------------------
# gene-rule shape

def _dnf(rule: str) -> Optional[List[List[str]]]:
    """Parse a rule that is an OR of ANDs (CarveMe's form) into [[genes]...]; None if it is nested otherwise."""
    rule = rule.strip()
    if not rule:
        return []
    alts = [a.strip() for a in re.split(r"\s+or\s+(?![^()]*\))", rule)]
    out = []
    for a in alts:
        inner = a[1:-1].strip() if a.startswith("(") and a.endswith(")") else a
        if "(" in inner or ")" in inner or re.search(r"\sor\s", inner):
            return None
        out.append([g.strip() for g in re.split(r"\s+and\s+", inner) if g.strip()])
    return out


def normalize_conjunction_rules(model: cobra.Model) -> List[dict]:
    """Drop a stand-alone alternative gene when the same gene is also a member of an AND-complex alternative
    of the same rule ('A or B or (A and B)' -> '(A and B)'; 'A or C or (A and B)' -> 'C or (A and B)')."""
    applied = []
    for r in model.reactions:
        dnf = _dnf(r.gene_reaction_rule)
        if not dnf:
            continue
        complex_members = {g for alt in dnf if len(alt) > 1 for g in alt}
        if not complex_members:
            continue
        keep = [alt for alt in dnf if len(alt) > 1 or alt[0] not in complex_members]
        # de-duplicate identical complexes
        seen, uniq = set(), []
        for alt in keep:
            key = tuple(sorted(alt))
            if key not in seen:
                seen.add(key)
                uniq.append(alt)
        if len(uniq) < len(dnf):
            before = r.gene_reaction_rule
            r.gene_reaction_rule = " or ".join(("(" + " and ".join(a) + ")") if len(a) > 1 and len(uniq) > 1 else " and ".join(a)
                                               for a in uniq)
            applied.append({"reaction": r.id, "gpr_before": before, "gpr_after": r.gene_reaction_rule})
    return applied


# ---------------------------------------------------------------------------------------------------------------------
# biomass: menaquinone

_MEN_STEPS = {
    "menF_isochorismate_synthase_for_menaquinone": re.compile(r"isochorismate synthase.*menaquinone|menaquinone.*isochorismate", re.I),
    "menD_SEPHCHC_synthase": re.compile(r"2-succinyl-5-enolpyruvyl-6-hydroxy-3-cyclohexene-1-carboxyl|sephchc synthase|\bmend\b", re.I),
    "menH_SHCHC_synthase": re.compile(r"2-succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxyl|shchc synthase|\bmenh\b", re.I),
    "menC_OSB_synthase": re.compile(r"o-succinylbenzoate synthase|o-succinylbenzoic acid synthase|osb synthase|\bmenc\b", re.I),
    "menE_OSB_CoA_ligase": re.compile(r"succinylbenzoate.{0,4}coa ligase|succinylbenzoic acid.{0,4}coa ligase|osb-coa ligase|\bmene\b", re.I),
    "menB_DHNA_CoA_synthase": re.compile(r"naphthoate synthase|naphthoyl-coa synthase|dihydroxynaphthoic acid synthetase|\bmenb\b", re.I),
    "menA_DHNA_prenyltransferase": re.compile(r"naphthoate (?:octa|poly)?prenyltransferase|dhna[- ]?(?:octa|poly)?prenyltransferase|\bmena\b", re.I),
    "mqn_futalosine_pathway": re.compile(r"futalosine|1,4-dihydroxy-6-naphtho|\bmqn[abcde]\b", re.I),
}


def menaquinone_pathway_evidence(genes: pd.DataFrame) -> Dict[str, List[str]]:
    ev: Dict[str, List[str]] = {}
    for _, r in genes.iterrows():
        d = clean_desc(r.get("desc", ""))
        if re.search(r"methyltransferase", d, re.I):     # ubiE/menG serves both quinones: not diagnostic
            continue
        for step, pat in _MEN_STEPS.items():
            if pat.search(d):
                ev.setdefault(step, []).append(r["sysName"] or r["locusId"])
    return ev


def remove_menaquinol_if_no_pathway(model: cobra.Model, genes: pd.DataFrame, biomass_id: str = "Growth",
                                    quinones: Sequence[str] = ("mql8_c", "mqn8_c", "2dmmql8_c")) -> List[dict]:
    ev = menaquinone_pathway_evidence(genes)
    if len(ev) >= 2:
        return [{"reaction": biomass_id, "kept_menaquinone": True, "evidence_steps": sorted(ev)}]
    bm = model.reactions.get_by_id(biomass_id)
    removed = {}
    for mid in quinones:
        if mid in model.metabolites and model.metabolites.get_by_id(mid) in bm.metabolites:
            met = model.metabolites.get_by_id(mid)
            coef = bm.metabolites[met]
            if coef < 0:
                removed[mid] = coef
                bm.add_metabolites({met: -coef})
    return [{"reaction": biomass_id, "biomass_removed": removed, "evidence_steps": sorted(ev)}] if removed else []


# ---------------------------------------------------------------------------------------------------------------------
# universal patches from the development cycles

def universal_reaction_patches(model: cobra.Model, universe_patches: Dict, org_label: str = "_transfer_") -> List[dict]:
    """Apply only patches without a taxon scope (scoped ones encode organism-specific development findings)."""
    unscoped = {"patches": [p for p in universe_patches["patches"] if not p.get("scope_orgs")]}
    return apply_universe_patches(model, org_label, unscoped)


def universal_model_additions(model: cobra.Model, model_patches: Dict, gm: GeneMap) -> List[dict]:
    universal = dict(model_patches)
    universal["patches"] = [p for p in model_patches["patches"]
                            if (p["org"] == "*" or (isinstance(p["org"], list) and "*" in p["org"]))]
    return apply_model_patches(model, "_transfer_", universal, gm, verbose=False)


# ---------------------------------------------------------------------------------------------------------------------
# experimental design without outcomes

def conditions_from_metadata(experiments: pd.DataFrame, cs_map: pd.DataFrame, media_map: pd.DataFrame,
                             allowed_condition_2: tuple = ("", "Dimethyl Sulfoxide")) -> List[dict]:
    """Carbon-source conditions (condition_1 x media) from experiment metadata alone, for steps that must not read
    fitness tables (reference-condition choice, candidate selection). Grouping matches carbon_source_conditions."""
    cm = cs_map.set_index("condition")
    known_media = set(media_map["media"])
    sel = experiments[(experiments["expGroup"] == "carbon source") & (experiments["condition_2"].isin(allowed_condition_2))]
    out: Dict[str, dict] = {}
    for _, r in sel.iterrows():
        name, media = r["condition_1"], r["media"]
        if media not in known_media:
            continue
        key = f"{name} | {media}"
        if key not in out:
            bigg = [x for x in cm.loc[name, "bigg_ids"].split(";") if x] if name in cm.index else []
            out[key] = {"key": key, "name": name, "media": media, "bigg_ids": bigg, "n_experiments": 0}
        out[key]["n_experiments"] += 1
    return list(out.values())


CENTRAL_SUBSTRATES = ("glc__D", "fru", "glyc", "lac__L", "pyr", "succ", "ac", "fum", "mal__L", "akg", "cit", "glu__L")


def choose_reference_condition(conditions: List[dict], universe_exchanges: Iterable[str]) -> Optional[dict]:
    """Frozen rule for the gap-fill reference condition, from experiment metadata only.

    1. The main medium is the base medium with the most carbon-source experiments (ties: medium name).
    2. Candidates are the main medium's conditions with a single mapped carbon source whose exchange is in the universe.
    3. Choose the first of CENTRAL_SUBSTRATES that is a candidate (D-glucose first); otherwise the candidate with the
       most experiments, ties broken by condition name.
    The rule formalises the hand-chosen references of the four development organisms (glucose for Btheta, Putida and
    Smeli; L-lactate for MR1, which has no glucose experiment)."""
    ux = set(universe_exchanges)
    per_medium: Dict[str, int] = {}
    for c in conditions:
        per_medium[c["media"]] = per_medium.get(c["media"], 0) + c["n_experiments"]
    if not per_medium:
        return None
    main = sorted(per_medium.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
    cands = [c for c in conditions if c["media"] == main and len(c["bigg_ids"]) == 1 and f"EX_{c['bigg_ids'][0]}_e" in ux]
    if not cands:
        return None
    for sub in CENTRAL_SUBSTRATES:
        hit = [c for c in cands if c["bigg_ids"][0] == sub]
        if hit:
            return sorted(hit, key=lambda c: (-c["n_experiments"], c["name"]))[0]
    return sorted(cands, key=lambda c: (-c["n_experiments"], c["name"]))[0]


# ---------------------------------------------------------------------------------------------------------------------
# blind adjudication decisions (docs/studies/transfer-adjudication-procedure-v1.md)

_TOKEN = re.compile(r"[A-Za-z][A-Za-z0-9_.\-']*")


def apply_decisions(model: cobra.Model, decisions: Dict, gm: GeneMap, genome_sysnames: Iterable[str]) -> List[dict]:
    """Apply R6/R2/R1 decisions. Locus tags are translated to the model's gene ids where the draft carries the gene;
    a locus tag of the genome that the draft lacks is added as a new gene under its locus tag. A decision naming an
    unknown gene, or a reaction absent from the model, is not applied and is reported."""
    inv = {b: m_ for m_, b in gm.model_to_browser.items()}
    genome = set(genome_sysnames)
    model_ids = {g.id for g in model.genes}
    out = []
    for d in decisions["decisions"]:
        kind = d.get("decision")
        if kind in (None, "abstain") or not d.get("new_rule"):
            continue
        rid = d["reaction"]
        if rid not in model.reactions:
            out.append({"reaction": rid, "decision": kind, "not_applied": "reaction absent from model"})
            continue
        toks = [t for t in _TOKEN.findall(d["new_rule"]) if t.lower() not in ("and", "or")]
        unknown = [t for t in toks if t not in inv and t not in model_ids and t not in genome]
        if unknown:
            out.append({"reaction": rid, "decision": kind, "not_applied": f"unknown genes {unknown}"})
            continue
        new_rule = _TOKEN.sub(lambda mm: mm.group(0) if mm.group(0).lower() in ("and", "or") or mm.group(0) in model_ids
                              else inv.get(mm.group(0), mm.group(0)), d["new_rule"])
        r = model.reactions.get_by_id(rid)
        before = r.gene_reaction_rule
        try:
            r.gene_reaction_rule = new_rule
        except Exception as e:  # noqa: BLE001  malformed boolean expression
            out.append({"reaction": rid, "decision": kind, "not_applied": f"unparseable rule: {e}"})
            r.gene_reaction_rule = before
            continue
        added = [t for t in toks if t not in inv and t not in model_ids]
        out.append({"reaction": rid, "decision": kind, "gpr_before": before, "gpr_after": r.gene_reaction_rule,
                    "genes_added": added})
    return out
