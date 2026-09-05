"""Link the parsed whole-body-model IEM biomarker table to Orphanet, OMIM, HPO and HMDB/ChEBI,
and cross-check each (IEM, metabolite, biofluid, direction) tuple against HPO phenotype annotations.

Inputs (all fetched from public, licence-clean sources; see provenance in the output):
  results/iem_ground_truth/iem_biomarkers_v0.tsv      (parsed from runIEM_HH.m, COBRA Toolbox)
  data/humangem_genes.tsv                              (Human-GEM: Entrez -> symbol)
  data/en_product6.xml, en_product1.xml, en_product4.xml (Orphanet: genes, cross-refs, HPO phenotypes; CC-BY-4.0)
  data/hp.obo                                          (HPO ontology)
  external/COBRA.models/mat/Recon3DModel_301.mat       (VMH metabolite names, HMDB/ChEBI ids)

Output: results/iem_ground_truth/iem_biomarkers_v0.1_linked.tsv, hpo_metabolite_tuples.tsv, linking_report.json
"""
from __future__ import annotations

import csv
import difflib
import json
import os
import re
import xml.etree.ElementTree as ET
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results", "iem_ground_truth")


# ---------------------------------------------------------------- loaders
def load_entrez_to_symbol():
    m = {}
    with open(os.path.join(DATA, "humangem_genes.tsv")) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            for e in str(r["geneEntrezID"]).split(";"):
                e = e.strip()
                if e:
                    m[e] = r["geneSymbols"]
    return m


def load_orphanet_genes():
    """symbol -> list of (orpha, disorder_name, assoc_type, assoc_status)"""
    tree = ET.parse(os.path.join(DATA, "en_product6.xml"))
    out = defaultdict(list)
    names = {}
    for d in tree.getroot().iter("Disorder"):
        orpha = d.findtext("OrphaCode")
        name = d.findtext("Name")
        names[orpha] = name
        for a in d.iter("DisorderGeneAssociation"):
            sym = a.findtext("Gene/Symbol")
            atype = a.findtext("DisorderGeneAssociationType/Name") or ""
            status = a.findtext("DisorderGeneAssociationStatus/Name") or ""
            out[sym].append((orpha, name, atype, status))
    return out, names


def load_orphanet_xrefs():
    """orpha -> {'OMIM': [...], 'ICD-10': [...], synonyms: [...]}"""
    tree = ET.parse(os.path.join(DATA, "en_product1.xml"))
    out = {}
    for d in tree.getroot().iter("Disorder"):
        orpha = d.findtext("OrphaCode")
        refs = defaultdict(list)
        for x in d.iter("ExternalReference"):
            refs[x.findtext("Source")].append(x.findtext("Reference"))
        syn = [s.text for s in d.iter("Synonym") if s.text]
        out[orpha] = {"xrefs": dict(refs), "synonyms": syn, "name": d.findtext("Name")}
    return out


def load_orphanet_hpo():
    """orpha -> list of (hpo_id, term, frequency)"""
    tree = ET.parse(os.path.join(DATA, "en_product4.xml"))
    out = defaultdict(list)
    for d in tree.getroot().iter("Disorder"):
        orpha = d.findtext("OrphaCode")
        for a in d.iter("HPODisorderAssociation"):
            out[orpha].append((a.findtext("HPO/HPOId"), a.findtext("HPO/HPOTerm"), a.findtext("HPOFrequency/Name") or ""))
    return out


def load_vmh_metabolites():
    import cobra
    m = cobra.io.load_matlab_model(os.path.join(ROOT, "external", "COBRA.models", "mat", "Recon3DModel_301.mat"))
    out = {}
    for met in m.metabolites:
        base = re.sub(r"\[[a-z]+\]$", "", met.id)
        if base in out:
            continue
        ann = met.annotation
        out[base] = {"name": met.name, "hmdb": ";".join(ann.get("hmdb", [])) if isinstance(ann.get("hmdb"), list) else str(ann.get("hmdb", "")),
                     "chebi": ";".join(ann.get("chebi", [])) if isinstance(ann.get("chebi"), list) else str(ann.get("chebi", "")),
                     "kegg": ";".join(ann.get("kegg.compound", [])) if isinstance(ann.get("kegg.compound"), list) else str(ann.get("kegg.compound", "")),
                     "formula": met.formula or ""}
    return out


# ------------------------------------------------- HPO term -> metabolite tuple
UP = r"(?:Elevated|Increased|High|Raised)"
DOWN = r"(?:Decreased|Reduced|Low|Diminished)"
FLUIDS = {"circulating": "blood", "plasma": "blood", "serum": "blood", "blood": "blood", "urinary": "urine", "urine": "urine",
          "csf": "csf", "cerebrospinal fluid": "csf"}


def hpo_term_to_tuple(term: str):
    """Heuristic parse of an HPO term name into (metabolite, biofluid, direction) or None."""
    t = term.strip()
    m = re.match(rf"^(?P<dir>{UP}|{DOWN}) (?P<fluid>circulating|plasma|serum|blood|urinary|urine|CSF|cerebrospinal fluid) (?P<met>.+?)(?: concentration| level| levels| excretion)?$", t, re.I)
    if m:
        d = "Increased" if re.match(UP, m.group("dir"), re.I) else "Decreased"
        return m.group("met").strip().lower(), FLUIDS[m.group("fluid").lower()], d
    m = re.match(rf"^(?P<dir>{UP}|{DOWN}) (?P<met>.+?) (?P<fluid>excretion in urine|in urine|in blood|in plasma|in serum|in CSF)$", t, re.I)
    if m:
        d = "Increased" if re.match(UP, m.group("dir"), re.I) else "Decreased"
        fl = "urine" if "urine" in m.group("fluid").lower() else ("csf" if "csf" in m.group("fluid").lower() else "blood")
        return m.group("met").strip().lower(), fl, d
    m = re.match(r"^(?P<met>[A-Za-z0-9,\-\(\) ]+?) ?aciduria$", t, re.I)
    if m:
        return (m.group("met").strip().lower() + " acid"), "urine", "Increased"
    m = re.match(r"^(?P<met>[A-Za-z0-9,\-\(\) ]+?)uria$", t, re.I)
    if m and len(m.group("met")) > 3 and not re.search(r"(hemat|prote|glucos|keton|album|micro|olig|poly|noct|hyposthen|isosthen|pne|py|dys)", m.group("met"), re.I):
        return m.group("met").strip().lower(), "urine", "Increased"
    m = re.match(r"^(?P<dir>Hyper|Hypo)(?P<met>[a-z0-9\-]+?)(?:aemia|emia)$", t, re.I)
    if m:
        stem = m.group("met").lower()
        return EMIA_STEMS.get(stem, stem), "blood", ("Increased" if m.group("dir").lower() == "hyper" else "Decreased")
    return None


EMIA_STEMS = {"glyc": "glucose", "ammon": "ammonia", "kal": "potassium", "natr": "sodium", "calc": "calcium", "phosphat": "phosphate",
              "uric": "uric acid", "lact": "lactate", "lactat": "lactate", "chol": "cholesterol", "cholesterol": "cholesterol",
              "triglycerid": "triglyceride", "bilirubin": "bilirubin", "glutamin": "glutamine", "argin": "arginine", "lysin": "lysine",
              "ornithin": "ornithine", "citrullin": "citrulline", "prolin": "proline", "methionin": "methionine", "phenylalanin": "phenylalanine",
              "tyrosin": "tyrosine", "homocystein": "homocysteine", "histidin": "histidine", "valin": "valine", "leucin": "leucine",
              "isoleucin": "isoleucine", "alanin": "alanine", "glycin": "glycine", "serin": "serine", "threonin": "threonine",
              "keton": "ketone bodies", "aminoacid": "amino acids", "insulin": "insulin", "magnes": "magnesium", "carotin": "carotene",
              "caroten": "carotene", "lipid": "lipids", "glycerol": "glycerol", "ferrit": "ferritin", "uricac": "uric acid", "oxal": "oxalate",
              "sarcosin": "sarcosine", "beta-alanin": "beta-alanine", "carnitin": "carnitine", "creatin": "creatine", "cystin": "cystine",
              "galactos": "galactose", "fructos": "fructose", "pyruv": "pyruvate", "succin": "succinate", "glutam": "glutamate",
              "orotic": "orotic acid", "xanthin": "xanthine", "guanidinoacet": "guanidinoacetate", "tryptophan": "tryptophan",
              "asparagin": "asparagine", "aspart": "aspartate", "thyroxin": "thyroxine", "bicarbonat": "bicarbonate", "chlor": "chloride"}


def norm_met(name: str) -> str:
    n = name.lower()
    n = re.sub(r"^(l|d|dl)-", "", n)
    n = n.replace("(r)-", "").replace("(s)-", "")
    n = re.sub(r"[^a-z0-9]+", "", n)
    return n


def met_match(vmh_name: str, hpo_met: str) -> bool:
    a, b = norm_met(vmh_name), norm_met(hpo_met)
    if not a or not b:
        return False
    if a == b or a in b or b in a:
        return True
    # acid/ate equivalence: 'glutaric acid' vs 'glutarate'
    a2 = re.sub(r"(ate|icacid|acid)$", "", a); b2 = re.sub(r"(ate|icacid|acid)$", "", b)
    return len(a2) > 3 and (a2 == b2 or a2 in b2 or b2 in a2)


# ------------------------------------------------------------------- main
def main() -> None:
    rows = list(csv.DictReader(open(os.path.join(OUT, "iem_biomarkers_v0.tsv")), delimiter="\t"))
    e2s = load_entrez_to_symbol()
    gene2dis, orpha_names = load_orphanet_genes()
    xrefs = load_orphanet_xrefs()
    orpha_hpo = load_orphanet_hpo()
    vmh = load_vmh_metabolites()

    # choose one Orphanet disorder per (IEM call): best name similarity among the gene's disorders
    iem_link = {}
    for r in rows:
        key = (r["iem_abbr"], r["call_index"])
        if key in iem_link:
            continue
        genes = [g.split(".")[0] for g in r["entrez_gene_ids"].split(";") if re.match(r"^\d+\.\d$", g.strip())]
        symbols = [e2s.get(g, "") for g in genes]
        cands = []
        for s in symbols:
            for orpha, dname, atype, status in gene2dis.get(s, []):
                syns = xrefs.get(orpha, {}).get("synonyms", [])
                score = max([difflib.SequenceMatcher(None, r["iem_name"].lower(), x.lower()).ratio() for x in [dname] + syns] or [0])
                cands.append((score, orpha, dname, atype, status, s))
        cands.sort(reverse=True)
        best = cands[0] if cands else None
        iem_link[key] = {"gene_symbols": ";".join(x for x in symbols if x), "unmapped_entrez": ";".join(g for g, s in zip(genes, symbols) if not s),
                         "orpha_code": best[1] if best else "", "orpha_name": best[2] if best else "", "orpha_match_score": round(best[0], 3) if best else "",
                         "orpha_assoc_type": best[3] if best else "", "orpha_assoc_status": best[4] if best else "",
                         "n_candidate_disorders": len(cands),
                         "omim_ids": ";".join(xrefs.get(best[1], {}).get("xrefs", {}).get("OMIM", [])) if best else "",
                         "icd10": ";".join(xrefs.get(best[1], {}).get("xrefs", {}).get("ICD-10", [])) if best else ""}

    # HPO metabolite tuples per chosen disorder
    hpo_tuples = defaultdict(list)   # orpha -> [(met, fluid, dir, hpo_id, term, freq)]
    for key, link in iem_link.items():
        o = link["orpha_code"]
        if not o:
            continue
        for hid, term, freq in orpha_hpo.get(o, []):
            tup = hpo_term_to_tuple(term)
            if tup:
                hpo_tuples[o].append((tup[0], tup[1], tup[2], hid, term, freq))

    out_rows = []
    stats = defaultdict(int)
    for r in rows:
        link = iem_link[(r["iem_abbr"], r["call_index"])]
        v = vmh.get(r["vmh_metabolite"], {})
        status, evid = "no_hpo_metabolite_term_for_disorder", ""
        o = link["orpha_code"]
        if o and hpo_tuples.get(o):
            status = "not_in_hpo_annotations"
            for met, fl, d, hid, term, freq in hpo_tuples[o]:
                if met_match(v.get("name", r["vmh_metabolite"]), met):
                    if fl == r["biofluid_tested"] and d == r["expected_direction"]:
                        status, evid = "corroborated", f"{hid} {term} [{freq}]"; break
                    elif d == r["expected_direction"]:
                        status, evid = "corroborated_other_biofluid", f"{hid} {term} [{freq}]"
                    elif fl == r["biofluid_tested"]:
                        status, evid = "direction_conflict", f"{hid} {term} [{freq}]"
                    elif status == "not_in_hpo_annotations":
                        status, evid = "metabolite_mentioned_other_biofluid_and_direction", f"{hid} {term} [{freq}]"
        stats[status] += 1
        out_rows.append({**r, **link, "vmh_name": v.get("name", ""), "hmdb": v.get("hmdb", ""), "chebi": v.get("chebi", ""),
                         "kegg": v.get("kegg", ""), "hpo_check": status, "hpo_evidence": evid})

    cols = list(out_rows[0].keys())
    with open(os.path.join(OUT, "iem_biomarkers_v0.1_linked.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out_rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # HPO-derived tuples not present in the lab table (candidates for extension)
    lab_keys = defaultdict(set)
    for r in out_rows:
        lab_keys[r["orpha_code"]].add((norm_met(r["vmh_name"] or r["vmh_metabolite"]), r["biofluid_tested"], r["expected_direction"]))
    extra = []
    for key, link in iem_link.items():
        o = link["orpha_code"]
        for met, fl, d, hid, term, freq in hpo_tuples.get(o, []):
            if not any(met_match(k[0], met) and k[1] == fl and k[2] == d for k in lab_keys[o]):
                extra.append({"iem_abbr": key[0], "orpha_code": o, "orpha_name": link["orpha_name"], "hpo_id": hid, "hpo_term": term,
                              "frequency": freq, "metabolite_parsed": met, "biofluid": fl, "direction": d})
    with open(os.path.join(OUT, "hpo_metabolite_tuples_not_in_lab_table.tsv"), "w") as fh:
        cols2 = list(extra[0].keys()) if extra else []
        fh.write("\t".join(cols2) + "\n")
        for r in extra:
            fh.write("\t".join(str(r[c]) for c in cols2) + "\n")

    report = {
        "n_rows": len(out_rows), "n_iem_calls": len(iem_link),
        "n_iem_calls_with_gene_symbol": sum(1 for l in iem_link.values() if l["gene_symbols"]),
        "n_iem_calls_with_orpha": sum(1 for l in iem_link.values() if l["orpha_code"]),
        "n_iem_calls_with_omim": sum(1 for l in iem_link.values() if l["omim_ids"]),
        "n_iem_calls_with_hpo_metabolite_terms": sum(1 for l in iem_link.values() if hpo_tuples.get(l["orpha_code"])),
        "orpha_match_score_below_0.5": [(k[0], l["orpha_name"], l["orpha_match_score"]) for k, l in iem_link.items() if l["orpha_code"] and l["orpha_match_score"] < 0.5],
        "unmapped_entrez": sorted({l["unmapped_entrez"] for l in iem_link.values() if l["unmapped_entrez"]}),
        "hpo_check_counts": dict(stats),
        "n_hpo_tuples_not_in_lab_table": len(extra),
        "vmh_metabolites_without_hmdb": sorted({r["vmh_metabolite"] for r in out_rows if not r["hmdb"]}),
        "sources": {"orphanet": "Orphadata (github.com/Orphanet/Orphadata_aggregated, CC-BY-4.0; product1/4/6 dated in file headers)",
                    "hpo": "HPO release 2026-09-01 (hp.obo) and phenotype.hpoa 2026-09-02",
                    "gene_mapping": "Human-GEM 2.0.0 model/genes.tsv (Entrez -> symbol)",
                    "metabolite_ids": "Recon3DModel_301.mat (opencobra/COBRA.models) metabolite annotations"},
        "method_caveats": ["Orphanet disorder chosen per IEM by best name similarity among disorders linked to the causal gene; low scores flagged.",
                           "HPO term parsing is heuristic (regex on term names); 'corroborated' means an HPO term for the chosen disorder names the same "
                           "metabolite with the same biofluid and direction; absence is not evidence against the lab's tuple.",
                           "HPO annotations describe patients; the model tests exchange/demand reactions — a biofluid mismatch is expected for some."],
    }
    with open(os.path.join(OUT, "linking_report.json"), "w") as fh:
        json.dump(report, fh, indent=2)
    print(json.dumps(report, indent=2)[:6000])


if __name__ == "__main__":
    main()
