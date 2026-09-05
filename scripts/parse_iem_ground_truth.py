"""Parse the whole-body-model IEM biomarker benchmark out of the COBRA Toolbox script runIEM_HH.m
(Thiele et al. 2020, Mol Syst Biol; 57 IEMs, 252 biofluid-specific biomarkers) into a
machine-readable tuple table.

The script is chunked at each `checkIEM_WBM(` call. Within a chunk: the last `%%` comment
line carrying a gene id is the IEM header (Entrez gene id(s) + name + abbreviation), the
`R = ...` assignments are the reaction patterns knocked out, and the BiomarkerRxns cell
arrays list (biomarker reaction, 'Direction (biofluids)') pairs — the first array is the
whole-body-model branch (urine exchanges + blood demand reactions), the second (if any)
the Recon3D branch.

Output: results/iem_ground_truth/iem_biomarkers_v0.tsv + provenance.json
"""
from __future__ import annotations

import json
import os
import re
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.join(ROOT, "external", "cobratoolbox")
SRC = os.path.join(REPO, "src", "analysis", "wholeBody", "PSCMToolbox", "runIEM_HH.m")
OUT = os.path.join(ROOT, "results", "iem_ground_truth")

PAIR = re.compile(r"'([^']+)'\s+'([^']+)'")
MET = re.compile(r"^(?:EX_|DM_)(?P<met>.+?)\[(?P<comp>[a-zA-Z]+)\]$")
FLUID = {"u": "urine", "bc": "blood", "csf": "csf", "bd": "bile", "fe": "feces", "d": "diet", "e": "extracellular", "c": "cytosol"}


def parse_header(h: str):
    h = h.strip().lstrip("%").strip()
    h = h.replace("gene ID:", "").replace("gene ID", "")
    genes = re.findall(r"\d+[.,]\d", h)
    genes = [g.replace(",", ".") for g in genes]
    name = re.sub(r"'?\d+[.,]\d'?", "", h)
    name = re.sub(r"^[\s:%\-']+", "", name).strip().strip("'").strip()
    name = re.sub(r"\s+", " ", name.replace("\t", " "))
    return genes, name


def main() -> None:
    lines = open(SRC, encoding="utf-8", errors="replace").read().splitlines()
    calls = [i for i, l in enumerate(lines) if "checkIEM_WBM(" in l]
    rows = []
    prev = 0
    n_iem_calls = 0
    for k, ci in enumerate(calls):
        chunk = lines[prev:ci + 1]
        abbr = re.search(r"IEMSol_(\w+)", lines[ci]).group(1)
        headers = [(prev + j + 1, l) for j, l in enumerate(chunk) if re.match(r"\s*%%", l) and re.search(r"\d+[.,]\d", l)]
        hdr_line, hdr = headers[-1] if headers else (None, "")
        genes, name = parse_header(hdr) if hdr else ([], "")
        R = []
        for l in chunk:
            m = re.match(r"\s*R\s*=\s*(\{[^}]*\}|'[^']*')", l)
            if m:
                R += re.findall(r"'([^']+)'", m.group(1))
        txt = "\n".join(chunk)
        sets = [PAIR.findall(m.group(1)) for m in re.finditer(r"BiomarkerRxns\s*=\s*\{(.*?)\}", txt, re.S)]
        wbm = sets[0] if sets else []
        recon = sets[1] if len(sets) > 1 else []
        recon_keys = {(r, l) for r, l in recon}
        n_iem_calls += 1
        seen = set()
        for rxn, lab in wbm:
            mm = MET.match(rxn.strip())
            met, comp = (mm.group("met"), mm.group("comp")) if mm else (rxn, "")
            m2 = re.match(r"\s*(Increased|Decreased|Unchanged|Normal|Unknown)\s*(?:\((.*?)\))?", lab, re.I)
            direction = m2.group(1).capitalize() if m2 else "Unparsed"
            fluids_raw = (m2.group(2) or "").strip().lower() if m2 else ""
            key = (met, comp, direction)
            if key in seen:
                continue
            seen.add(key)
            rows.append({
                "iem_abbr": abbr, "call_index": k + 1, "iem_name": name, "entrez_gene_ids": ";".join(genes),
                "knockout_reaction_patterns": ";".join(dict.fromkeys(R)), "biomarker_reaction": rxn.strip(),
                "vmh_metabolite": met, "model_compartment": comp, "biofluid_tested": FLUID.get(comp, comp),
                "expected_direction": direction, "reported_biofluids_raw": fluids_raw, "label_raw": lab.strip(),
                "also_in_recon3d_branch": (rxn, lab) in recon_keys,
                "source_header_line": hdr_line, "source_call_line": ci + 1,
            })
        prev = ci + 1

    os.makedirs(OUT, exist_ok=True)
    cols = list(rows[0].keys())
    with open(os.path.join(OUT, "iem_biomarkers_v0.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    try:
        commit = subprocess.check_output(["git", "-C", REPO, "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        commit = "unknown"
    abbrs = [r["iem_abbr"] for r in rows]
    prov = {
        "source_file": "src/analysis/wholeBody/PSCMToolbox/runIEM_HH.m", "repository": "github.com/opencobra/cobratoolbox",
        "commit": commit, "script_author_note": "Ines Thiele 2018-2019 (file header)",
        "primary_publication": "Thiele et al. 2020, Mol Syst Biol 16:e8982 — 57 IEMs, 252 biofluid-specific biomarkers, "
                               "85.3% (Harvetta) / 84.9% (Harvey) direction agreement",
        "licence": "COBRA Toolbox: GPLv3 (code). The biomarker knowledge was compiled from the literature by the authors; "
                   "cite the paper when reusing.",
        "n_checkIEM_calls": len(calls), "n_rows": len(rows), "n_unique_iem_abbr": len(set(abbrs)),
        "n_unique_metabolites": len({r["vmh_metabolite"] for r in rows}),
        "by_biofluid_tested": {f: sum(1 for r in rows if r["biofluid_tested"] == f) for f in sorted({r["biofluid_tested"] for r in rows})},
        "by_direction": {d: sum(1 for r in rows if r["expected_direction"] == d) for d in sorted({r["expected_direction"] for r in rows})},
        "rows_per_iem": {a: abbrs.count(a) for a in sorted(set(abbrs))},
        "caveats": [
            "The paper reports 57 IEMs and 252 biofluid-specific biomarkers; the script contains 63 checkIEM_WBM calls "
            "(some IEMs run twice, e.g. MMA with two causal genes; a few abbreviations repeat). Reconcile against the "
            "paper's supplementary table before quoting counts.",
            "Rows are unique per (metabolite, model compartment, direction) within a call; [u] = urine exchange, "
            "[bc] = blood demand reaction, [csf] = CSF.",
            "Gene ids are Entrez ids as written in the script headers (e.g. '4967.1'); at least one header id looks "
            "mis-assigned (4967.1 used for both AKGD and ADSL) — flagged for curation, not corrected.",
            "Not yet linked to OMIM/Orphanet/HPO identifiers; VMH metabolite ids need mapping to HMDB/ChEBI.",
        ],
    }
    with open(os.path.join(OUT, "provenance.json"), "w") as fh:
        json.dump(prov, fh, indent=2)
    print(json.dumps({k: v for k, v in prov.items() if k != "rows_per_iem"}, indent=2))


if __name__ == "__main__":
    main()
