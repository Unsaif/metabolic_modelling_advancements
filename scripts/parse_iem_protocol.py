"""Extract the per-IEM simulation protocol (include/exclude reaction patterns, per-block bound
tweaks, biomarker reactions with expected direction) from runIEM_HH.m into JSON for the Python port.
"""
from __future__ import annotations

import json
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "external", "cobratoolbox", "src", "analysis", "wholeBody", "PSCMToolbox", "runIEM_HH.m")
OUT = os.path.join(ROOT, "results", "iem_ground_truth", "iem_protocol_v0.json")

PAIR = re.compile(r"'([^']+)'\s+'([^']+)'")
STRFIND_VAR = re.compile(r"^\s*(\w+)\s*=\s*model\.rxns\(.*strfind\(model\.rxns,\s*'([^']+)'\)")
SETDIFF = re.compile(r"IEMRxns\s*=\s*setdiff\((\w+)\s*,\s*(\w+)\)")
BOUND = re.compile(r"^\s*model(?:IEM)?\.(lb|ub)\(ismember\(model(?:IEM)?\.rxns,\s*(\w+)\)\)\s*=\s*([-\d.eE]+)\s*;")


def main() -> None:
    lines = open(SRC, encoding="utf-8", errors="replace").read().splitlines()
    calls = [i for i, l in enumerate(lines) if "checkIEM_WBM(" in l]
    prev = 0
    out = []
    for k, ci in enumerate(calls):
        chunk = lines[prev:ci + 1]
        abbr = re.search(r"IEMSol_(\w+)", lines[ci]).group(1)
        include = []
        txt0 = "\n".join(chunk)
        r_assignments = list(re.finditer(r"^\s*R\s*=\s*(\{.*?\}|'[^']*')\s*;", txt0, re.S | re.M))
        if r_assignments:   # the last assignment in the chunk is the block's own (the first chunk also contains the global preamble)
            include += re.findall(r"'([^']+)'", r_assignments[-1].group(1))
        strvars = {}
        for l in chunk:
            m = STRFIND_VAR.match(l)
            if m:
                strvars.setdefault(m.group(1), []).append(m.group(2))
        exclude = []
        for l in chunk:
            m = SETDIFF.search(l)
            if m and m.group(2) in strvars:
                exclude += strvars[m.group(2)]
        tweaks = []
        for l in chunk:
            m = BOUND.match(l)
            if m and m.group(2) in strvars:
                for pat in strvars[m.group(2)]:
                    tweaks.append({"pattern": pat, "bound": m.group(1), "value": float(m.group(3))})
        txt = "\n".join(chunk)
        sets = [PAIR.findall(m.group(1)) for m in re.finditer(r"BiomarkerRxns\s*=\s*\{(.*?)\}", txt, re.S)]
        wbm = [(r.strip(), l.strip()) for r, l in (sets[0] if sets else [])]
        # if the R patterns are empty the block used an explicit list (e.g. BTD, FED, TETB): capture IEMRxns = {...}
        if not include:
            m = re.search(r"IEMRxns\s*=\s*\{([^}]*)\}", txt, re.S)
            if m:
                include = re.findall(r"'([^']+)'", m.group(1))
        out.append({"iem": abbr, "call_index": k + 1, "include_patterns": list(dict.fromkeys(include)),
                    "exclude_patterns": list(dict.fromkeys(exclude)), "bound_tweaks": tweaks, "biomarkers": wbm,
                    "source_call_line": ci + 1})
        prev = ci + 1
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2)
    n_empty = [o["iem"] for o in out if not o["include_patterns"]]
    print(f"{len(out)} IEM calls; {sum(len(o['biomarkers']) for o in out)} biomarker entries; empty include lists: {n_empty}")
    for o in out[:6]:
        print(o["iem"], o["include_patterns"], "excl", o["exclude_patterns"], "tweaks", o["bound_tweaks"], "n_bm", len(o["biomarkers"]))


if __name__ == "__main__":
    main()
