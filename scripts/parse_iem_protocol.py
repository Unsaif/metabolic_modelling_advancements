"""Extract the per-IEM simulation protocol (include/exclude reaction patterns, per-block bound
tweaks, biomarker reactions with expected direction) from runIEM_HH.m into JSON for the Python port.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "external", "cobratoolbox", "src", "analysis", "wholeBody", "PSCMToolbox", "runIEM_HH.m")
OUT = os.path.join(ROOT, "results", "iem_ground_truth", "iem_protocol_v0.2.json")

PAIR = re.compile(r"'([^']+)'\s+'([^']+)'")
STRFIND_VAR = re.compile(r"^\s*(\w+)\s*=\s*model\.rxns\(.*strfind\(model\.rxns,\s*'([^']+)'\)")
SETDIFF = re.compile(r"IEMRxns\s*=\s*setdiff\((\w+)\s*,\s*(\w+)\)")
BOUND = re.compile(r"^\s*model(?:IEM)?\.(lb|ub)\((?:find\()?ismember\(model(?:IEM)?\.rxns,\s*(\w+)\)\)?\)\s*=\s*([-\d.eE]+)\s*;")


def strip_matlab_comments(text: str) -> str:
    """Drop percent comments outside quoted strings while preserving line numbers."""
    cleaned = []
    block_comment = False
    for line in text.splitlines():
        if line.strip() == "%{":
            block_comment = True
        if block_comment:
            if line.strip() == "%}":
                block_comment = False
            cleaned.append("")
            continue
        quoted = False
        i = 0
        while i < len(line):
            if line[i] == "'":
                if quoted and i + 1 < len(line) and line[i + 1] == "'":
                    i += 2
                    continue
                quoted = not quoted
            elif line[i] == "%" and not quoted:
                line = line[:i]
                break
            i += 1
        cleaned.append(line)
    return "\n".join(cleaned)


def active_wbm_lines(text: str) -> list[str]:
    """Filter literal if 0/1 and Recon3D branches for the whole-body protocol.

    This is deliberately a source extractor, not a MATLAB interpreter. Unknown
    runtime conditions are retained; loop/end nesting prevents a disabled IEM's
    inner end from accidentally re-enabling its remaining lines.
    """
    lines = strip_matlab_comments(text).splitlines()
    stack = []
    out = []

    def condition(expr):
        expr = expr.strip().rstrip(";").strip()
        if expr in {"0", "1"}:
            return expr == "1"
        if re.fullmatch(r"~strcmp\(modelName,\s*'Recon3D'\)", expr):
            return True
        if re.fullmatch(r"strcmp\(modelName,\s*'Recon3D'\)", expr):
            return False
        return None

    for line in lines:
        word = re.match(r"^\s*(if|for|parfor|while|switch|try|elseif|else|end)\b(.*)", line)
        before = all(frame[1] for frame in stack)
        if word:
            kind, rest = word.groups()
            if kind in {"if", "for", "parfor", "while", "switch", "try"}:
                known = condition(rest) if kind == "if" else None
                stack.append([known, known is not False])
            elif kind in {"else", "elseif"} and stack:
                previous = stack[-1][0]
                if previous is not None:
                    new = (not previous) if kind == "else" else (not previous and condition(rest) is not False)
                    stack[-1][1] = new
                    if kind == "elseif":
                        stack[-1][0] = previous or new
            elif kind == "end" and stack:
                stack.pop()
        enabled = before and all(frame[1] for frame in stack)
        out.append(line if enabled else "")
    return out


def extract_protocol(text: str) -> list[dict]:
    original = strip_matlab_comments(text).splitlines()
    source_calls = [i for i, line in enumerate(original) if "checkIEM_WBM(" in line]
    lines = active_wbm_lines(text)
    calls = [i for i, line in enumerate(lines) if "checkIEM_WBM(" in line]
    prev = 0
    out = []
    for k, ci in enumerate(calls):
        chunk = lines[prev:ci + 1]
        resets = [i for i, line in enumerate(chunk) if re.match(r"\s*model\s*=\s*modelO\s*;", line)]
        if resets:
            chunk = chunk[resets[-1]:]
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
        # Source blocks also constrain their selected IEM reactions directly,
        # e.g. EP makes _XYLUR irreversible using find(ismember(...)).
        strvars["IEMRxns"] = include
        r2_assignment = re.search(r"^\s*R2\s*=\s*\{(.*?)\}\s*;", txt0, re.S | re.M)
        if r2_assignment and re.search(r"\bX\s*=\s*unique\(RxnsAll2\)", txt0):
            strvars["X"] = re.findall(r"'([^']+)'", r2_assignment.group(1))
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
        out.append({"iem": abbr, "call_index": source_calls.index(ci) + 1, "include_patterns": list(dict.fromkeys(include)),
                    "exclude_patterns": list(dict.fromkeys(exclude)), "bound_tweaks": tweaks, "biomarkers": wbm,
                    "demand_metabolites": list(dict.fromkeys(re.findall(r"addDemandReaction\(model,\s*'([^']+)'", txt))),
                    "source_call_line": ci + 1})
        prev = ci + 1
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", default=SRC)
    ap.add_argument("--output", default=OUT)
    args = ap.parse_args()
    with open(args.source, "rb") as fh:
        source = fh.read()
    out = extract_protocol(source.decode("utf-8"))
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(out, fh, indent=2)
    # Git blob identity remains useful even when the source comes from a local checkout.
    provenance = {
        "source_url": "https://github.com/opencobra/cobratoolbox/blob/master/src/analysis/wholeBody/PSCMToolbox/runIEM_HH.m",
        "source_git_blob_sha1": hashlib.sha1(b"blob " + str(len(source)).encode() + b"\0" + source).hexdigest(),
        "source_sha256": hashlib.sha256(source).hexdigest(),
        "parser_version": "0.2",
        "n_iem_calls": len(out), "n_biomarkers": sum(len(p["biomarkers"]) for p in out),
        "scope": "Whole-body branches; literal disabled blocks and comments excluded. Physiological/diet preprocessing still not ported.",
    }
    provenance["source_blob_api_url"] = "https://api.github.com/repos/opencobra/cobratoolbox/git/blobs/" + provenance["source_git_blob_sha1"]
    with open(os.path.splitext(args.output)[0] + ".provenance.json", "w") as fh:
        json.dump(provenance, fh, indent=2)
    print(json.dumps(provenance, indent=2))
    empty = [o["iem"] for o in out if not o["include_patterns"]]
    if empty:
        raise ValueError(f"No IEM reactions extracted for {empty}")


if __name__ == "__main__":
    main()
