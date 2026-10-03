"""Fitness Browser (FEBA) medium recipes -> BiGG medium rows, by a fixed rule (transfer study v1).

The FEBA code repository (bitbucket.org/berkeleylab/feba) describes every medium in metadata/media and every stock
mix in metadata/mixes, as blocks separated by blank lines:

    Media<TAB>name
    Description<TAB>...
    Minimal<TAB>TRUE|FALSE          (media)   or   X<TAB>100   (mixes: an N-fold stock, used in media at "1 X")
    Controlled vocabulary<TAB>Concentration<TAB>Units
    <component><TAB><concentration><TAB><units>
    ...

`medium_row` turns one medium into a row of the media table read by gembench.fitness_browser.base_medium
(media, aerobic, bigg_ids, note, trace_components) using a component dictionary
(data/studies/transfer_v1/feba_component_bigg.tsv) whose classes are

    inorganic      salts, acids and bases: mapped to their ions (counter-ions and water of hydration only as ions)
    trace_organic  vitamins, cofactors, hemin, nucleobases, reductants: mapped and limited to trace uptake
    organic        any other defined organic compound: mapped and unlimited (it can serve as a carbon source)
    ignore         buffers and chelators without a BiGG metabolite (MOPS, PIPES, HEPES, tricine, NTA, EDTA, ...)
    undefined      complex ingredients (yeast extract, peptone, tryptone, casamino acids, ...): the medium is
                   undefined and is not mapped, so its experiments are excluded

Fixed additions (the convention of the development media table): h2o, h and co2 always; the full trace-metal set
(gembench.fitness_browser.TRACE_METALS, capped at trace-metal uptake by base_medium) always; o2 when the carbon-source
experiments on the medium are recorded as aerobic.
"""
from __future__ import annotations

import os
from collections import Counter
from typing import Dict, List, Optional, Tuple

import pandas as pd

from .fitness_browser import TRACE_METALS

ALWAYS = ("h2o", "h", "co2")


def parse_blocks(path: str) -> Dict[str, dict]:
    """Parse a FEBA media or mixes file into {name: {description, minimal, x, components: [(name, conc, units)]}}."""
    out: Dict[str, dict] = {}
    cur: Optional[dict] = None
    in_components = False
    for raw in open(path, encoding="utf-8"):
        cells = [c.strip() for c in raw.rstrip("\n").split("\t")]
        if not any(cells):
            cur, in_components = None, False
            continue
        key = cells[0]
        if key == "Media" and len(cells) > 1 and cells[1]:
            cur = {"name": cells[1], "description": "", "minimal": None, "x": None, "components": []}
            out[cells[1]] = cur
            in_components = False
        elif cur is None:
            continue
        elif key == "Description" and not in_components:
            cur["description"] = cells[1] if len(cells) > 1 else ""
        elif key == "Minimal" and not in_components:
            cur["minimal"] = (cells[1] if len(cells) > 1 else "").upper() == "TRUE"
        elif key == "X" and not in_components:
            cur["x"] = cells[1] if len(cells) > 1 else ""
        elif key == "Controlled vocabulary":
            in_components = True
        elif in_components:
            conc = cells[1] if len(cells) > 1 else ""
            units = cells[2] if len(cells) > 2 else ""
            cur["components"].append((key, conc, units))
    return out


def expand(name: str, media: Dict[str, dict], mixes: Dict[str, dict], _depth: int = 0) -> List[Tuple[str, str, str, str]]:
    """Components of a medium with mixes (units 'X') expanded recursively: (component, conc, units, via)."""
    if _depth > 5:
        raise ValueError(f"mix nesting too deep at {name}")
    block = media.get(name) or mixes.get(name)
    if block is None:
        raise KeyError(f"medium or mix '{name}' not found")
    rows = []
    for comp, conc, units in block["components"]:
        if units.upper() == "X" and (comp in mixes or comp in media):
            for c, cc, u, via in expand(comp, media, mixes, _depth + 1):
                rows.append((c, cc, u, comp if via == "" else f"{comp} > {via}"))
        else:
            rows.append((comp, conc, units, ""))
    return rows


def load_component_table(paths) -> pd.DataFrame:
    """One or more dictionary files (later files may add components, never redefine one)."""
    paths = [paths] if isinstance(paths, str) else list(paths)
    t = pd.concat([pd.read_table(p, dtype=str, keep_default_na=False, comment="#") for p in paths], ignore_index=True)
    t["key"] = t["component"].str.strip().str.lower()
    dup = t["key"][t["key"].duplicated()]
    if len(dup):
        raise ValueError(f"components defined more than once: {sorted(set(dup))}")
    return t.set_index("key")


def aerobic_flag(experiments: pd.DataFrame, media_name: str) -> Optional[bool]:
    """Majority value of the 'aerobic' column over the carbon-source experiments on this medium (None if absent)."""
    if "aerobic" not in experiments.columns:
        return None
    sel = experiments[(experiments["media"] == media_name) & (experiments["expGroup"] == "carbon source")]
    vals = Counter(v.strip().lower() for v in sel["aerobic"] if v.strip())
    if not vals:
        return None
    top = vals.most_common(1)[0][0]
    return top == "aerobic"


def medium_row(name: str, media: Dict[str, dict], mixes: Dict[str, dict], components: pd.DataFrame,
               aerobic: Optional[bool]) -> dict:
    """Map one FEBA medium to a media-table row, or report why it cannot be mapped.

    Returns {"media", "status" ('mapped' | 'undefined' | 'unmapped_components' | 'not_found' | 'aerobic_unknown'),
    "aerobic", "bigg_ids", "trace_components", "note", "unmapped", "undefined", "ignored"}."""
    out = {"media": name, "status": "mapped", "aerobic": "yes" if aerobic else "no", "bigg_ids": "",
           "trace_components": "", "note": "", "unmapped": [], "undefined": [], "ignored": []}
    try:
        rows = expand(name, media, mixes)
    except KeyError:
        out["status"] = "not_found"
        return out
    if aerobic is None:
        out["status"] = "aerobic_unknown"
    ids: List[str] = list(ALWAYS)
    trace_moiety, bulk = set(), set(ALWAYS)
    for comp, conc, units, via in rows:
        key = comp.strip().lower()
        if key not in components.index:
            out["unmapped"].append(comp)
            continue
        r = components.loc[key]
        cls = r["class"]
        if cls == "undefined":
            out["undefined"].append(comp)
            continue
        own = [b for b in r["bigg_ids"].split(";") if b]
        ions = [b for b in r.get("counter_ions", "").split(";") if b]
        if cls == "ignore":
            out["ignored"].append(comp)          # buffer or chelator: only its counter-ions are added
            own = []
        for b in own + ions:
            if b not in ids:
                ids.append(b)
        # a trace organic's own metabolites are limited to trace uptake unless a bulk component also supplies them
        (trace_moiety if cls == "trace_organic" else bulk).update(own)
        bulk.update(ions)
    trace = [b for b in ids if b in trace_moiety and b not in bulk]
    for m in TRACE_METALS:
        if m not in ids:
            ids.append(m)
    if out["undefined"]:
        out["status"] = "undefined"
    elif out["unmapped"]:
        out["status"] = "unmapped_components"
    out["bigg_ids"] = ";".join(ids)
    out["trace_components"] = ";".join(trace)
    desc = (media.get(name) or {}).get("description", "")
    out["note"] = (f"FEBA recipe '{name}' ({desc}) mapped by gembench.feba_media with "
                   f"data/studies/transfer_v1/feba_component_bigg.tsv; ignored (buffers/chelators): "
                   f"{', '.join(sorted(set(out['ignored']))) or 'none'}; trace metals always present (capped)")
    return out


def media_table_rows(names: List[str], media_path: str, mixes_path: str, component_paths,
                     experiments: Dict[str, pd.DataFrame]) -> List[dict]:
    """Rows for several media; `experiments` maps orgId -> experiment metadata (for the aerobic flag)."""
    media, mixes = parse_blocks(media_path), parse_blocks(mixes_path)
    comps = load_component_table(component_paths)
    rows = []
    for n in names:
        flags = {aerobic_flag(e, n) for e in experiments.values()} - {None}
        aerobic = flags.pop() if len(flags) == 1 else None
        rows.append(medium_row(n, media, mixes, comps, aerobic))
    return rows
