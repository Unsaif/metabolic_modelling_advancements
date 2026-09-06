"""Offline SBML screens; these detect model representations, not organism biology."""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from collections import Counter


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _metabolite(identifier: str) -> tuple[str, str]:
    identifier = identifier.removeprefix("M_")
    base, sep, compartment = identifier.rpartition("_")
    return (base, compartment) if sep else (identifier, "")


def _ion_coupled_atp_candidate(reaction: ET.Element) -> bool:
    """Recognize BiGG-style ATP/ADP/Pi plus transmembrane H or Na equations.

    Exclude ATP-powered cotransport of other solutes. This is a candidate screen:
    it neither checks usable flux bounds nor establishes physiological direction.
    """
    stoich: dict[tuple[str, str], float] = {}
    for side in reaction:
        sign = {"listOfReactants": -1, "listOfProducts": 1}.get(_local(side.tag))
        if sign is None:
            continue
        for ref in side:
            if _local(ref.tag) != "speciesReference":
                continue
            met = _metabolite(ref.attrib["species"])
            stoich[met] = stoich.get(met, 0.0) + sign * float(ref.get("stoichiometry", "1"))
    stoich = {met: coefficient for met, coefficient in stoich.items() if coefficient}
    if any(met not in {"atp", "adp", "pi", "h2o", "h", "na1"} for met, _ in stoich):
        return False
    for (base, inside), atp in stoich.items():
        if base != "atp" or not inside:
            continue
        if any(abs(stoich.get((met, inside), 0.0) + atp) > 1e-8 for met in ("adp", "pi")):
            continue
        for (ion, outside), coefficient in stoich.items():
            if ion not in {"h", "na1"} or outside == inside or not outside:
                continue
            # Orient to ATP synthesis; ions must enter rather than leave.
            if coefficient / atp >= 0:
                continue
            inside_coefficient = stoich.get((ion, inside), 0.0)
            expected = -coefficient - atp if ion == "h" else -coefficient
            if abs(inside_coefficient - expected) < 1e-8:
                return True
    return False


def screen_sbml(xml: bytes) -> dict:
    """Inspect actual reaction elements, never arbitrary IDs in notes or genes."""
    root = ET.fromstring(xml)
    if _local(root.tag) != "sbml":
        raise ValueError("Downloaded document is not SBML")
    reactions = [element for element in root.iter() if _local(element.tag) == "reaction"]
    if not reactions:
        raise ValueError("SBML document contains no reactions")
    ids, names, equations = [], [], []
    for reaction in reactions:
        identifier = reaction.attrib["id"]
        if re.match(r"^(?:R_)?ATPS", identifier, re.IGNORECASE):
            ids.append(identifier)
        if re.search(r"\batp\s+synth(?:ase|etase)\b", reaction.get("name", ""), re.IGNORECASE):
            names.append(identifier)
        if _ion_coupled_atp_candidate(reaction):
            equations.append(identifier)
    return {
        "atps_ids": ";".join(sorted(set(ids))),
        "atp_synthase_name_candidates": ";".join(sorted(set(names))),
        "ion_coupled_atp_candidates": ";".join(sorted(set(equations))),
        "n_reactions": len(reactions),
        "has_atps": bool(ids),
    }


def summarize_rows(rows: list[dict]) -> dict:
    """Summarize new surveys and both archived TSV schemas, excluding failures."""
    present = absent = failed = 0
    identifiers = []
    for row in rows:
        field = "atps_ids" if "atps_ids" in row else "atp_synthase"
        if field not in row:
            raise ValueError("Survey has no atps_ids or atp_synthase column")
        value = str(row.get(field) or "")
        invalid = row.get("status", "ok") != "ok" or value.startswith("DOWNLOAD_FAILED")
        if "bytes" in row:
            invalid |= int(row.get("bytes") or 0) <= 0
        count = row.get("n_reactions", row.get("n_rxns"))
        if count is not None:
            invalid |= int(count or 0) <= 0
        if invalid:
            failed += 1
        elif any(part.strip() for part in value.split(";")):
            present += 1
        else:
            absent += 1
        identifiers.append(row.get("assembly") or row.get("model") or "")
    successes = present + absent
    return {
        "rows": len(rows), "successful_models": successes, "failed_models": failed,
        "with_atps_identifier": present, "without_atps_identifier": absent,
        "fraction_without_atps_identifier": absent / successes if successes else None,
        "duplicate_model_identifiers": sorted(k for k, n in Counter(identifiers).items() if k and n > 1),
        "interpretation": "Absence of an ATPS reaction identifier in the sampled model files; not biological absence, current CarveMe performance, or a demonstrated reconstruction failure mechanism.",
    }
