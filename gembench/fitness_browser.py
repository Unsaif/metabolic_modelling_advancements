"""Fitness Browser (Price et al. 2018) RB-TnSeq data for any organism, with carbon-source
conditions mapped to BiGG exchanges and base media mapped to BiGG components.

Files per organism (downloaded 5 September 2026, see data/fitness_browser/PROVENANCE.md):
  experiments.tsv   experiment metadata (expName, expGroup, media, condition_1, units_1, concentration_1, condition_2, ...)
  fit_logratios.tsv gene fitness (log2), columns '<expName> <short description>'
  fit_t.tsv         t-like statistics (optional)
  genes.tsv         locus table (locusId, sysName, ...)

Reference tables (curated here, data/reference/):
  fitness_browser_carbon_sources_bigg.tsv   condition -> BiGG metabolite ids (';' separated) or none
  fitness_browser_media_bigg.tsv            media -> aerobic flag and BiGG components
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .media import Medium, bigg_exchange

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FB_DIR = os.path.join(ROOT, "data", "fitness_browser")
REF_DIR = os.path.join(ROOT, "data", "reference")


@dataclass
class Condition:
    name: str                      # Fitness Browser condition_1 (controlled vocabulary)
    media: str                     # Fitness Browser media name
    bigg_ids: List[str]            # BiGG metabolite ids opened as carbon source(s); empty if unmapped
    experiments: List[str]         # expName values averaged
    concentration: str
    mapping_confidence: str
    note: str = ""

    @property
    def key(self) -> str:
        return f"{self.name} | {self.media}"

    @property
    def exchanges(self) -> List[str]:
        return [bigg_exchange(m) for m in self.bigg_ids]


@dataclass
class FitnessBrowserOrganism:
    org_id: str
    genes: pd.DataFrame
    experiments: pd.DataFrame
    fitness: pd.DataFrame          # index sysName, columns expName
    t_scores: Optional[pd.DataFrame]
    provenance: Dict[str, str] = field(default_factory=dict)


def load_organism(org_id: str, data_dir: str = FB_DIR) -> FitnessBrowserOrganism:
    d = os.path.join(data_dir, org_id)
    genes = pd.read_table(os.path.join(d, "genes.tsv"), dtype=str, keep_default_na=False)
    exp = pd.read_table(os.path.join(d, "experiments.tsv"), dtype=str, keep_default_na=False)
    fit = pd.read_table(os.path.join(d, "fit_logratios.tsv"), dtype=str, keep_default_na=False)
    meta_cols = [c for c in ["orgId", "locusId", "sysName", "geneName", "desc"] if c in fit.columns]
    exp_cols = [c for c in fit.columns if c not in meta_cols]
    colmap = {c.split(" ")[0]: c for c in exp_cols}
    fitness = fit[exp_cols].apply(pd.to_numeric, errors="coerce")
    fitness.columns = [c.split(" ")[0] for c in exp_cols]
    fitness.index = fit["sysName"].where(fit["sysName"] != "", fit["locusId"]).values
    t = None
    tp = os.path.join(d, "fit_t.tsv")
    if os.path.exists(tp):
        tt = pd.read_table(tp, dtype=str, keep_default_na=False)
        tcols = [c for c in tt.columns if c not in meta_cols]
        t = tt[tcols].apply(pd.to_numeric, errors="coerce")
        t.columns = [c.split(" ")[0] for c in tcols]
        t.index = fitness.index
    prov = {"dataset": f"Fitness Browser RB-TnSeq gene fitness, orgId '{org_id}'",
            "primary_source": "Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov",
            "download": "5 September 2026 via createFitData.cgi / createExpData.cgi / orgGenes.cgi (see data/fitness_browser/PROVENANCE.md)",
            "n_genes_with_fitness": str(fitness.shape[0]), "n_experiments": str(fitness.shape[1])}
    return FitnessBrowserOrganism(org_id=org_id, genes=genes, experiments=exp, fitness=fitness, t_scores=t, provenance=prov)


def load_carbon_source_map(path: str = os.path.join(REF_DIR, "fitness_browser_carbon_sources_bigg.tsv")) -> pd.DataFrame:
    return pd.read_table(path, dtype=str, keep_default_na=False)


def load_media_map(path: str = os.path.join(REF_DIR, "fitness_browser_media_bigg.tsv")) -> pd.DataFrame:
    return pd.read_table(path, dtype=str, keep_default_na=False)


TRACE_METALS = ("fe2", "fe3", "mn2", "zn2", "cu2", "cobalt2", "mobd", "ni2", "sel", "slnt", "tungs")


def base_medium(media_name: str, media_map: Optional[pd.DataFrame] = None, uptake: float = -1000.0,
                trace_uptake: float = -0.001, trace_metal_uptake: float = -0.1) -> Medium:
    """Medium of a Fitness Browser experiment as BiGG exchange bounds.

    Bulk components are unlimited; organic trace components (vitamins, hemin, reductant amino acids, nucleobases)
    are limited to `trace_uptake` so that they cannot serve as carbon sources; trace metals are limited to
    `trace_metal_uptake` (0.1 mmol/gDW/h, ten times the biomass demand at a growth rate of 1/h) so that an
    unlimited metal cannot act as an unlimited electron acceptor — with unlimited Fe(III) a draft carrying a
    gene-less extracellular ferric reductase respires iron and, once it has an ATP synthase, grows on the
    proton gradient that produces (decision D17, Sprint 3 note)."""
    mm = media_map if media_map is not None else load_media_map()
    row = mm[mm["media"] == media_name]
    if row.empty:
        raise KeyError(f"no BiGG mapping for medium '{media_name}'")
    row = row.iloc[0]
    comps = [c for c in row["bigg_ids"].split(";") if c]
    aerobic = row["aerobic"].strip().lower() in ("yes", "true", "1")
    if aerobic and "o2" not in comps:
        comps.append("o2")
    trace = set(str(row.get("trace_components", "")).split(";")) - {""}

    def bound(c: str) -> float:
        if c in trace:
            return trace_uptake
        if c in TRACE_METALS:
            return trace_metal_uptake
        return uptake

    return Medium(name=media_name, description=row["note"],
                  uptakes={bigg_exchange(c): bound(c) for c in comps},
                  provenance="Fitness Browser media definition (bitbucket.org/berkeleylab/feba metadata/media, mixes) "
                             "mapped to BiGG ids in data/reference/fitness_browser_media_bigg.tsv",
                  notes=[f"aerobic={aerobic}", f"trace organics {trace_uptake}, trace metals {trace_metal_uptake} mmol/gDW/h"])


def carbon_source_conditions(org: FitnessBrowserOrganism, cs_map: Optional[pd.DataFrame] = None,
                             allowed_condition_2: tuple = ("", "Dimethyl Sulfoxide"),
                             media_map: Optional[pd.DataFrame] = None) -> List[Condition]:
    """Carbon-source experiments grouped into conditions (condition_1 x media).

    Experiments with a second condition are excluded unless it is empty or DMSO (the solvent used
    for hydrophobic carbon sources). Conditions whose base medium has no BiGG mapping are skipped.
    """
    cm = cs_map if cs_map is not None else load_carbon_source_map()
    cm = cm.set_index("condition")
    mm = media_map if media_map is not None else load_media_map()
    known_media = set(mm["media"])
    exp = org.experiments
    sel = exp[(exp["expGroup"] == "carbon source") & (exp["condition_2"].isin(allowed_condition_2))]
    out: Dict[str, Condition] = {}
    for _, r in sel.iterrows():
        name, media = r["condition_1"], r["media"]
        if media not in known_media:
            continue
        if r["expName"] not in org.fitness.columns:
            continue
        key = f"{name} | {media}"
        if key not in out:
            if name in cm.index:
                row = cm.loc[name]
                bigg = [x for x in row["bigg_ids"].split(";") if x]
                conf, note = row["confidence"], row["note"]
            else:
                bigg, conf, note = [], "none", "condition not in the curated mapping table"
            out[key] = Condition(name=name, media=media, bigg_ids=bigg, experiments=[],
                                 concentration=f"{r['concentration_1']} {r['units_1']}".strip(),
                                 mapping_confidence=conf, note=note)
        out[key].experiments.append(r["expName"])
    return list(out.values())


def replicate_averaged_fitness(org: FitnessBrowserOrganism, conditions: List[Condition]) -> pd.DataFrame:
    """genes x conditions matrix of replicate-averaged log2 fitness (NaN-aware mean)."""
    cols = {}
    for c in conditions:
        cols[c.key] = org.fitness[c.experiments].mean(axis=1, skipna=True)
    return pd.DataFrame(cols, index=org.fitness.index)
