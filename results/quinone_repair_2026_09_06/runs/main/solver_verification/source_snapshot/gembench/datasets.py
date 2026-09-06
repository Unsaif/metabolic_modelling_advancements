"""Ground-truth datasets with provenance.

Each loader returns a plain dataclass holding the data plus a `provenance`
dictionary that is copied verbatim into benchmark cards.
"""
from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXTERNAL = os.path.join(ROOT, "external")


def _git_head(path: str) -> Optional[str]:
    try:
        return subprocess.check_output(["git", "-C", path, "rev-parse", "HEAD"],
                                       text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return None


# ---------------------------------------------------------------------------
# E. coli BW25113 RB-TnSeq fitness across carbon sources (Fitness Browser)
# ---------------------------------------------------------------------------

@dataclass
class CarbonFitnessDataset:
    genes: List[str]                       # b-numbers (sysName)
    carbon_sources: List[str]              # BiGG metabolite ids (unique, sorted)
    experiments: Dict[str, List[str]]      # carbon source -> Fitness Browser experiment names
    fitness: np.ndarray                    # genes x carbon sources, replicate-averaged
    fitness_per_experiment: pd.DataFrame   # genes x experiments (raw)
    provenance: Dict[str, str] = field(default_factory=dict)


def load_ecoli_bw25113_carbon_fitness(
    repo_dir: str = os.path.join(EXTERNAL, "E_coli_GEM_validation"),
    selection: str = "bernstein2023",
) -> CarbonFitnessDataset:
    """Load the Fitness Browser E. coli BW25113 ('Keio') carbon-source fitness data.

    selection:
      'bernstein2023' — the 27 single-component carbon sources on M9 minimal medium
         that Bernstein et al. selected by row index (D,L-lactate and casamino acids
         excluded because they are not single BiGG components; MOPS-medium experiments
         excluded). Sucrose and 'man' are removed later by the protocol, not here.
      'all_m9_single' — every 'carbon source' experiment on 'M9 minimal media_noCarbon'
         with a single BiGG component mapping (superset check of the above).
    """
    fdir = os.path.join(repo_dir, "Fitness_Data", "E_coli_BW25113")
    exp = pd.read_table(os.path.join(fdir, "exp_organism_Keio_Mapped.txt"), encoding="latin-1")
    fit = pd.read_table(os.path.join(fdir, "fit_organism_Keio.tsv"))

    mask = (exp["expGroup"] == "carbon source") & (exp["media"] == "M9 minimal media_noCarbon")
    comp = exp["BiGG Component"].astype(str)
    single = comp.notna() & (comp != "nan") & (~comp.str.contains(";"))
    sel = exp[mask & single].copy()
    if selection == "bernstein2023":
        # Bernstein et al. read rows 1-28, 31-50, 53-58 of the metadata file (0 = header).
        idx = list(range(0, 28)) + list(range(30, 50)) + list(range(52, 58))
        sel = exp.iloc[idx].copy()
        sel = sel[sel["BiGG Component"].notna() & (~sel["BiGG Component"].astype(str).str.contains(";"))]
    elif selection != "all_m9_single":
        raise ValueError(selection)

    # fitness columns are named '<expName> <expDesc>'
    colmap = {c.split(" ")[0]: c for c in fit.columns[5:]}
    genes = fit["sysName"].astype(str).tolist()
    experiments: Dict[str, List[str]] = {}
    for _, row in sel.iterrows():
        cs = str(row["BiGG Component"]).strip()
        experiments.setdefault(cs, []).append(row["expName"])
    carbon_sources = sorted(experiments)
    per_exp_cols = [colmap[e] for cs in carbon_sources for e in experiments[cs] if e in colmap]
    per_exp = fit[["sysName"] + per_exp_cols].set_index("sysName")
    mat = np.zeros((len(genes), len(carbon_sources)))
    for j, cs in enumerate(carbon_sources):
        cols = [colmap[e] for e in experiments[cs] if e in colmap]
        mat[:, j] = fit[cols].to_numpy(dtype=float).mean(axis=1)

    prov = {
        "dataset": "Fitness Browser RB-TnSeq gene fitness, E. coli BW25113 (orgId 'Keio')",
        "primary_source": "Price et al. 2018, Nature 557:503-509, https://fit.genomics.lbl.gov",
        "copy_obtained_from": "github.com/dbernste/E_coli_GEM_validation (MIT-licensed code; data redistributed there)",
        "copy_commit": _git_head(repo_dir) or "unknown",
        "selection": selection,
        "n_genes": str(len(genes)),
        "n_carbon_sources": str(len(carbon_sources)),
        "n_experiments": str(sum(len(v) for v in experiments.values())),
        "units": "log2 gene fitness (Fitness Browser convention); replicates averaged per carbon source",
        "licence_note": "Fitness Browser data are publicly available; check fit.genomics.lbl.gov terms before redistribution.",
    }
    return CarbonFitnessDataset(genes=genes, carbon_sources=carbon_sources, experiments=experiments,
                                fitness=mat, fitness_per_experiment=per_exp, provenance=prov)


# ---------------------------------------------------------------------------
# S. cerevisiae deletion viability (Stanford deletion collection via yeast-GEM)
# ---------------------------------------------------------------------------

@dataclass
class YeastEssentialityDataset:
    inviable_orfs: set
    verified_orfs: set
    provenance: Dict[str, str] = field(default_factory=dict)


def load_yeast_essentiality(repo_dir: str = os.path.join(EXTERNAL, "yeast-GEM")) -> YeastEssentialityDataset:
    d = os.path.join(repo_dir, "data", "essentialGenes")
    inv = {l.strip() for l in open(os.path.join(d, "inviable_orfs.txt")) if l.strip()}
    ver = {l.strip() for l in open(os.path.join(d, "verified_orfs.txt")) if l.strip()}
    ver_txt = open(os.path.join(repo_dir, "version.txt")).read().strip() if os.path.exists(os.path.join(repo_dir, "version.txt")) else "unknown"
    prov = {
        "dataset": "Stanford yeast deletion project inviable ORFs (14 Aug 2011 snapshot) and SGD verified ORFs (27 Aug 2013)",
        "copy_obtained_from": f"github.com/SysBioChalmers/yeast-GEM data/essentialGenes (yeast-GEM version {ver_txt})",
        "copy_commit": _git_head(repo_dir) or "unknown",
        "n_inviable_unique": str(len(inv)),
        "n_verified": str(len(ver)),
        "caveat": "Screened in complex media supplemented with auxotrophic markers; the yeast-GEM test compares "
                  "against a synthetic complete medium, so the reference is imperfect but stable (yeast-GEM README).",
    }
    return YeastEssentialityDataset(inviable_orfs=inv, verified_orfs=ver, provenance=prov)
