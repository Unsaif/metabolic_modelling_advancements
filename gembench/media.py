"""Explicit media definitions and their application to COBRA models.

A Medium is a named mapping from exchange-reaction identifiers to lower bounds
(negative = uptake allowed).  Applying a medium first closes every exchange
(lower bound 0, upper bound +1000) and then opens the listed uptakes, so that
the state of the model after `apply_medium` depends only on the medium and
not on whatever bounds the SBML file shipped with.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import cobra


@dataclass
class Medium:
    name: str
    description: str
    uptakes: Dict[str, float]           # exchange id -> lower bound (e.g. -1000)
    provenance: str = ""
    notes: List[str] = field(default_factory=list)

    def with_carbon_source(self, exchange_id: str, lower_bound: float = -10.0,
                           name: Optional[str] = None) -> "Medium":
        upt = dict(self.uptakes)
        upt[exchange_id] = lower_bound
        return Medium(name=name or f"{self.name}+{exchange_id}",
                      description=f"{self.description}; carbon source {exchange_id} at {lower_bound}",
                      uptakes=upt, provenance=self.provenance, notes=list(self.notes))


def bigg_exchange(component: str) -> str:
    """BiGG convention: metabolite 'glc__D' -> exchange 'EX_glc__D_e'."""
    return f"EX_{component}_e"


# Base M9 minimal medium (no carbon) as encoded for BiGG iML1515 in the
# Fitness Browser experiment mapping used by Bernstein et al. 2023
# (Fitness_Data/E_coli_BW25113/exp_organism_Keio_Mapped_Media.txt).
_M9_BIGG_COMPONENTS = [
    "pi", "co2", "fe3", "h", "mn2", "fe2", "zn2", "mg2", "ca2", "ni2", "cu2",
    "sel", "cobalt2", "h2o", "mobd", "so4", "nh4", "k", "na1", "cl", "o2",
    "tungs", "slnt",
]

M9_MINIMAL_NO_CARBON = Medium(
    name="M9_minimal_noCarbon_BiGG",
    description="M9 minimal medium without a carbon source; all listed non-carbon "
                "components unlimited (-1000 mmol/gDW/h); aerobic.",
    uptakes={bigg_exchange(c): -1000.0 for c in _M9_BIGG_COMPONENTS},
    provenance="Fitness Browser 'M9 minimal media_noCarbon' mapped to BiGG iML1515 "
               "exchange ids by Bernstein et al. 2023 (github.com/dbernste/E_coli_GEM_validation).",
    notes=["Carbon sources are added separately at -10 mmol/gDW/h (Bernstein et al. 2023)."],
)


def apply_medium(model: cobra.Model, medium: Medium, close_all: bool = True) -> List[str]:
    """Apply `medium` to `model` in place. Returns exchange ids that were not found."""
    if close_all:
        for ex in model.exchanges:
            ex.lower_bound = 0.0
            ex.upper_bound = 1000.0
    missing: List[str] = []
    for ex_id, lb in medium.uptakes.items():
        try:
            model.reactions.get_by_id(ex_id).lower_bound = lb
        except KeyError:
            missing.append(ex_id)
    return missing
