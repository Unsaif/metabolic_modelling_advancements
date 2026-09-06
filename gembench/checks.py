"""Structural sanity checks used as verifier gates for patches.

`energy_from_nothing` closes every exchange and maximises ATP hydrolysis (and, optionally, NADH/NADPH oxidation
and proton-gradient dissipation): any positive value is an energy-generating cycle (Fritzemeier et al. 2017), the
classic way a new reaction — an ATP synthase, a transporter, a reversible step — can silently break a model.
"""
from __future__ import annotations

from typing import Dict

import cobra

_CURRENCIES = {
    "atp": ({"atp_c": -1.0, "h2o_c": -1.0, "adp_c": 1.0, "pi_c": 1.0, "h_c": 1.0}, "ATP hydrolysis"),
    "nadh": ({"nadh_c": -1.0, "nad_c": 1.0, "h_c": 1.0}, "NADH oxidation"),
    "nadph": ({"nadph_c": -1.0, "nadp_c": 1.0, "h_c": 1.0}, "NADPH oxidation"),
    "q8h2": ({"q8h2_c": -1.0, "q8_c": 1.0, "h_c": 2.0}, "ubiquinol oxidation"),
    "h_p": ({"h_p": -1.0, "h_c": 1.0}, "proton gradient dissipation"),
}


def energy_from_nothing(model: cobra.Model, currencies=("atp", "nadh", "nadph", "q8h2", "h_p"), tol: float = 1e-6) -> Dict[str, float]:
    """Maximal production of each energy currency with every exchange closed (uptake and secretion).
    Returns {currency: flux}; values above `tol` indicate an energy-generating cycle."""
    out: Dict[str, float] = {}
    with model:
        for ex in model.exchanges:
            ex.bounds = (0.0, 0.0)
        for rid in ("ATPM",):
            if rid in model.reactions:
                model.reactions.get_by_id(rid).bounds = (0.0, 1000.0)
        for key in currencies:
            stoich, name = _CURRENCIES[key]
            if any(m not in model.metabolites for m in stoich):
                continue
            with model:
                dm = cobra.Reaction(f"EGC_{key}", name=name, lower_bound=0.0, upper_bound=1000.0)
                dm.add_metabolites({model.metabolites.get_by_id(m): c for m, c in stoich.items()})
                model.add_reactions([dm])
                model.objective = dm
                v = model.slim_optimize()
                out[key] = 0.0 if v is None or v != v else float(v)
    return out
