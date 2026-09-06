"""Structural sanity checks used as verifier gates for patches.

`energy_from_nothing` closes every boundary reaction and maximises ATP hydrolysis (and, optionally, NADH/NADPH oxidation
and proton-gradient dissipation): any positive value is an energy-generating cycle (Fritzemeier et al. 2017), the
classic way a new reaction — an ATP synthase, a transporter, a reversible step — can silently break a model.
"""
from __future__ import annotations

from typing import Dict
import math

import cobra

_CURRENCIES = {
    "atp": ({"atp_c": -1.0, "h2o_c": -1.0, "adp_c": 1.0, "pi_c": 1.0, "h_c": 1.0}, "ATP hydrolysis"),
    "nadh": ({"nadh_c": -1.0, "nad_c": 1.0, "h_c": 1.0}, "NADH oxidation"),
    "nadph": ({"nadph_c": -1.0, "nadp_c": 1.0, "h_c": 1.0}, "NADPH oxidation"),
    "q8h2": ({"q8h2_c": -1.0, "q8_c": 1.0, "h_c": 2.0}, "ubiquinol oxidation"),
    "h_p": ({"h_p": -1.0, "h_c": 1.0}, "proton gradient dissipation"),
}


def energy_from_nothing(model: cobra.Model, currencies=("atp", "nadh", "nadph", "q8h2", "h_p"), tol: float = 1e-6) -> Dict[str, float]:
    """Maximal dissipation with all exchanges, demands and sinks closed.

    Mandatory fluxes are relaxed to include zero while retaining their allowed
    direction: maintenance or biomass requirements must not make a closed model
    infeasible and thereby conceal a cycle. Missing currencies are omitted, not
    certified as absent. A failed solve raises rather than certifying zero flux.
    Returns {currency: flux}; values above `tol` indicate an energy-generating cycle.
    """
    out: Dict[str, float] = {}
    with model:
        for reaction in model.reactions:
            if reaction.boundary:
                reaction.bounds = (0.0, 0.0)
            else:
                reaction.bounds = (min(0.0, reaction.lower_bound), max(0.0, reaction.upper_bound))
        for key in currencies:
            stoich, name = _CURRENCIES[key]
            if any(m not in model.metabolites for m in stoich):
                continue
            with model:
                rid = f"EGC_{key}"
                while rid in model.reactions:
                    rid += "_test"
                dm = cobra.Reaction(rid, name=name, lower_bound=0.0, upper_bound=1000.0)
                dm.add_metabolites({model.metabolites.get_by_id(m): c for m, c in stoich.items()})
                model.add_reactions([dm])
                model.objective = dm
                model.objective_direction = "max"
                v = model.slim_optimize()
                if model.solver.status != "optimal" or v is None or not math.isfinite(v):
                    raise RuntimeError(f"Energy-cycle check for {key} failed: solver status {model.solver.status}")
                out[key] = float(v)
    return out
