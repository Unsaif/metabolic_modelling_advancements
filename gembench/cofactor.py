"""Explicit cofactor-pool certificates and production probes for ordinary FBA.

Conserved pools and demand maximization are established methods. These helpers
make the declared pool, exact row sum, model constraints and failed solves
reviewable. They neither discover chemical moieties nor certify biological
essentiality, native cofactor composition, or physiological dilution rates.
"""
from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
import math

import cobra
import numpy as np


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be a finite positive number')
    return float(value)


def _metabolites(model, identifiers):
    if not identifiers or len(set(identifiers)) != len(identifiers):
        raise ValueError('Declare a nonempty set of distinct metabolite identifiers')
    if any(not isinstance(mid, str) or not mid for mid in identifiers):
        raise ValueError('Metabolite identifiers must be nonempty strings')
    missing = [mid for mid in identifiers if mid not in model.metabolites]
    if missing:
        raise ValueError(f'Metabolites absent from the model; no result can be inferred: {missing}')
    return [model.metabolites.get_by_id(mid) for mid in identifiers]


def _assert_stoichiometric_balance(met):
    """Reject custom solver rows that do not express the stored S row."""
    if met.constraint.lb != 0 or met.constraint.ub != 0:
        raise ValueError(f'{met.id} does not have a zero steady-state balance')
    expected = {}
    for reaction in met.reactions:
        coefficient = float(reaction.metabolites[met])
        if not math.isfinite(coefficient):
            raise ValueError(f'Nonfinite stoichiometry in {reaction.id}/{met.id}')
        if coefficient:
            expected[reaction.forward_variable] = coefficient
            expected[reaction.reverse_variable] = -coefficient
    actual = {term: float(coefficient)
              for term, coefficient in met.constraint.expression.as_coefficients_dict().items()
              if coefficient != 0}
    if actual != expected:
        raise ValueError(f'{met.id} solver balance differs from stored stoichiometry; '
                         'a pool certificate for the operative constraints cannot be inferred')


def pool_balance_certificate(model: cobra.Model, weights: Mapping[str, float]):
    """Return the exact weighted stoichiometric row sum for a declared pool.

    Exactness refers to rational arithmetic on the stored floating-point numbers,
    not to chemically verified stoichiometry. The declared metabolite solver
    rows must match their stored stoichiometry. A nonzero coefficient is never
    rounded away. Direction-compatible supply reactions are possibilities, not
    claims that their substrates are available or that they can carry flux.
    """
    if not isinstance(weights, Mapping):
        raise ValueError('Pool weights must map metabolite identifiers to positive numbers')
    metabolites = _metabolites(model, list(weights))
    normalized = {mid: _positive(weight, f'Weight of {mid}') for mid, weight in weights.items()}
    for met in metabolites:
        _assert_stoichiometric_balance(met)
    incident = sorted(set.union(*(set(m.reactions) for m in metabolites)), key=lambda r: r.id)
    rows = []
    for reaction in incident:
        total = Fraction(0)
        terms = {}
        for met in metabolites:
            coefficient = float(reaction.metabolites.get(met, 0.0))
            if not math.isfinite(coefficient):
                raise ValueError(f'Nonfinite stoichiometry in {reaction.id}/{met.id}')
            if coefficient:
                terms[met.id] = coefficient
                total += Fraction.from_float(normalized[met.id]) * Fraction.from_float(coefficient)
        rows.append({'reaction': reaction.id, 'pool_terms': terms,
                     'net_pool_coefficient': float(total),
                     'exact_numerator': str(total.numerator), 'exact_denominator': str(total.denominator),
                     'bounds': list(reaction.bounds), 'boundary': bool(reaction.boundary),
                     'can_supply_by_direction': bool((total > 0 and reaction.upper_bound > 0) or
                                                     (total < 0 and reaction.lower_bound < 0)),
                     'can_consume_by_direction': bool((total < 0 and reaction.upper_bound > 0) or
                                                      (total > 0 and reaction.lower_bound < 0))})
    conserved = all(row['exact_numerator'] == '0' for row in rows)
    return {'pool_weights': normalized, 'exactly_conserved_in_stoichiometry': conserved,
            'n_incident_reactions': len(rows), 'reactions': rows,
            'direction_compatible_supply_reactions': [r['reaction'] for r in rows if r['can_supply_by_direction']],
            'scope': 'Declared metabolite pool and ordinary zero-balance FBA rows; no automatic moiety discovery or biological validation.'}


def probe_metabolite_production(model: cobra.Model, identifiers, *, capacity=1000.0, threshold=1e-6):
    """Maximize separate outward demands under the current model constraints.

    Medium, maintenance, growth bounds, and any additional constraints remain in
    force. Each demand is removed before probing the next metabolite. An optimal
    zero denotes blocked production under these constraints, not biological
    absence. Nonoptimal/nonfinite solutions raise instead of returning zero.
    The artificial demand is a measurement device, never a source or a repair.
    """
    capacity = _positive(capacity, 'Demand capacity')
    threshold = _positive(threshold, 'Production threshold')
    if threshold > capacity:
        raise ValueError('Production threshold exceeds the demand capacity')
    metabolites = _metabolites(model, list(identifiers))
    results = []
    for met in metabolites:
        with model:
            rid = f'COFACTOR_PROBE_{met.id}'
            while rid in model.reactions:
                rid += '_temporary'
            demand = cobra.Reaction(rid, lower_bound=0, upper_bound=capacity)
            demand.add_metabolites({met: -1})
            model.add_reactions([demand])
            model.objective = demand
            model.objective_direction = 'max'
            solution = model.optimize()
            if solution.status != 'optimal' or solution.objective_value is None or not math.isfinite(solution.objective_value):
                raise RuntimeError(f'Production probe for {met.id} failed: {solution.status}')
            values = solution.fluxes.reindex([r.id for r in model.reactions]).to_numpy()
            if not np.isfinite(values).all():
                raise RuntimeError(f'Production probe for {met.id} returned nonfinite fluxes')
            rate = float(solution.objective_value)
            results.append({'metabolite': met.id, 'status': solution.status, 'maximum_demand_flux': rate,
                            'producible_at_threshold': bool(rate >= threshold), 'threshold': threshold,
                            'probe_capacity': capacity,
                            'capacity_reached': bool(rate >= capacity or math.isclose(rate, capacity, rel_tol=1e-9, abs_tol=0.0)),
                            'solver': model.solver.interface.__name__})
    return results
