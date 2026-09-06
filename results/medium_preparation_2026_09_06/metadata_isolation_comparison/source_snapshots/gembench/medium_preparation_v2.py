"""Version 2: explicit medium preparation with declared exchange identity bindings.

This version leaves historical preparation code unchanged. It accepts only
one-metabolite, coefficient -1 exchanges; other orientations need an explicit
adapter. Bounds are flux capacities, not nutrient concentrations. Completion
is the same biological assumption as before: declared components may receive
gene-less, uptake-only carriers, except explicitly excluded components.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import math
from typing import Iterable

import cobra

from .media import Medium


@dataclass(frozen=True)
class CompartmentPolicy:
    external: str
    cytoplasm: str
    external_aliases: tuple[str, ...] = ()
    cytoplasm_aliases: tuple[str, ...] = ()
    completion_exclude: tuple[str, ...] = ('pnto__R', 'fol', 'hco3')
    secretion_capacity: float = 1000.
    exchange_metabolites: dict[str, str] = field(default_factory=dict)


@dataclass
class PreparedMedium:
    model: cobra.Model
    report: dict


def _require(test, message):
    if not test:
        raise ValueError(message)


def _finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _component(exchange):
    _require(isinstance(exchange, str) and exchange.startswith('EX_') and exchange.endswith('_e')
             and len(exchange) > 5, f'Expected canonical EX_<component>_e identifier: {exchange}')
    return exchange[3:-2]


def _validate_medium(medium):
    for rid, lower in medium.uptakes.items():
        _component(rid)
        _require(_finite(lower) and lower <= 0, f'Expected finite nonpositive uptake bound: {rid}')


def _sbo(reaction):
    values = reaction.annotation.get('sbo', '')
    values = [values] if isinstance(values, str) else values
    _require(isinstance(values, (list, tuple)) and all(isinstance(v, str) for v in values),
             f'Invalid SBO annotation: {reaction.id}')
    return {value.upper() for value in values}


def _exchange(reaction, policy):
    """Validate role, identity and sign independently of COBRA's heuristic."""
    stem = _component(reaction.id)
    _require(len(reaction.metabolites) == 1, f'Exchange is not a single-metabolite boundary: {reaction.id}')
    met, coefficient = next(iter(reaction.metabolites.items()))
    expected_metabolite = policy.exchange_metabolites.get(reaction.id, stem + '_e')
    _require(met.id == expected_metabolite and coefficient == -1,
             f'Noncanonical exchange identity/orientation: {reaction.id}')
    _require(met.compartment == policy.external, f'Exchange compartment mismatch: {reaction.id}')
    sbo = _sbo(reaction)
    _require(not any(term in ('SBO:0000628', 'SBO:0000632') for term in sbo),
             f'Exchange conflicts with annotated demand/sink role: {reaction.id}')


def _inventory(model, policy):
    exchanges = []
    for reaction in model.reactions:
        external_boundary = reaction.boundary and any(m.compartment == policy.external for m in reaction.metabolites)
        tagged_exchange = 'SBO:0000627' in _sbo(reaction)
        if reaction.id.startswith('EX_') or external_boundary or tagged_exchange:
            _exchange(reaction, policy)
            exchanges.append(reaction)
    return sorted(exchanges, key=lambda r: r.id)


def _carrier(reaction, extracellular, cytoplasm, capacity):
    _require({m.id: value for m, value in reaction.metabolites.items()} == {extracellular: -1., cytoplasm: 1.}
             and reaction.bounds == (0., capacity) and not reaction.gene_reaction_rule,
             f'Existing completion carrier conflicts with declared recipe: {reaction.id}')


def _chemistry(extracellular, cytoplasm):
    a, b = extracellular.elements, cytoplasm.elements
    if a and b:
        _require(a == b, f'Completion metabolite formula mismatch: {extracellular.id}/{cytoplasm.id}')
    if extracellular.charge is not None and cytoplasm.charge is not None:
        _require(extracellular.charge == cytoplasm.charge,
                 f'Completion metabolite charge mismatch: {extracellular.id}/{cytoplasm.id}')
    return {'extracellular': extracellular.id, 'cytoplasm': cytoplasm.id,
            'formula_metadata_known': bool(a and b),
            'charge_metadata_known': extracellular.charge is not None and cytoplasm.charge is not None,
            'note': 'Known metadata must agree; agreement or placeholders do not establish chemically correct annotation.'}


def prepare_medium(model: cobra.Model, medium: Medium, *, policy: CompartmentPolicy,
                   completion_media: Iterable[Medium] = (), missing_policy: str = 'error') -> PreparedMedium:
    """Return an isolated prepared copy and a complete boundary-bound record.

    Explicit exchange identity bindings permit known source identifier exceptions;
    they do not rename reactions, create nutrients or relax shape/sign validation.
    Existing internal demands and sinks retain their bounds and are reported,
    including any able to supply material. This function does not infer whether
    that internal supply is physiologically justified. Missing components and
    existing exchanges lacking transport are explicit; no transporter is added
    merely because its existing exchange cannot carry flux.

    On any error, the caller's input remains unchanged. External aliases may be
    normalized only when explicitly listed, and only on metabolites ending _e.
    Cytoplasmic aliases are similarly explicit and restricted to _c identifiers.
    No solver is called. Apply strain deletions AFTER preparation: by design,
    medium preparation resets environmental exchange bounds.
    Missing medium components raise unless missing_policy='report' is explicit.
    """
    _require(isinstance(policy.external, str) and policy.external and isinstance(policy.cytoplasm, str)
             and policy.cytoplasm and policy.external != policy.cytoplasm, 'Declare distinct compartment identifiers')
    _require(_finite(policy.secretion_capacity) and policy.secretion_capacity > 0, 'Invalid secretion capacity')
    _require(missing_policy in ('error', 'report'), 'Unknown missing-component policy')
    _require(len(set(policy.external_aliases)) == len(policy.external_aliases)
             and all(isinstance(a, str) and a and a not in (policy.external, policy.cytoplasm) for a in policy.external_aliases),
             'Invalid or ambiguous external compartment aliases')
    _require(len(set(policy.cytoplasm_aliases)) == len(policy.cytoplasm_aliases)
             and all(isinstance(a, str) and a and a not in (policy.external, policy.cytoplasm) for a in policy.cytoplasm_aliases)
             and not set(policy.external_aliases) & set(policy.cytoplasm_aliases),
             'Invalid or ambiguous cytoplasmic compartment aliases')
    media = list(completion_media)
    _validate_medium(medium)
    for item in media:
        _validate_medium(item)
    _require(all(gene.functional for gene in model.genes),
             'Prepare the medium before gene knockouts; inactive gene flags are present')
    _require(isinstance(policy.exchange_metabolites, dict), 'Declare exchange identity bindings as a mapping')
    for rid, mid in policy.exchange_metabolites.items():
        _component(rid)
        _require(rid in model.reactions and isinstance(mid, str) and mid.endswith('_e'),
                 f'Identity binding requires an existing exchange and external metabolite identifier: {rid}')
    copy = model.copy()
    renamed = []
    for met in copy.metabolites:
        if met.compartment in policy.external_aliases:
            _require(met.id.endswith('_e'), f'Alias includes a non-extracellular identifier: {met.id}')
            renamed.append({'metabolite': met.id, 'before': met.compartment, 'after': policy.external})
            met.compartment = policy.external
        if met.compartment in policy.cytoplasm_aliases:
            _require(met.id.endswith('_c'), f'Cytoplasmic alias includes a non-cytoplasmic identifier: {met.id}')
            renamed.append({'metabolite': met.id, 'before': met.compartment, 'after': policy.cytoplasm})
            met.compartment = policy.cytoplasm
        if met.id.endswith('_e'):
            _require(met.compartment == policy.external, f'Undeclared external compartment alias: {met.id}/{met.compartment}')
        if met.id.endswith('_c'):
            _require(met.compartment == policy.cytoplasm, f'Cytoplasmic compartment mismatch: {met.id}')
    _inventory(copy, policy)
    requested = sorted({rid for item in media for rid in item.uptakes})
    added, skipped, chemistry = [], [], []
    for rid in requested:
        stem = _component(rid)
        ext, cyt, carrier_id = stem + '_e', stem + '_c', 'MEDt_' + stem
        if rid in copy.reactions:
            if carrier_id in copy.reactions:
                _carrier(copy.reactions.get_by_id(carrier_id), ext, cyt, policy.secretion_capacity)
                chemistry.append(_chemistry(copy.metabolites.get_by_id(ext), copy.metabolites.get_by_id(cyt)))
            skipped.append({'exchange': rid, 'reason': 'exchange_already_exists'})
            continue
        if stem in policy.completion_exclude or cyt not in copy.metabolites:
            skipped.append({'exchange': rid, 'reason': 'excluded' if stem in policy.completion_exclude else 'no_cytoplasmic_metabolite'})
            continue
        # A name collision must not be silently treated as completed transport.
        _require(carrier_id not in copy.reactions, f'Orphan/colliding completion carrier: {carrier_id}')
        cyt_met = copy.metabolites.get_by_id(cyt)
        if ext not in copy.metabolites:
            copy.add_metabolites([cobra.Metabolite(ext, name=cyt_met.name, compartment=policy.external,
                                                formula=cyt_met.formula, charge=cyt_met.charge)])
        ext_met = copy.metabolites.get_by_id(ext)
        chemistry.append(_chemistry(ext_met, cyt_met))
        exchange = cobra.Reaction(rid, name=f'{stem} exchange (medium preparation v1)',
                                  lower_bound=0., upper_bound=policy.secretion_capacity)
        exchange.add_metabolites({ext_met: -1.})
        carrier = cobra.Reaction(carrier_id, name=f'{stem} uptake (medium preparation v1, no gene)',
                                 lower_bound=0., upper_bound=policy.secretion_capacity)
        carrier.add_metabolites({ext_met: -1., cyt_met: 1.})
        copy.add_reactions([exchange, carrier])
        added.append(rid)
    exchanges = _inventory(copy, policy)
    for reaction in exchanges:
        reaction.bounds = (0., policy.secretion_capacity)
    missing = []
    for rid, lower in medium.uptakes.items():
        if rid not in copy.reactions:
            missing.append(rid)
        else:
            _exchange(copy.reactions.get_by_id(rid), policy)
            copy.reactions.get_by_id(rid).lower_bound = lower
    _require(not missing or missing_policy == 'report', f'Missing medium components: {sorted(missing)}')
    ex_ids = {r.id for r in exchanges}
    boundaries, internal_supply = {}, []
    for reaction in sorted(copy.reactions, key=lambda r: r.id):
        if reaction.boundary:
            boundaries[reaction.id] = {'bounds': list(reaction.bounds),
                'stoichiometry': {m.id: value for m, value in reaction.metabolites.items()},
                'role': 'exchange' if reaction.id in ex_ids else 'internal_boundary'}
            if reaction.id not in ex_ids:
                for coefficient in reaction.metabolites.values():
                    if coefficient * reaction.lower_bound > 0 or coefficient * reaction.upper_bound > 0:
                        internal_supply.append(reaction.id)
    no_connection = []
    for reaction in exchanges:
        met = next(iter(reaction.metabolites))
        if not any(not other.boundary for other in met.reactions):
            no_connection.append(reaction.id)
    report = {'schema_version': 1, 'policy': asdict(policy), 'missing_policy': missing_policy, 'medium': asdict(medium),
              'completion_media': [asdict(item) for item in media], 'added_exchanges': added,
              'completion_skips': skipped, 'normalized_compartments': sorted(renamed, key=lambda r: r['metabolite']),
              'completion_chemistry': chemistry,
              'missing_medium': sorted(missing), 'exchanges_without_network_connection': no_connection,
              'exchange_ids': sorted(ex_ids), 'boundary_reactions': boundaries,
              'internal_supply_capable_boundaries': sorted(set(internal_supply)),
              'scope': 'Canonical BiGG exchange preparation; internal boundaries retained and disclosed. No biological transport validation, optimization or concentration conversion.'}
    return PreparedMedium(copy, report)
