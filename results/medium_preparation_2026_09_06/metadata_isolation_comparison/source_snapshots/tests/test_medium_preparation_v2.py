"""An explicit source exchange/metabolite exception must not relax other guards."""
from dataclasses import replace
import json

import cobra
import pytest

from gembench.media import Medium
from gembench import medium_preparation as V1
from gembench import medium_preparation_v2 as V2


MEDIUM = Medium('synthetic', '', {'EX_glc__D_e': -2.})


def add(model, rid, stoichiometry, bounds=(0., 1000.), rule=''):
    reaction = cobra.Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
    reaction.add_metabolites({model.metabolites.get_by_id(mid): value for mid, value in stoichiometry.items()})
    reaction.gene_reaction_rule = rule
    model.add_reactions([reaction])
    return reaction


@pytest.fixture
def toy():
    model = cobra.Model('existing_exchange_exception_toy')
    for mid, compartment, formula in [('glc__D_e', 'e', 'C6H12O6'), ('glc__D_c', 'c', 'C6H12O6'),
                                      ('2ameph_e', 'e', 'C2H8NO3P'), ('2ameph_c', 'c', 'C2H8NO3P'),
                                      ('btn_c', 'c', 'C10H16N2O3S')]:
        model.add_metabolites([cobra.Metabolite(mid, compartment=compartment, formula=formula, charge=0)])
    add(model, 'EX_glc__D_e', {'glc__D_e': -1.}, (-3., 7.))
    add(model, 'EX_AEP_e', {'2ameph_e': -1.}, (-.2, 4.), 'g_source_exchange')
    add(model, 'TGLC', {'glc__D_e': -1., 'glc__D_c': 1.}, rule='g_transport')
    add(model, 'TAEP', {'2ameph_e': -1., '2ameph_c': 1.})
    add(model, 'BIOMASS_TEST', {'glc__D_c': -1.})
    model.objective = 'BIOMASS_TEST'
    return model


def policy(**kwargs):
    return V2.CompartmentPolicy(external='e', cytoplasm='c', exchange_metabolites={'EX_AEP_e': '2ameph_e'}, **kwargs)


def snapshot(model):
    return json.dumps({'model': cobra.io.dict.model_to_dict(model),
                       'functional': {g.id: g.functional for g in model.genes},
                       'objective': str(model.objective.expression), 'sense': model.objective_direction}, sort_keys=True)


def signature(model):
    return ({r.id: ({m.id: v for m, v in r.metabolites.items()}, r.bounds, r.gene_reaction_rule)
             for r in model.reactions},
            {m.id: (m.compartment, m.formula, m.charge) for m in model.metabolites},
            {r.id: v for r, v in cobra.util.solver.linear_reaction_coefficients(model).items()}, model.objective_direction)


def test_noncanonical_binding_requires_explicit_exception(toy):
    before = snapshot(toy)
    with pytest.raises(ValueError, match='identity|binding|metabolite'):
        V2.prepare_medium(toy, MEDIUM, policy=V2.CompartmentPolicy(external='e', cytoplasm='c'))
    assert snapshot(toy) == before


def test_exact_exception_preserves_source_ids_and_opens_no_unrequested_supply(toy):
    before = snapshot(toy)
    answer = V2.prepare_medium(toy, MEDIUM, policy=policy())
    assert snapshot(toy) == before
    assert {r.id for r in answer.model.reactions} == {r.id for r in toy.reactions}
    assert {m.id for m in answer.model.metabolites} == {m.id for m in toy.metabolites}
    retained = answer.model.reactions.EX_AEP_e
    assert {m.id: v for m, v in retained.metabolites.items()} == {'2ameph_e': -1.}
    assert retained.gene_reaction_rule == toy.reactions.EX_AEP_e.gene_reaction_rule
    assert retained.bounds == (0., 1000.)
    assert answer.report['added_exchanges'] == []
    assert answer.report['policy']['exchange_metabolites'] == {'EX_AEP_e': '2ameph_e'}
    assert answer.report['boundary_reactions']['EX_AEP_e']['stoichiometry'] == {'2ameph_e': -1.}


@pytest.mark.parametrize('binding', [
    {'EX_AEP_e': 'AEP_e'}, {'EX_AEP_e': 'glc__D_e'}, {'EX_AEP_e': 'not_present_e'},
    {'EX_not_present_e': '2ameph_e'}, {'EX_AEP_e': '2ameph_c'}, {'TAEP': '2ameph_e'}])
def test_wrong_or_missing_bindings_fail_without_mutation(toy, binding):
    before = snapshot(toy)
    with pytest.raises(ValueError):
        V2.prepare_medium(toy, MEDIUM, policy=replace(policy(), exchange_metabolites=binding))
    assert snapshot(toy) == before


@pytest.mark.parametrize('shape', ['positive', 'scaled', 'multiple'])
def test_binding_exception_never_waives_boundary_shape_or_orientation(toy, shape):
    reaction = toy.reactions.EX_AEP_e
    if shape == 'positive':
        reaction.add_metabolites({toy.metabolites.get_by_id('2ameph_e'): 2.})
    elif shape == 'scaled':
        reaction.add_metabolites({toy.metabolites.get_by_id('2ameph_e'): -1.})
    else:
        reaction.add_metabolites({toy.metabolites.get_by_id('2ameph_c'): 1.})
    before = snapshot(toy)
    with pytest.raises(ValueError):
        V2.prepare_medium(toy, MEDIUM, policy=policy())
    assert snapshot(toy) == before


def test_exception_cannot_move_an_existing_internal_metabolite_outside(toy):
    toy.metabolites.get_by_id('2ameph_e').compartment = 'c'
    before = snapshot(toy)
    with pytest.raises(ValueError):
        V2.prepare_medium(toy, MEDIUM, policy=policy())
    assert snapshot(toy) == before


@pytest.mark.parametrize('sbo', ['SBO:0000628', 'sbo:0000632'])
def test_exception_does_not_override_demand_or_sink_annotation(toy, sbo):
    toy.reactions.EX_AEP_e.annotation['sbo'] = sbo
    before = snapshot(toy)
    with pytest.raises(ValueError, match='role|demand|sink'):
        V2.prepare_medium(toy, MEDIUM, policy=policy())
    assert snapshot(toy) == before


def test_override_cannot_be_reserved_for_a_future_completion_reaction(toy):
    before = snapshot(toy)
    extra = {'EX_AEP_e': '2ameph_e', 'EX_btn_e': 'btn_e'}
    vitamin = Medium('vitamin', '', {**MEDIUM.uptakes, 'EX_btn_e': -.001})
    with pytest.raises(ValueError):
        V2.prepare_medium(toy, vitamin, policy=replace(policy(), exchange_metabolites=extra), completion_media=[vitamin])
    assert snapshot(toy) == before
    assert 'EX_btn_e' not in toy.reactions


def test_exception_cannot_rebind_source_stoichiometry_even_if_metabolite_exists(toy):
    toy.add_metabolites([cobra.Metabolite('AEP_e', compartment='e', formula='C2H8NO3P', charge=0)])
    before = snapshot(toy)
    with pytest.raises(ValueError):
        V2.prepare_medium(toy, MEDIUM, policy=replace(policy(), exchange_metabolites={'EX_AEP_e': 'AEP_e'}))
    assert snapshot(toy) == before


@pytest.mark.parametrize('malformed', [[], None, {'EX_AEP_e': False}])
def test_identity_bindings_require_a_mapping_of_identifiers(toy, malformed):
    before = snapshot(toy)
    with pytest.raises(ValueError):
        V2.prepare_medium(toy, MEDIUM, policy=replace(policy(), exchange_metabolites=malformed))
    assert snapshot(toy) == before


def test_one_valid_exception_does_not_relax_another_exchange_identity(toy):
    add(toy, 'EX_unrelated_e', {'glc__D_e': -1.})
    before = snapshot(toy)
    with pytest.raises(ValueError, match='identity|orientation'):
        V2.prepare_medium(toy, MEDIUM, policy=policy())
    assert snapshot(toy) == before


def test_empty_overrides_reproduce_v1_formulation_with_completion(toy):
    toy.reactions.EX_AEP_e.id = 'EX_2ameph_e'
    before = snapshot(toy)
    vitamin = Medium('vitamin', '', {**MEDIUM.uptakes, 'EX_btn_e': -.001})
    first = V1.prepare_medium(toy, vitamin, policy=V1.CompartmentPolicy(external='e', cytoplasm='c'), completion_media=[vitamin])
    second = V2.prepare_medium(toy, vitamin, policy=V2.CompartmentPolicy(external='e', cytoplasm='c'), completion_media=[vitamin])
    assert snapshot(toy) == before
    assert signature(first.model) == signature(second.model)
    assert first.report['boundary_reactions'] == second.report['boundary_reactions']
    assert first.report['added_exchanges'] == second.report['added_exchanges']


def test_explicit_exchange_bounds_reset_on_reuse_without_renaming(toy):
    supplied = Medium('supplied', '', {**MEDIUM.uptakes, 'EX_AEP_e': -.1})
    started = V2.prepare_medium(toy, supplied, policy=policy())
    assert started.model.reactions.EX_AEP_e.lower_bound == -.1
    reused = V2.prepare_medium(started.model, MEDIUM, policy=policy())
    fresh = V2.prepare_medium(toy, MEDIUM, policy=policy())
    assert signature(reused.model) == signature(fresh.model)
    assert reused.model.reactions.EX_AEP_e.bounds == (0., 1000.)
