"""Preparation safety contracts, exercised only with small synthetic models."""
from dataclasses import replace
import json

import cobra
import pytest

from gembench.media import Medium
from gembench.medium_preparation import CompartmentPolicy, prepare_medium


POLICY = CompartmentPolicy(external='C_e', cytoplasm='C_c')
CARBON = Medium('carbon', 'Synthetic glucose medium', {'EX_glc__D_e': -2.})
SUPPLEMENTED = Medium('supplemented', 'Synthetic glucose and vitamin medium',
                      {'EX_glc__D_e': -2., 'EX_btn_e': -.5})


def reaction(model, rid, stoichiometry, bounds=(0., 1000.), rule=''):
    value = cobra.Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
    value.add_metabolites({model.metabolites.get_by_id(mid): coefficient
                           for mid, coefficient in stoichiometry.items()})
    value.gene_reaction_rule = rule
    model.add_reactions([value])
    return value


@pytest.fixture
def toy():
    model = cobra.Model('explicit_medium_toy')
    for mid, compartment, formula, charge in [
        ('glc__D_e', 'C_e', 'C6H12O6', 0), ('glc__D_c', 'C_c', 'C6H12O6', 0),
        ('co2_e', 'C_e', 'CO2', 0), ('btn_c', 'C_c', 'C10H16N2O3S', 0)]:
        model.add_metabolites([cobra.Metabolite(mid, compartment=compartment, formula=formula, charge=charge)])
    reaction(model, 'EX_glc__D_e', {'glc__D_e': -1.}, (-3., 7.))
    reaction(model, 'EX_co2_e', {'co2_e': -1.}, (0., 4.))
    reaction(model, 'TGLC', {'glc__D_e': -1., 'glc__D_c': 1.}, rule='g_transport')
    reaction(model, 'Growth', {'glc__D_c': -1., 'btn_c': -1.})
    model.objective = 'Growth'
    model.annotation = {'audit_fixture': ['synthetic']}
    return model


def snapshot(model):
    """Include genotype and objective state omitted by ordinary model JSON."""
    return json.dumps({'model': cobra.io.dict.model_to_dict(model),
        'gene_functional': {g.id: g.functional for g in model.genes},
        'objective_direction': model.objective_direction,
        'objective': str(model.objective.expression),
        'constraints': {c.name: [c.lb, c.ub, str(c.expression)] for c in model.constraints}},
        sort_keys=True, allow_nan=False)


def lp_signature(model):
    return ({r.id: ({m.id: float(v) for m, v in r.metabolites.items()}, r.bounds, r.gene_reaction_rule)
             for r in model.reactions},
            {m.id: (m.compartment, m.formula, m.charge) for m in model.metabolites},
            {r.id: float(v) for r, v in cobra.util.solver.linear_reaction_coefficients(model).items()},
            model.objective_direction)


def completed(model, medium=SUPPLEMENTED, **kwargs):
    return prepare_medium(model, medium, policy=POLICY, completion_media=[SUPPLEMENTED], **kwargs)


def test_completion_uses_explicit_compartment_and_isolated_copy(toy):
    before = snapshot(toy)
    answer = completed(toy)
    assert snapshot(toy) == before
    assert answer.model is not toy
    assert answer.model.metabolites.btn_e.compartment == 'C_e'
    assert answer.model.metabolites.btn_e.formula == toy.metabolites.btn_c.formula
    assert answer.model.metabolites.btn_e.charge == toy.metabolites.btn_c.charge
    assert answer.model.reactions.EX_btn_e in answer.model.exchanges
    assert answer.model.reactions.MEDt_btn.bounds == (0., 1000.)
    assert answer.model.reactions.MEDt_btn.check_mass_balance() == {}
    assert not answer.model.reactions.MEDt_btn.gene_reaction_rule
    assert answer.model.slim_optimize() == pytest.approx(.5)
    assert snapshot(toy) == before
    answer.model.metabolites.btn_c.name = 'changed in copy'
    answer.model.genes.g_transport.knock_out()
    assert snapshot(toy) == before


def test_medium_switching_is_independent_of_prior_supplement_and_repeat_preparation(toy):
    supplied = completed(toy)
    supplied_before = snapshot(supplied.model)
    switched = completed(supplied.model, CARBON)
    fresh = completed(toy, CARBON)
    again = completed(switched.model, CARBON)
    assert snapshot(supplied.model) == supplied_before
    assert lp_signature(switched.model) == lp_signature(fresh.model) == lp_signature(again.model)
    assert switched.model.reactions.EX_btn_e.bounds == (0., 1000.)
    assert switched.model.slim_optimize() == pytest.approx(0.)
    assert completed(switched.model).model.slim_optimize() == pytest.approx(.5)
    assert switched.report['added_exchanges'] == []


def test_switch_without_repeating_completion_still_closes_completed_boundary(toy):
    supplied = completed(toy).model
    switched = prepare_medium(supplied, CARBON, policy=POLICY)
    assert switched.model.reactions.EX_btn_e.bounds == (0., 1000.)
    assert switched.model.slim_optimize() == pytest.approx(0.)


def test_explicit_legacy_alias_normalization_closes_hidden_exchange(toy):
    legacy = completed(toy).model
    legacy.metabolites.btn_e.compartment = 'e'
    before = snapshot(legacy)
    with pytest.raises(ValueError, match='alias'):
        completed(legacy, CARBON)
    assert snapshot(legacy) == before
    answer = prepare_medium(legacy, CARBON, policy=replace(POLICY, external_aliases=('e',)),
                            completion_media=[SUPPLEMENTED])
    assert snapshot(legacy) == before
    assert answer.report['normalized_compartments'] == [{'metabolite': 'btn_e', 'before': 'e', 'after': 'C_e'}]
    assert answer.model.reactions.EX_btn_e.lower_bound == 0
    assert answer.model.reactions.EX_btn_e in answer.model.exchanges
    assert answer.model.slim_optimize() == pytest.approx(0.)


def test_alias_cannot_relabel_an_internal_metabolite(toy):
    toy.metabolites.btn_c.compartment = 'e'
    before = snapshot(toy)
    with pytest.raises(ValueError, match='non-extracellular'):
        prepare_medium(toy, CARBON, policy=replace(POLICY, external_aliases=('e',)))
    assert snapshot(toy) == before


def test_cytoplasmic_alias_normalizes_only_when_declared_and_preserves_chemistry(toy):
    toy.metabolites.btn_c.compartment = 'c'
    before = snapshot(toy)
    with pytest.raises(ValueError, match='Cytoplasmic compartment mismatch'):
        completed(toy)
    answer = prepare_medium(toy, SUPPLEMENTED, policy=replace(POLICY, cytoplasm_aliases=('c',)),
                            completion_media=[SUPPLEMENTED])
    assert snapshot(toy) == before
    assert answer.model.metabolites.btn_c.compartment == 'C_c'
    assert answer.model.metabolites.btn_c.formula == toy.metabolites.btn_c.formula
    assert answer.model.metabolites.btn_c.charge == toy.metabolites.btn_c.charge
    assert answer.report['normalized_compartments'] == [{'metabolite': 'btn_c', 'before': 'c', 'after': 'C_c'}]
    assert answer.model.slim_optimize() == pytest.approx(.5)


@pytest.mark.parametrize('invalid', ['overlapping_roles', 'external_identifier'])
def test_cytoplasmic_alias_cannot_overlap_roles_or_relabel_external_identifier(toy, invalid):
    policy = replace(POLICY, cytoplasm_aliases=('c',))
    if invalid == 'overlapping_roles':
        policy = replace(policy, external_aliases=('c',))
    else:
        toy.metabolites.glc__D_e.compartment = 'c'
    before = snapshot(toy)
    with pytest.raises(ValueError, match='alias'):
        prepare_medium(toy, CARBON, policy=policy)
    assert snapshot(toy) == before


@pytest.mark.parametrize('change', [
    {'external': ''}, {'cytoplasm': 'C_e'}, {'external_aliases': ('e', 'e')},
    {'external_aliases': ('C_c',)}, {'external_aliases': ('C_e',)},
    {'secretion_capacity': 0.}, {'secretion_capacity': float('inf')}, {'secretion_capacity': True}])
def test_invalid_policy_does_not_touch_input(toy, change):
    before = snapshot(toy)
    with pytest.raises(ValueError):
        prepare_medium(toy, CARBON, policy=replace(POLICY, **change))
    assert snapshot(toy) == before


@pytest.mark.parametrize('bound', [1., float('nan'), float('inf'), float('-inf'), True, '-1'])
def test_invalid_uptake_never_partially_mutates_input(toy, bound):
    before = snapshot(toy)
    medium = Medium('bad', '', {'EX_glc__D_e': -5., 'EX_co2_e': bound})
    with pytest.raises(ValueError, match='finite nonpositive'):
        prepare_medium(toy, medium, policy=POLICY)
    assert snapshot(toy) == before


@pytest.mark.parametrize('identifier', ['TGLC', 'EX_glc__D_c', 'EX__e', 'EX_glc__D', 'glc__D_e'])
def test_invalid_medium_identifier_cannot_change_internal_reactions(toy, identifier):
    before = snapshot(toy)
    with pytest.raises(ValueError, match='canonical'):
        prepare_medium(toy, Medium('bad', '', {identifier: -1.}), policy=POLICY)
    assert snapshot(toy) == before


@pytest.mark.parametrize('kind', ['reversed', 'scaled', 'multiple', 'wrong_identity', 'sink_role', 'demand_role', 'lowercase_demand_role'])
def test_malformed_exchange_shape_or_role_rejected(toy, kind):
    ex = toy.reactions.EX_glc__D_e
    if kind == 'reversed':
        ex.add_metabolites({toy.metabolites.glc__D_e: 2.})
    elif kind == 'scaled':
        ex.add_metabolites({toy.metabolites.glc__D_e: -1.})
    elif kind == 'multiple':
        ex.add_metabolites({toy.metabolites.glc__D_c: 1.})
    elif kind == 'wrong_identity':
        ex.add_metabolites({toy.metabolites.glc__D_e: 1., toy.metabolites.co2_e: -1.})
    elif kind == 'lowercase_demand_role':
        ex.annotation['sbo'] = 'sbo:0000628'
    else:
        ex.annotation['sbo'] = 'SBO:0000632' if kind == 'sink_role' else ['SBO:0000628']
    before = snapshot(toy)
    with pytest.raises(ValueError, match='Exchange|exchange'):
        prepare_medium(toy, CARBON, policy=POLICY)
    assert snapshot(toy) == before


@pytest.mark.parametrize('annotation', ['SBO:0000627', 'sbo:0000627'])
def test_exchange_annotation_on_internal_reaction_is_rejected(toy, annotation):
    toy.reactions.TGLC.annotation['sbo'] = annotation
    assert toy.reactions.TGLC in toy.exchanges
    before = snapshot(toy)
    with pytest.raises(ValueError, match='canonical'):
        completed(toy)
    assert snapshot(toy) == before


def test_unrecognized_external_boundary_cannot_remain_hidden(toy):
    toy.reactions.EX_co2_e.id = 'arbitrary_source'
    before = snapshot(toy)
    with pytest.raises(ValueError, match='canonical'):
        prepare_medium(toy, CARBON, policy=POLICY)
    assert snapshot(toy) == before


def test_inactive_gene_flags_require_preparation_before_deletion(toy):
    toy.genes.g_transport.knock_out()
    before = snapshot(toy)
    with pytest.raises(ValueError, match='before gene knockouts'):
        completed(toy)
    assert snapshot(toy) == before
    assert toy.reactions.TGLC.bounds == (0., 0.)


def test_orphan_carrier_collision_rejected_without_added_exchange(toy):
    reaction(toy, 'MEDt_btn', {'btn_c': -1., 'glc__D_c': 1.})
    before = snapshot(toy)
    with pytest.raises(ValueError, match='colliding'):
        completed(toy)
    assert snapshot(toy) == before
    assert 'EX_btn_e' not in toy.reactions


@pytest.mark.parametrize('change', ['reverse_bound', 'wrong_coefficient', 'gene_rule'])
def test_existing_managed_carrier_must_match_recipe(toy, change):
    prepared = completed(toy).model
    carrier = prepared.reactions.MEDt_btn
    if change == 'reverse_bound':
        carrier.lower_bound = -1000.
    elif change == 'wrong_coefficient':
        carrier.add_metabolites({prepared.metabolites.btn_c: 1.})
    else:
        carrier.gene_reaction_rule = 'new_carrier_gene'
    before = snapshot(prepared)
    with pytest.raises(ValueError, match='carrier conflicts'):
        completed(prepared)
    assert snapshot(prepared) == before


@pytest.mark.parametrize('field,bad_value', [('formula', 'C2H4'), ('charge', -1)])
def test_existing_extracellular_metadata_cannot_create_unbalanced_completion(toy, field, bad_value):
    cyt = toy.metabolites.btn_c
    external = cobra.Metabolite('btn_e', compartment='C_e', formula=cyt.formula, charge=cyt.charge)
    setattr(external, field, bad_value)
    toy.add_metabolites([external])
    before = snapshot(toy)
    with pytest.raises(ValueError, match='metadata|formula|charge|chemistry|balance'):
        completed(toy)
    assert snapshot(toy) == before


def test_existing_exchange_without_connection_is_reported_without_invented_carrier(toy):
    cyt = toy.metabolites.btn_c
    toy.add_metabolites([cobra.Metabolite('btn_e', compartment='C_e', formula=cyt.formula, charge=cyt.charge)])
    reaction(toy, 'EX_btn_e', {'btn_e': -1.})
    answer = completed(toy)
    assert 'EX_btn_e' in answer.report['exchanges_without_network_connection']
    assert 'MEDt_btn' not in answer.model.reactions
    assert answer.model.slim_optimize() == pytest.approx(0.)


def test_chemistry_compares_composition_and_reports_unknown_metadata(toy):
    toy.add_metabolites([cobra.Metabolite('btn_e', compartment='C_e', formula='H16C10O3N2S', charge=0)])
    answer = completed(toy)
    assert answer.model.reactions.MEDt_btn.check_mass_balance() == {}
    assert answer.report['completion_chemistry'][0]['formula_metadata_known']
    toy.metabolites.btn_e.formula = None
    toy.metabolites.btn_e.charge = None
    answer = completed(toy)
    assert not answer.report['completion_chemistry'][0]['formula_metadata_known']
    assert not answer.report['completion_chemistry'][0]['charge_metadata_known']


def test_late_missing_component_error_does_not_leak_completed_copy_to_caller(toy):
    before = snapshot(toy)
    bad = Medium('bad', '', {**SUPPLEMENTED.uptakes, 'EX_absent_e': -1.})
    with pytest.raises(ValueError, match='[Mm]issing'):
        prepare_medium(toy, bad, policy=POLICY, completion_media=[SUPPLEMENTED])
    assert snapshot(toy) == before
    assert 'btn_e' not in toy.metabolites
    assert 'EX_btn_e' not in toy.reactions


def test_all_completion_media_are_validated_before_preparation(toy):
    before = snapshot(toy)
    bad = Medium('bad', '', {'EX_btn_e': float('nan')})
    with pytest.raises(ValueError, match='finite nonpositive'):
        prepare_medium(toy, CARBON, policy=POLICY, completion_media=[SUPPLEMENTED, bad])
    assert snapshot(toy) == before


def test_missing_components_error_by_default_report_only_if_requested(toy):
    before = snapshot(toy)
    with pytest.raises(ValueError, match='[Mm]issing'):
        prepare_medium(toy, SUPPLEMENTED, policy=POLICY)
    assert snapshot(toy) == before
    answer = prepare_medium(toy, SUPPLEMENTED, policy=POLICY, missing_policy='report')
    assert answer.report['missing_medium'] == ['EX_btn_e']
    assert answer.model.reactions.EX_glc__D_e.lower_bound == -2.
    assert snapshot(toy) == before
    with pytest.raises(ValueError, match='policy'):
        prepare_medium(toy, CARBON, policy=POLICY, missing_policy='ignore')


def test_completion_exclusions_and_absent_cytoplasm_are_explicit(toy):
    for stem in ('pnto__R', 'fol', 'hco3'):
        toy.add_metabolites([cobra.Metabolite(stem + '_c', compartment='C_c')])
    missing = Medium('unavailable', '', {f'EX_{stem}_e': -1. for stem in ('pnto__R', 'fol', 'hco3', 'unknown')})
    answer = prepare_medium(toy, missing, policy=POLICY, completion_media=[missing], missing_policy='report')
    assert answer.report['added_exchanges'] == []
    assert {row['exchange']: row['reason'] for row in answer.report['completion_skips']} == {
        'EX_pnto__R_e': 'excluded', 'EX_fol_e': 'excluded', 'EX_hco3_e': 'excluded',
        'EX_unknown_e': 'no_cytoplasmic_metabolite'}


def test_internal_boundaries_retained_and_supply_directions_disclosed(toy):
    reaction(toy, 'sink_btn_c', {'btn_c': -1.}, (-.25, .1))
    reaction(toy, 'source_glc', {'glc__D_c': 1.}, (0., .4))
    reaction(toy, 'DM_glc', {'glc__D_c': -1.}, (0., .2))
    answer = completed(toy)
    assert answer.report['internal_supply_capable_boundaries'] == ['sink_btn_c', 'source_glc']
    assert answer.model.reactions.sink_btn_c.bounds == (-.25, .1)
    assert answer.model.reactions.source_glc.bounds == (0., .4)
    assert answer.model.reactions.DM_glc.bounds == (0., .2)
    assert set(answer.report['boundary_reactions']) == {r.id for r in answer.model.reactions if r.boundary}
    for rid in ('sink_btn_c', 'source_glc', 'DM_glc'):
        assert answer.report['boundary_reactions'][rid]['role'] == 'internal_boundary'


def test_declared_export_capacity_is_explicit_and_applies_to_completed_carriers(toy):
    answer = prepare_medium(toy, SUPPLEMENTED, policy=replace(POLICY, secretion_capacity=5.),
                            completion_media=[SUPPLEMENTED])
    assert answer.model.reactions.EX_glc__D_e.bounds == (-2., 5.)
    assert answer.model.reactions.EX_co2_e.bounds == (0., 5.)
    assert answer.model.reactions.MEDt_btn.bounds == (0., 5.)
    assert answer.report['policy']['secretion_capacity'] == 5.


def test_completed_model_survives_sbml_roundtrip_and_subsequent_medium_switch(toy, tmp_path):
    prepared = completed(toy).model
    path = tmp_path / 'completed.xml.gz'
    cobra.io.write_sbml_model(prepared, str(path))
    restored = cobra.io.read_sbml_model(str(path))
    assert lp_signature(restored) == lp_signature(prepared)
    assert restored.reactions.EX_btn_e in restored.exchanges
    switched = completed(restored, CARBON)
    assert switched.model.reactions.EX_btn_e.lower_bound == 0.
    assert switched.model.slim_optimize() == pytest.approx(0.)


def test_nondefault_objective_and_additional_constraints_are_preserved(toy):
    toy.objective = {toy.reactions.Growth: 2., toy.reactions.TGLC: 3.}
    toy.objective_direction = 'min'
    constraint = toy.problem.Constraint(toy.reactions.TGLC.flux_expression + toy.reactions.Growth.flux_expression,
                                        lb=.1, ub=.7, name='declared_coupling')
    toy.add_cons_vars([constraint])
    before = snapshot(toy)
    answer = completed(toy)
    assert snapshot(toy) == before
    assert answer.model.objective_direction == 'min'
    assert {r.id: float(v) for r, v in cobra.util.solver.linear_reaction_coefficients(answer.model).items()} == {
        'Growth': 2., 'TGLC': 3.}
    copied = answer.model.constraints.declared_coupling
    assert (copied.lb, copied.ub, str(copied.expression)) == (constraint.lb, constraint.ub, str(constraint.expression))
