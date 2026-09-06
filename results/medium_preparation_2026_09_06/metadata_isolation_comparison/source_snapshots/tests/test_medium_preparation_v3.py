"""Metadata isolation and unchanged preparation, using synthetic models only."""
from copy import deepcopy
from dataclasses import asdict

import cobra
import pytest

from gembench.media import Medium
from gembench import medium_preparation_v2 as V2
from gembench import medium_preparation_v3 as V3


MEDIUM = Medium('synthetic', '', {'EX_a_e': -2., 'EX_b_e': -.5}, notes=['original'])
POLICY = V3.CompartmentPolicy(external='e', cytoplasm='c',
    external_aliases=('legacy_e',), cytoplasm_aliases=('legacy_c',),
    exchange_metabolites={'EX_alias_e': 'd_e'})


def metadata_objects(model):
    return {'model': model, **{f'{kind}:{item.id}': item
        for kind, items in (('metabolite', model.metabolites), ('reaction', model.reactions),
                            ('gene', model.genes), ('group', model.groups)) for item in items}}


def algebra(model):
    return {'reactions': {r.id: ({m.id: v for m, v in r.metabolites.items()}, r.bounds,
                                r.gene_reaction_rule) for r in model.reactions},
        'metabolites': {m.id: (m.compartment, m.formula, m.charge) for m in model.metabolites},
        'objective': (str(model.objective.expression), model.objective_direction),
        'variables': {v.name: (v.lb, v.ub, v.type) for v in model.variables},
        'constraints': {c.name: (c.lb, c.ub, str(c.expression)) for c in model.constraints},
        'functional': {g.id: g.functional for g in model.genes}}


def snapshot(model):
    return {'algebra': algebra(model), 'compartments': model.compartments,
        'metadata': deepcopy({key: (obj.annotation, obj.notes)
                              for key, obj in metadata_objects(model).items()}),
        'groups': {g.id: sorted((type(member).__name__, member.id) for member in g.members)
                   for g in model.groups}}


@pytest.fixture
def toy():
    model = cobra.Model('metadata_isolation_toy')
    model.solver = 'glpk'
    for mid, compartment in [('a_e', 'legacy_e'), ('a_c', 'legacy_c'),
                             ('b_c', 'legacy_c'), ('d_e', 'legacy_e')]:
        model.add_metabolites([cobra.Metabolite(mid, compartment=compartment, formula='C', charge=0)])
    for rid, coefficients, bounds, rule in [
        ('EX_a_e', {'a_e': -1.}, (-3., 7.), ''),
        ('EX_alias_e', {'d_e': -1.}, (-.1, 4.), ''),
        ('T_a', {'a_e': -1., 'a_c': 1.}, (0., 1000.), 'g_transport'),
        ('BIOMASS_TEST', {'a_c': -1., 'b_c': -1.}, (0., 1000.), '')]:
        reaction = cobra.Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
        reaction.add_metabolites({model.metabolites.get_by_id(mid): value
                                 for mid, value in coefficients.items()})
        reaction.gene_reaction_rule = rule
        model.add_reactions([reaction])
    model.objective = 'BIOMASS_TEST'
    model.compartments = {'legacy_e': 'source extracellular', 'legacy_c': 'source cytoplasm',
                          'e': 'declared extracellular', 'c': 'declared cytoplasm'}
    group = cobra.core.Group('bundle', members=[model.metabolites.a_e, model.reactions.T_a,
                                               model.genes.g_transport])
    outer = cobra.core.Group('outer', members=[group])
    model.add_groups([group, outer])
    for obj in metadata_objects(model).values():
        obj.annotation = {'audit': {'tags': ['original']}}
        obj.notes = {'audit': {'tags': ['original']}}
    return model


def prepared(model, medium=MEDIUM, **kwargs):
    return V3.prepare_medium(model, medium, policy=POLICY, completion_media=[MEDIUM], **kwargs)


@pytest.mark.parametrize('target', ['model', 'metabolite:a_e', 'reaction:T_a', 'gene:g_transport',
                                   'group:bundle', 'group:outer'])
@pytest.mark.parametrize('field', ['annotation', 'notes'])
@pytest.mark.parametrize('direction', ['returned_to_caller', 'caller_to_returned'])
def test_nested_metadata_mutations_are_isolated_both_directions(toy, target, field, direction):
    answer = prepared(toy)
    changed, untouched = (answer.model, toy) if direction == 'returned_to_caller' else (toy, answer.model)
    before = snapshot(untouched)
    metadata = getattr(metadata_objects(changed)[target], field)
    metadata['audit']['tags'].append('changed')
    metadata['new'] = ['new value']
    assert snapshot(untouched) == before


@pytest.mark.parametrize('direction', ['returned_to_caller', 'caller_to_returned'])
def test_compartment_descriptions_are_isolated_both_directions(toy, direction):
    answer = prepared(toy)
    changed, untouched = (answer.model, toy) if direction == 'returned_to_caller' else (toy, answer.model)
    before = snapshot(untouched)
    changed.compartments = {'e': 'changed extracellular', 'legacy_e': 'changed source extracellular'}
    assert snapshot(untouched) == before


def test_nested_groups_and_solver_references_belong_to_returned_model(toy):
    answer = prepared(toy)
    copy = answer.model
    assert copy.groups.outer.members == {copy.groups.bundle}
    assert copy.groups.bundle.members == {copy.metabolites.a_e, copy.reactions.T_a, copy.genes.g_transport}
    for obj in metadata_objects(copy).values():
        if obj is not copy:
            assert obj._model is copy
    for reaction in copy.reactions:
        assert reaction.forward_variable is copy.variables[reaction.id]
        assert reaction.reverse_variable is copy.variables[reaction.reverse_id]
        assert all(m._model is copy for m in reaction.metabolites)
    before = snapshot(toy)
    copy.groups.bundle.remove_members([copy.metabolites.a_e])
    copy.groups.outer.add_members([copy.reactions.BIOMASS_TEST])
    assert snapshot(toy) == before
    before = snapshot(copy)
    toy.groups.bundle.add_members([toy.reactions.BIOMASS_TEST])
    toy.groups.outer.remove_members([toy.groups.bundle])
    assert snapshot(copy) == before


def test_v3_matches_v2_algebra_report_and_repeat_preparation(toy):
    before = snapshot(toy)
    v2 = V2.prepare_medium(toy, MEDIUM, policy=POLICY, completion_media=[MEDIUM])
    v3 = prepared(toy)
    assert snapshot(toy) == before
    assert snapshot(v3.model) == snapshot(v2.model)
    assert v3.report == v2.report
    withdrawn = Medium('withdrawn', '', {'EX_a_e': -2.})
    reused = prepared(v3.model, withdrawn)
    fresh = prepared(toy, withdrawn)
    assert snapshot(reused.model) == snapshot(fresh.model)
    assert reused.model.reactions.EX_b_e.bounds == (0., 1000.)


def test_report_does_not_share_medium_or_policy_data(toy):
    medium, policy = deepcopy(MEDIUM), deepcopy(POLICY)
    before = (asdict(medium), asdict(policy))
    answer = V3.prepare_medium(toy, medium, policy=policy, completion_media=[medium])
    answer.report['medium']['notes'].append('report change')
    answer.report['completion_media'][0]['uptakes']['EX_a_e'] = -7.
    answer.report['policy']['exchange_metabolites']['EX_alias_e'] = 'changed_e'
    assert (asdict(medium), asdict(policy)) == before
    report_before = deepcopy(answer.report)
    medium.notes.append('caller change')
    policy.exchange_metabolites['EX_alias_e'] = 'other_e'
    assert answer.report == report_before


@pytest.mark.parametrize('failure', [False, True])
def test_caller_context_and_bounds_survive_success_and_failure(toy, failure):
    original = snapshot(toy)
    with toy:
        toy.reactions.T_a.upper_bound = .25
        inside = snapshot(toy)
        context = toy._contexts[-1]
        pending = context.size()
        if failure:
            absent = Medium('missing', '', {**MEDIUM.uptakes, 'EX_absent_e': -1.})
            with pytest.raises(ValueError, match='Missing medium components'):
                prepared(toy, absent)
        else:
            answer = prepared(toy)
            assert answer.model._contexts == []
            with answer.model:
                answer.model.reactions.T_a.upper_bound = .1
            assert answer.model.reactions.T_a.upper_bound == .25
            answer.model.reactions.T_a.upper_bound = .2
        assert snapshot(toy) == inside
        assert toy._contexts[-1] is context
        assert context.size() == pending
    assert snapshot(toy) == original


def test_deepcopy_failure_does_not_expose_partial_result_or_mutate_caller(toy, monkeypatch):
    before = snapshot(toy)
    def fail_copy(value):
        raise RuntimeError('synthetic copy failure')
    monkeypatch.setattr(V3, 'deepcopy', fail_copy)
    with pytest.raises(RuntimeError, match='synthetic copy failure'):
        prepared(toy)
    assert snapshot(toy) == before


def test_preparation_does_not_optimize(toy, monkeypatch):
    def fail_solve(*args, **kwargs):
        raise AssertionError('preparation must not optimize')
    monkeypatch.setattr(cobra.Model, 'optimize', fail_solve)
    monkeypatch.setattr(cobra.Model, 'slim_optimize', fail_solve)
    monkeypatch.setattr(type(toy.solver), 'optimize', fail_solve)
    answer = prepared(toy)
    assert answer.model.reactions.EX_a_e.bounds == (-2., 1000.)


def test_deepcopied_solver_solves_ordinary_synthetic_lp_independently(toy):
    before = snapshot(toy)
    answer = prepared(toy)
    copy = answer.model
    assert copy.solver is not toy.solver
    assert copy.solver.configuration is not toy.solver.configuration
    assert copy.slim_optimize(error_value=None) == pytest.approx(.5)
    with copy:
        copy.reactions.T_a.upper_bound = .2
        assert copy.slim_optimize(error_value=None) == pytest.approx(.2)
    assert copy.slim_optimize(error_value=None) == pytest.approx(.5)
    assert snapshot(toy) == before
    toy.reactions.T_a.upper_bound = .1
    assert copy.reactions.T_a.upper_bound == 1000.
    assert copy.slim_optimize(error_value=None) == pytest.approx(.5)


def test_nondefault_objective_and_custom_constraints_match_v2(toy):
    toy.objective = {toy.reactions.BIOMASS_TEST: 2., toy.reactions.T_a: 3.}
    toy.objective_direction = 'min'
    constraint = toy.problem.Constraint(toy.reactions.T_a.flux_expression +
        toy.reactions.BIOMASS_TEST.flux_expression, lb=.1, ub=.7, name='declared_coupling')
    toy.add_cons_vars([constraint])
    before = snapshot(toy)
    v2 = V2.prepare_medium(toy, MEDIUM, policy=POLICY, completion_media=[MEDIUM])
    v3 = prepared(toy)
    assert algebra(v3.model) == algebra(v2.model)
    assert v3.report == v2.report
    v3.model.constraints.declared_coupling.ub = .6
    assert snapshot(toy) == before
    assert v2.model.constraints.declared_coupling.ub == .7
