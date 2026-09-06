from copy import deepcopy
from types import SimpleNamespace

import cobra
import pytest

from gembench.gene_mapping import GeneMap
from scripts import run_ppnp_repair as R


@pytest.fixture
def chemistry():
    specs = {'adn_c': ('C10H13N5O4', 0), 'ade_c': ('C5H5N5', 0),
             'ins_c': ('C10H12N4O5', 0), 'hxan_c': ('C5H4N4O', 0),
             'pi_c': ('HO4P', -2), 'r1p_c': ('C5H9O8P', -2)}
    source = cobra.Model('source')
    mets = {mid: cobra.Metabolite(mid, formula=f, charge=c, compartment='C_c') for mid, (f, c) in specs.items()}
    for rid, nucleoside, base in [('PUNP1', 'adn_c', 'ade_c'), ('PUNP5', 'ins_c', 'hxan_c')]:
        r = cobra.Reaction(rid, lower_bound=-1000, upper_bound=1000)
        r.add_metabolites({mets[nucleoside]: -1, mets['pi_c']: -1, mets[base]: 1, mets['r1p_c']: 1})
        source.add_reactions([r])
    target = cobra.Model('parent')
    target.add_metabolites([m.copy() for m in source.metabolites])
    for m in target.metabolites:
        m.charge = 0
    for mid in specs:
        r = cobra.Reaction(f'EX_{mid}', lower_bound=0, upper_bound=1000)
        r.add_metabolites({target.metabolites.get_by_id(mid): -1})
        target.add_reactions([r])
    growth = cobra.Reaction('Growth', lower_bound=0, upper_bound=1000)
    growth.add_metabolites({target.metabolites.r1p_c: -1})
    target.add_reactions([growth])
    target.objective = 'Growth'
    candidate = {'reaction_records': [R.reaction_record(r) for r in source.reactions],
        'all_involved_metabolites': {m.id: R.met_record(m) for m in source.metabolites},
        'target_metabolite_map': {m.id: m.id for m in source.metabolites},
        'expected_target_metabolites': {m.id: R.met_record(m) for m in target.metabolites},
        'proposed_reactions': [dict(id=r.id, source_id=r.id, substrate=s,
            metabolites={m.id: c for m, c in r.metabolites.items()}) for r, s in zip(source.reactions, ['adenosine', 'inosine'])]}
    return target, source, candidate, GeneMap('Putida', {}, [])


def applied(chemistry, direction='forward'):
    base, source, candidate, mapping = chemistry
    config = dict(id='toy', substrates=['adenosine', 'inosine'], direction=direction, close_rhcys=False)
    return R.apply_arm(base, source, candidate, mapping, config, {'PP_4248'})[0]


def test_phosphorolysis_requires_phosphate_and_candidate_gene(chemistry):
    model = applied(chemistry)
    model.reactions.EX_adn_c.lower_bound = -10
    assert model.slim_optimize() == pytest.approx(0)
    model.reactions.EX_pi_c.lower_bound = -10
    assert model.slim_optimize() == pytest.approx(10)
    with model:
        model.genes.PP_4248.knock_out()
        assert model.slim_optimize() == pytest.approx(0)
    assert model.slim_optimize() == pytest.approx(10)


def test_forward_primary_does_not_silently_allow_reverse_synthesis(chemistry):
    for direction, expected in [('forward', 0), ('reversible', 10)]:
        model = applied(chemistry, direction)
        model.reactions.EX_r1p_c.lower_bound = -10
        model.reactions.EX_ade_c.lower_bound = -10
        model.objective = 'EX_adn_c'
        assert model.slim_optimize() == pytest.approx(expected)


def test_no_free_ribose_phosphorylation_is_invented(chemistry):
    model = applied(chemistry)
    ribose = cobra.Metabolite('rib__D_c', formula='C5H10O5', charge=0, compartment='C_c')
    feed = cobra.Reaction('EX_rib__D_c', lower_bound=-10, upper_bound=1000)
    feed.add_metabolites({ribose: -1})
    model.add_reactions([feed])
    model.reactions.EX_pi_c.lower_bound = -10
    assert model.slim_optimize() == pytest.approx(0)


@pytest.mark.parametrize('change', ['formula', 'charge', 'compartment', 'stoichiometry', 'duplicate', 'equivalent', 'identity'])
def test_conflicts_rejected_before_parent_mutation(chemistry, change):
    base, source, candidate, mapping = chemistry
    if change in {'formula', 'charge', 'compartment'}:
        setattr(base.metabolites.adn_c, change, {'formula': 'H', 'charge': 1, 'compartment': 'external'}[change])
    elif change == 'stoichiometry':
        candidate['proposed_reactions'][0]['metabolites']['pi_c'] = -2
    elif change in {'duplicate', 'equivalent'}:
        r = source.reactions.PUNP1.copy()
        r.id = 'PUNP1' if change == 'duplicate' else 'existing_renamed'
        base.add_reactions([r])
    else:
        candidate['target_metabolite_map']['adn_c'] = 'ins_c'
    before = [R.reaction_record(r) for r in base.reactions]
    with pytest.raises(ValueError):
        applied(chemistry)
    assert [R.reaction_record(r) for r in base.reactions] == before


def test_source_parent_and_biomass_remain_unchanged(chemistry):
    base, source, candidate, _ = chemistry
    old_source = [R.reaction_record(r) for r in source.reactions]
    old_base = [R.reaction_record(r) for r in base.reactions]
    model = applied(chemistry)
    assert [R.reaction_record(r) for r in source.reactions] == old_source
    assert [R.reaction_record(r) for r in base.reactions] == old_base
    assert R.reaction_record(model.reactions.Growth) == R.reaction_record(base.reactions.Growth)
    assert {m.id: R.met_record(m) for m in model.metabolites} == candidate['expected_target_metabolites']


def test_failed_solve_is_not_zero_growth(chemistry, monkeypatch):
    model = applied(chemistry)
    monkeypatch.setattr(model, 'optimize', lambda: SimpleNamespace(status='undefined', objective_value=None))
    with pytest.raises(RuntimeError, match='Failed/nonfinite'):
        R.witness(model, {'physical_settings': {'residual_tolerance': 1e-8}, 'witness_reactions': []})
