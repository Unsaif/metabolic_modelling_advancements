import cobra
import pytest

from scripts.run_quinone_biomass_sensitivity import apply_biomass_intervention


def model():
    result = cobra.Model('quinone')
    carbon, mql8, q8h2 = [cobra.Metabolite(name, compartment='c') for name in ('carbon_c', 'mql8_c', 'q8h2_c')]
    result.add_metabolites([mql8, q8h2])
    biomass = cobra.Reaction('Growth')
    biomass.add_metabolites({carbon: -1})
    result.add_reactions([biomass])
    return result


@pytest.mark.parametrize('arm,metabolite,expected', [
    ('no_quinone_demand', None, 0),
    ('ubiquinol_template_amount', 'q8h2_c', -0.0001),
    ('ubiquinol_curated_amount', 'q8h2_c', -0.000223),
    ('menaquinol_restored_control', 'mql8_c', -0.0001),
])
def test_absent_biomass_substrate_can_be_added_without_changing_other_components(arm, metabolite, expected):
    instance = model()
    coefficient = apply_biomass_intervention(instance, {'id': arm, 'metabolite': metabolite}, -0.0001, -0.000223)
    assert coefficient == expected
    composition = {m.id: value for m, value in instance.reactions.Growth.metabolites.items()}
    assert composition == ({'carbon_c': -1, metabolite: expected} if metabolite else {'carbon_c': -1})


def test_refuses_to_overwrite_an_existing_quinone_requirement():
    instance = model()
    instance.reactions.Growth.add_metabolites({instance.metabolites.mql8_c: -0.001})
    with pytest.raises(ValueError, match='no biomass quinone demand'):
        apply_biomass_intervention(instance, {'id': 'ubiquinol_template_amount', 'metabolite': 'q8h2_c'}, -0.0001, -0.000223)
