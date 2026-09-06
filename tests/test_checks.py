from types import SimpleNamespace

import cobra
import pandas as pd
import pytest

from gembench.checks import energy_from_nothing
from gembench.protocols.carbon_fitness_generic import _growth_by_gene, run


def energy_model():
    model = cobra.Model("energy")
    mets = {mid: cobra.Metabolite(mid, compartment="c")
            for mid in ("atp_c", "adp_c", "pi_c", "h_c", "h2o_c")}
    reaction = cobra.Reaction("FREE_ATP")
    reaction.add_metabolites({mets["atp_c"]: 1, mets["h2o_c"]: 1,
                             mets["adp_c"]: -1, mets["pi_c"]: -1, mets["h_c"]: -1})
    model.add_reactions([reaction])
    model.objective = reaction
    return model


def test_mandatory_flux_cannot_hide_energy_cycle_and_model_is_restored():
    model = energy_model()
    forced = cobra.Reaction("FORCED", lower_bound=1)
    forced.add_metabolites({cobra.Metabolite("x_c", compartment="c"): -1,
                            cobra.Metabolite("y_c", compartment="c"): 1})
    model.add_reactions([forced])
    before = {r.id: r.bounds for r in model.reactions}
    model.objective_direction = "min"
    assert energy_from_nothing(model, currencies=("atp",))["atp"] == pytest.approx(1000)
    assert {r.id: r.bounds for r in model.reactions} == before
    assert model.objective_direction == "min"


def test_source_capable_sinks_are_closed_for_energy_diagnostic():
    model = energy_model()
    model.reactions.FREE_ATP.bounds = (0, 0)
    for met in list(model.metabolites):
        sink = cobra.Reaction(f"arbitrary_source_{met.id}", lower_bound=-1000)
        sink.add_metabolites({met: -1})
        model.add_reactions([sink])
    assert energy_from_nothing(model, currencies=("atp",)) == {"atp": 0}


def test_failed_diagnostic_is_not_a_clean_zero():
    model = energy_model()
    model.add_cons_vars(model.problem.Constraint(0, lb=1, name="impossible"))
    with pytest.raises(RuntimeError, match="Energy-cycle check.*infeasible"):
        energy_from_nothing(model, currencies=("atp",))


def test_existing_diagnostic_id_does_not_silently_replace_probe():
    model = energy_model()
    existing = cobra.Reaction("EGC_atp", upper_bound=0)
    existing.add_metabolites({model.metabolites.atp_c: -1})
    model.add_reactions([existing])
    assert energy_from_nothing(model, currencies=("atp",))["atp"] == pytest.approx(1000)
    assert model.reactions.EGC_atp.bounds == (0, 0)


def test_protocol_stops_on_energy_cycle_before_scoring():
    org = SimpleNamespace(org_id="test", fitness=pd.DataFrame())
    gene_map = SimpleNamespace(model_to_browser={})
    with pytest.raises(ValueError, match="Energy-generating cycle"):
        run(energy_model(), org, [], gene_map, verbose=False)


@pytest.mark.parametrize("status", ["time_limit", "undefined", "optimal"])
def test_knockout_failure_or_nonfinite_optimum_is_not_essentiality(status):
    result = pd.DataFrame({"ids": [{"g"}], "growth": [float("nan")], "status": [status]})
    with pytest.raises(RuntimeError, match="Growth optimization failed"):
        _growth_by_gene(result, ["g"])


def test_infeasible_knockout_counts_as_no_growth_but_missing_result_fails():
    result = pd.DataFrame({"ids": [{"g"}], "growth": [float("nan")], "status": ["infeasible"]})
    assert _growth_by_gene(result, ["g"]) == {"g": 0}
    with pytest.raises(RuntimeError, match="Missing gene deletion results"):
        _growth_by_gene(result, ["g", "missing"])
