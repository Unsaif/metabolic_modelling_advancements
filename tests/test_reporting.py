from gembench.cards import carbon_fitness_leakage
from scripts.run_carbon_fitness_generic import load_model
import pytest


def test_patched_draft_discloses_fitness_reuse():
    card = carbon_fitness_leakage("MR1", "gapfilled", patched=True, medium_completion=True)
    assert card.ground_truth_used_in_model_curation.startswith("yes:")
    assert "No independent held-out" in card.notes[0]
    assert any("gene-less transport" in note for note in card.notes)


def test_gapfill_discloses_observed_growth():
    assert "observed wild-type growth" in carbon_fitness_leakage("MR1", "gapfilled", patched=False).ground_truth_used_in_model_curation


def test_unavailable_curated_arm_must_not_silently_load_a_draft():
    with pytest.raises(ValueError, match="No curated model"):
        load_model("Smeli", "curated")
