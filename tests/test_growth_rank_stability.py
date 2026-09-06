import numpy as np
import pytest

from scripts.audit_growth_rank_stability import transformed
from scripts.verify_quinone_repair import metric_points


def test_near_zero_ordering_changes_ap_without_changing_growth_classification():
    fit = np.array([-3., 0., -3., 0.])
    left = np.array([-1e-11, 1e-11, -1e-12, 1e-12])
    right = -left
    a, b = (metric_points(v, fit, .001, -2) for v in (left, right))
    assert a['aucpr_standard'] != b['aucpr_standard']
    assert np.array_equal(left >= .001, right >= .001)
    assert np.array_equal(transformed(left, 1e-9, 'zero_band'), transformed(right, 1e-9, 'zero_band'))
    assert metric_points(transformed(left, 1e-9, 'zero_band'), fit, .001, -2)['aucpr_standard'] == pytest.approx(.5)


def test_precision_policy_preserves_missing_values_and_original_array():
    original = np.array([np.nan, np.inf, -np.inf, 1e-11, 2.])
    before = original.copy()
    for policy in ('zero_band', 'round_to_quantum'):
        result = transformed(original, 1e-9, policy)
        assert np.array_equal(result[:3], before[:3], equal_nan=True)
        assert result[3] == 0 and result[4] == 2
        assert np.array_equal(original, before, equal_nan=True)
