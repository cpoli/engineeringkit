import numpy as np
import pytest

from engineeringkit import reliability as r


def test_rates_counts_cost_and_maintenance_limits():
    assert r.exponential_mttf(0.5) == 2
    assert r.exponential_hazard(0.5) == 0.5
    assert r.failure_rate_from_mtbf(10) == 0.1
    assert r.maintainability(0, 0.5) == 0
    assert np.isclose(r.maintainability(2, 0.5), 1 - 1 / np.e)
    assert r.repair_rate(2) == 0.5
    assert r.expected_failures(0.5, 2) == 1
    assert r.reliability_growth_duane(4, 2, 0.5) == 4
    assert r.bathtub_hazard(0) == 1.01
    assert r.safety_integrity_pfd(0.001, 10) == 0.005
    assert r.renewal_interval_cost_rate(10, 20, 0.1, 2) == 7
    with pytest.raises(ValueError):
        r.birnbaum_importance_series([0.9], 1)


def test_contract_errors_for_k_of_n_and_pfd():
    with pytest.raises(ValueError):
        r.k_of_n_reliability(2, 3, 0.9)
    with pytest.raises(ValueError):
        r.safety_integrity_pfd(1e-3, 1000)
