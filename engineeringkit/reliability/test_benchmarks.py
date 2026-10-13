import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import binom, poisson, weibull_min

from engineeringkit import reliability as r


@pytest.mark.parametrize("beta", [0.5, 1, 2, 4])
def test_weibull_against_distribution_and_survival_integral(beta):
    eta = 3.0
    times = np.linspace(0, 10, 21)
    np.testing.assert_allclose(r.weibull_reliability(times, eta, beta), weibull_min.sf(times, beta, scale=eta))
    numerical = quad(lambda t: r.weibull_reliability(t, eta, beta), 0, np.inf)[0]
    assert np.isclose(r.weibull_mttf(eta, beta), numerical, rtol=1e-8)


def test_hazard_shape_and_zero_time_limits():
    assert np.isinf(r.weibull_hazard(0, 2, 0.5))
    assert r.weibull_hazard(0, 2, 1) == 0.5
    assert r.weibull_hazard(0, 2, 2) == 0
    assert r.exponential_failure_probability(1e-12, 1) > 0
    assert r.weibull_cdf(1e-12, 1, 1) > 0


def test_system_reliability_by_exhaustive_states():
    probs = [0.7, 0.8, 0.9]
    states = [tuple((code >> i) & 1 for i in range(3)) for code in range(8)]
    masses = [np.prod([p if state[i] else 1 - p for i, p in enumerate(probs)]) for state in states]
    series = sum(mass for state, mass in zip(states, masses) if all(state))
    parallel = sum(mass for state, mass in zip(states, masses) if any(state))
    assert np.isclose(r.series_reliability(*probs), series)
    assert np.isclose(r.parallel_reliability(*probs), parallel)
    assert np.isclose(r.k_of_n_reliability(5, 3, 0.8), binom.sf(2, 5, 0.8))


@pytest.mark.parametrize("k,mean", [(0, 0), (0, 2), (3, 2), (1000, 1000)])
def test_poisson_extreme_counts(k, mean):
    assert np.isclose(r.poisson_failures_probability(k, mean, 1), poisson.pmf(k, mean))


def test_availability_and_birnbaum_sensitivity():
    assert r.availability(100, 0) == 1
    assert r.availability(100, 25) + r.unavailability(100, 25) == 1
    x = [0.7, 0.8, 0.9]
    eps = 1e-5
    left = x.copy()
    right = x.copy()
    left[1] -= eps
    right[1] += eps
    derivative = (r.series_reliability(*right) - r.series_reliability(*left)) / (2 * eps)
    assert np.isclose(derivative, r.birnbaum_importance_series(x, 1))
