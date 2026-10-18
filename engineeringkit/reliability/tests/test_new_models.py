"""Analytical and Monte Carlo checks for the models added in the reliability rebuild."""

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import norm

from engineeringkit import reliability as r


def test_life_table_hazard_recovers_a_constant_hazard():
    ages = np.array([0.0, 1.0, 3.0, 6.0])
    survivors = 1000 * np.exp(-0.2 * ages)
    np.testing.assert_allclose(r.life_table_hazard(ages, survivors), 0.2)


@pytest.mark.parametrize(
    "ages,survivors",
    [([0, 1], [10, 5, 2]), ([0, 0], [10, 5]), ([0, 1], [5, 10]), ([0, 1], [10, 0]), ([0], [10])],
)
def test_life_table_hazard_rejects_bad_tables(ages, survivors):
    with pytest.raises(ValueError):
        r.life_table_hazard(ages, survivors)


def test_kaplan_meier_without_censoring_is_the_empirical_survival_function():
    times = np.array([5.0, 1.0, 3.0, 3.0, 8.0])
    km = r.kaplan_meier(times)
    np.testing.assert_array_equal(km.times, [1, 3, 5, 8])
    np.testing.assert_allclose(km.survival, [4 / 5, 2 / 5, 1 / 5, 0])
    np.testing.assert_array_equal(km.at_risk, [5, 4, 2, 1])
    np.testing.assert_array_equal(km.failures, [1, 2, 1, 1])
    np.testing.assert_allclose(km.survival_at([0.5, 3.0, 4.0, 100.0]), [1, 0.4, 0.4, 0])
    # Greenwood's formula reduces to the binomial standard error without censoring.
    np.testing.assert_allclose(km.std_error[:3], np.sqrt(km.survival[:3] * (1 - km.survival[:3]) / 5))
    assert km.std_error[-1] == 0
    assert km.summary() == {"n": 5, "failures": 5, "censored": 0, "final_survival": 0.0}


def test_kaplan_meier_hand_worked_censored_example():
    # A unit censored at a failure time is still at risk for that failure.
    km = r.kaplan_meier([2, 3, 3, 5, 6, 7], [1, 1, 0, 0, 1, 0])
    np.testing.assert_array_equal(km.at_risk, [6, 5, 2])
    expected = np.cumprod([5 / 6, 4 / 5, 1 / 2])
    np.testing.assert_allclose(km.survival, expected)
    greenwood = expected * np.sqrt(np.cumsum([1 / (6 * 5), 1 / (5 * 4), 1 / (2 * 1)]))
    np.testing.assert_allclose(km.std_error, greenwood)
    assert km.summary()["censored"] == 3


def test_kaplan_meier_all_censored_and_plot():
    km = r.kaplan_meier([1.0, 2.0], [0, 0])
    assert km.survival.size == 0
    assert km.summary()["final_survival"] == 1.0
    fig, ax = r.kaplan_meier([1.0, 2.0, 4.0]).plot()
    assert ax.get_ylabel() == "survival probability"
    fig.clf()


@pytest.mark.parametrize("failed", [[1, 0], [1, 2, 0], [[1, 0, 1]]])
def test_kaplan_meier_rejects_bad_flags(failed):
    with pytest.raises(ValueError):
        r.kaplan_meier([1.0, 2.0, 3.0], failed)


def test_mtbf_bounds_cover_the_true_mtbf_at_the_stated_rate():
    rng = np.random.default_rng(1953)
    theta, n, failures, hits = 1000.0, 10, 4, 0
    trials = 4000
    for _ in range(trials):
        life = np.sort(rng.exponential(theta, n))
        total = life[:failures].sum() + (n - failures) * life[failures - 1]
        low, high = r.mtbf_confidence_bounds(total, failures, 0.9)
        hits += low <= theta <= high
    assert abs(hits / trials - 0.9) < 0.02


def test_mtbf_bounds_time_terminated_and_errors():
    low, high = r.mtbf_confidence_bounds(1000.0, 0, 0.9, time_terminated=True)
    assert math.isclose(low, 1000 / -math.log(0.05)) and high == math.inf
    low_t, _ = r.mtbf_confidence_bounds(1000.0, 3, time_terminated=True)
    low_f, _ = r.mtbf_confidence_bounds(1000.0, 3)
    assert low_t < low_f
    with pytest.raises(ValueError):
        r.mtbf_confidence_bounds(1000.0, 0)
    with pytest.raises(ValueError):
        r.mtbf_confidence_bounds(1000.0, 2, confidence=1.0)


def test_weibull_mttf_is_vectorized_over_the_shape():
    beta = np.array([0.5, 1.0, 2.0])
    np.testing.assert_allclose(r.weibull_mttf(3.0, beta), [3.0 * math.gamma(1 + 1 / b) for b in beta])


def test_stress_strength_matches_monte_carlo():
    rng = np.random.default_rng(1947)
    strength, stress = rng.normal(500, 30, 400_000), rng.normal(400, 40, 400_000)
    assert np.isclose(r.stress_strength_reliability(500, 30, 400, 40), np.mean(strength > stress), atol=2e-3)
    assert np.isclose(r.stress_strength_reliability(500, 30, 400, 40), norm.cdf(2.0))
    with pytest.raises(ValueError):
        r.stress_strength_reliability(500, 0, 400, 0)


def test_normal_index_conversions_round_trip():
    beta = np.array([-1.0, 0.0, 3.0, 5.0])
    np.testing.assert_allclose(r.reliability_index_from_probability(r.failure_probability_from_index(beta)), beta)
    with pytest.raises(ValueError):
        r.reliability_index_from_probability(0.0)


def test_hasofer_lind_linear_limit_state_is_exact():
    result = r.hasofer_lind_index(lambda x: x[0] - x[1], [500, 400], [30, 40])
    assert np.isclose(result.beta, 2.0)
    np.testing.assert_allclose(result.alpha, [-0.6, 0.8], atol=1e-8)
    assert np.isclose(result.design_point[0], result.design_point[1])
    assert np.isclose(result.failure_probability, 1 - r.stress_strength_reliability(500, 30, 400, 40))
    assert result.summary()["iterations"] == result.iterations


def test_hasofer_lind_is_invariant_to_the_limit_state_formulation():
    linear = r.hasofer_lind_index(lambda x: x[0] - x[1], [500, 400], [30, 40])
    ratio = r.hasofer_lind_index(lambda x: x[0] / x[1] - 1, [500, 400], [30, 40])
    assert np.isclose(linear.beta, ratio.beta, rtol=1e-6)


def test_hasofer_lind_nonlinear_matches_known_design_point():
    # g = u1 + u2**2/4 ... written in physical variables; the nearest point is analytic.
    result = r.hasofer_lind_index(lambda x: 3 - x[0] - 0.1 * x[1] ** 2, [0, 0], [1, 1])
    u = result.design_point
    assert np.isclose(3 - u[0] - 0.1 * u[1] ** 2, 0, atol=1e-8)
    assert result.beta <= 3


def test_hasofer_lind_negative_index_and_errors():
    assert r.hasofer_lind_index(lambda x: x[0] - x[1], [400, 500], [30, 40]).beta < 0
    with pytest.raises(ValueError):
        r.hasofer_lind_index("g", [1], [1])
    with pytest.raises(ValueError):
        r.hasofer_lind_index(lambda x: x[0], [1, 2], [1])
    with pytest.raises(ValueError):
        r.hasofer_lind_index(lambda x: 1.0, [1], [1])
    with pytest.raises(ValueError):
        r.hasofer_lind_index(lambda x: 3 - x[0] - 0.1 * x[1] ** 2, [0, 0], [1, 1], max_iter=1)


def test_fault_tree_gate_algebra_and_repeated_events():
    p = {"A": 0.1, "B": 0.2, "C": 0.3}
    assert np.isclose(r.fault_tree_probability(("or", "A", ("and", "B", "C")), p), 1 - 0.9 * (1 - 0.06))
    assert r.fault_tree_probability("A", p) == 0.1
    shared = ("and", ("or", "A", "B"), ("or", "A", "C"))
    assert np.isclose(r.fault_tree_probability(shared, p), 0.1 + 0.9 * 0.2 * 0.3)


def test_fault_tree_matches_monte_carlo():
    rng = np.random.default_rng(1961)
    p = {"A": 0.3, "B": 0.4, "C": 0.5, "D": 0.2}
    tree = ("or", ("and", "A", "B"), ("and", "B", "C", "D"), ("and", "A", "D"))
    s = {name: rng.random(400_000) < q for name, q in p.items()}
    top = (s["A"] & s["B"]) | (s["B"] & s["C"] & s["D"]) | (s["A"] & s["D"])
    assert np.isclose(r.fault_tree_probability(tree, p), top.mean(), atol=3e-3)


@pytest.mark.parametrize(
    "tree,p",
    [
        (("xor", "A", "B"), {"A": 0.1, "B": 0.1}),
        (("and",), {}),
        (["and", "A", "B"], {"A": 0.1, "B": 0.1}),
        (("and", "A", "B"), {"A": 0.1}),
        (("and", "A", "B"), [0.1, 0.1]),
        (("and", "A", "B"), {"A": 0.1, "B": 1.5}),
        (("and", "E0", *[f"E{i}" for i in range(21)]), {f"E{i}": 0.1 for i in range(21)}),
    ],
)
def test_fault_tree_rejects_bad_input(tree, p):
    with pytest.raises(ValueError):
        r.fault_tree_probability(tree, p)


def test_event_tree_sequences_partition_the_initiator():
    result = r.event_tree_probability(1e-3, [0.01, 0.1, 0.5])
    assert result.outcomes.shape == (8, 3)
    assert not result.outcomes[0].any() and result.outcomes[-1].all()
    assert np.isclose(result.frequencies.sum(), 1e-3)
    assert np.isclose(result.frequency([0]), 1e-3 * 0.01)
    assert np.isclose(result.frequency([1, 2]), 1e-3 * 0.05)
    assert np.isclose(result.frequency([]), 1e-3)
    assert result.summary() == {
        "sequences": 8,
        "total_frequency": pytest.approx(1e-3),
        "all_success_frequency": pytest.approx(1e-3 * 0.99 * 0.9 * 0.5),
    }
    with pytest.raises(ValueError):
        result.frequency([3])
    with pytest.raises(ValueError):
        r.event_tree_probability(1.0, [0.1] * 21)


def test_safety_integrity_level_bands():
    assert r.safety_integrity_level(8.76e-3) == 2
    assert r.safety_integrity_level(1e-3) == 2
    np.testing.assert_array_equal(r.safety_integrity_level([1e-6, 5e-4, 0.05, 0.5]), [4, 3, 1, 0])
    with pytest.raises(ValueError):
        r.safety_integrity_level(1.5)


def test_safety_integrity_pfd_is_vectorized_and_near_the_exact_average():
    lam, T = 2e-6, np.array([100.0, 8760.0, 40000.0])
    exact = 1 - (1 - np.exp(-lam * T)) / (lam * T)
    np.testing.assert_allclose(r.safety_integrity_pfd(lam, T), exact, rtol=0.03)
    assert isinstance(r.safety_integrity_pfd(lam, 8760.0), float)
    with pytest.raises(ValueError):
        r.safety_integrity_pfd(lam, np.array([8760.0, 1e5]))


def test_markov_availability_matches_the_generator_solution():
    lam, mu = 1e-3, 0.05
    t = np.array([0.0, 10.0, 100.0, 1e4])
    Q = [[-lam, lam], [mu, -mu]]
    p = r.markov_state_probabilities(Q, [1, 0], t)
    assert p.shape == (4, 2)
    np.testing.assert_allclose(r.markov_availability(t, lam, mu), p[:, 0])
    assert np.isclose(r.markov_availability(1e9, lam, mu), r.availability(1 / lam, 1 / mu))
    np.testing.assert_allclose(r.markov_state_probabilities(Q, [0.5, 0.5], 0.0), [0.5, 0.5])
    with pytest.raises(ValueError):
        r.markov_availability(1.0, 0.0, 0.0)


@pytest.mark.parametrize(
    "Q,p0",
    [
        ([1.0, -1.0], [1.0]),
        ([[-1.0, 1.0], [-1.0, 1.0]], [1, 0]),
        ([[-1.0, 0.5], [1.0, -1.0]], [1, 0]),
        ([[-1.0, 1.0], [1.0, -1.0]], [1, 0, 0]),
        ([[-1.0, 1.0], [1.0, -1.0]], [0.5, 0.6]),
        ([[0.0, 0.0], [0.0, 0.0]], [1.0, -0.0001]),
    ],
)
def test_markov_state_probabilities_rejects_bad_input(Q, p0):
    with pytest.raises(ValueError):
        r.markov_state_probabilities(Q, p0, 1.0)


def test_age_replacement_cost_rate_against_quadrature_and_optimality():
    cp, cf, eta, beta = 1.0, 10.0, 1000.0, 2.5

    def direct(T):
        R = r.weibull_reliability(T, eta, beta)
        return (cp * R + cf * (1 - R)) / quad(lambda t: r.weibull_reliability(t, eta, beta), 0, T)[0]

    for T in (100.0, 500.0, 2000.0):
        assert np.isclose(r.age_replacement_cost_rate(T, cp, cf, eta, beta), direct(T))
    grid = np.linspace(50, 3000, 4000)
    T_star = grid[np.argmin(r.age_replacement_cost_rate(grid, cp, cf, eta, beta))]
    # Barlow-Proschan optimality: h(T) * integral(R) - F(T) = cp / (cf - cp).
    lhs = r.weibull_hazard(T_star, eta, beta) * quad(lambda t: r.weibull_reliability(t, eta, beta), 0, T_star)[0]
    assert np.isclose(lhs - r.weibull_cdf(T_star, eta, beta), cp / (cf - cp), rtol=1e-2)


def test_age_replacement_never_pays_for_a_constant_hazard():
    T = np.geomspace(1, 1e6, 50)
    cost = r.age_replacement_cost_rate(T, 1.0, 10.0, 1000.0, 1.0)
    assert np.all(np.diff(cost) <= 0) and cost[0] > 100 * cost[-1]
    assert np.isclose(cost[-1], r.renewal_interval_cost_rate(0.0, 10.0, 1e-3, 1.0))
