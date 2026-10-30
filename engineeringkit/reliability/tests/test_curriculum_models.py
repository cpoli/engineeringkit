"""Analytical, distribution and Monte Carlo checks for the models added by ROADMAP item D9.

Lognormal and gamma lives, median ranks and Weibull fitting, accelerated
life testing, block-diagram cut sets, cold standby, FMEA, Monte Carlo
failure probability and Bayesian updating of a failure rate.
"""

import itertools
import math

import matplotlib.pyplot as plt
import numpy as np
import pytest
from scipy import optimize, stats
from scipy.constants import eV, k
from scipy.integrate import quad

from engineeringkit import reliability as r

BRIDGE = {"A": ("s", 1), "B": ("s", 2), "C": (1, 2), "D": (1, "t"), "E": (2, "t")}


# Lognormal and gamma lifetimes ------------------------------------------------


def test_lognormal_matches_scipy_and_its_mean_integral():
    mu, sigma = np.log(2000.0), 0.6
    t = np.linspace(0, 20000, 41)
    dist = stats.lognorm(sigma, scale=np.exp(mu))
    np.testing.assert_allclose(r.lognormal_reliability(t, mu, sigma), dist.sf(t), atol=1e-15)
    np.testing.assert_allclose(r.lognormal_cdf(t, mu, sigma), dist.cdf(t), atol=1e-15)
    positive = t[1:]
    np.testing.assert_allclose(r.lognormal_hazard(positive, mu, sigma), dist.pdf(positive) / dist.sf(positive), rtol=1e-9)
    assert r.lognormal_hazard(0.0, mu, sigma) == 0
    assert r.lognormal_reliability(0.0, mu, sigma) == 1
    assert np.isclose(r.lognormal_mttf(mu, sigma), quad(lambda x: r.lognormal_reliability(x, mu, sigma), 0, np.inf)[0], rtol=1e-7)
    assert np.isclose(r.lognormal_mttf(mu, sigma), dist.mean())


def test_lognormal_hazard_stays_finite_far_into_the_tail():
    # Where the survival underflows the hazard tends to z/(sigma t) (Mills ratio).
    t, sigma = np.exp(40 * 0.5), 0.5
    z = 40.0
    h = r.lognormal_hazard(t, 0.0, sigma)
    assert np.isfinite(h) and np.isclose(h, z / (sigma * t), rtol=1e-2)
    assert r.lognormal_reliability(t, 0.0, sigma) == 0


@pytest.mark.parametrize("shape", [0.5, 1.0, 2.5, 4.0])
def test_gamma_matches_scipy_and_its_mean_integral(shape):
    scale = 300.0
    t = np.linspace(1, 5000, 50)
    dist = stats.gamma(shape, scale=scale)
    np.testing.assert_allclose(r.gamma_reliability(t, shape, scale), dist.sf(t), rtol=1e-12)
    np.testing.assert_allclose(r.gamma_cdf(t, shape, scale), dist.cdf(t), rtol=1e-12)
    np.testing.assert_allclose(r.gamma_hazard(t, shape, scale), dist.pdf(t) / dist.sf(t), rtol=1e-9)
    assert np.isclose(r.gamma_mttf(shape, scale), quad(lambda x: r.gamma_reliability(x, shape, scale), 0, np.inf)[0], rtol=1e-7)


def test_gamma_hazard_limits():
    assert np.isinf(r.gamma_hazard(0.0, 0.5, 2.0))
    assert r.gamma_hazard(0.0, 1.0, 2.0) == 0.5
    assert r.gamma_hazard(0.0, 3.0, 2.0) == 0
    # Every shape tends to the stage rate 1/scale.
    assert np.isclose(r.gamma_hazard(200.0, 3.0, 2.0), 0.5, rtol=0.02)


# Median ranks and Weibull fitting ---------------------------------------------


def test_median_ranks_are_beta_medians_and_benard_is_close():
    n = 12
    exact = r.median_ranks(n)
    np.testing.assert_allclose(exact, [stats.beta.median(i, n - i + 1) for i in range(1, n + 1)])
    assert np.max(np.abs(r.median_ranks(n, method="benard") - exact)) < 0.003
    assert r.median_ranks(5)[2] == pytest.approx(0.5)
    np.testing.assert_allclose(r.median_ranks(4, [1.5], "benard"), [1.2 / 4.4])


@pytest.mark.parametrize("args", [(0,), (5, [0.5]), (5, [6]), (5, None, "mean")])
def test_median_ranks_reject_bad_input(args):
    with pytest.raises(ValueError):
        r.median_ranks(*args)


def test_adjusted_ranks_without_censoring_and_hand_example():
    times, order = r.adjusted_ranks([30, 10, 20])
    np.testing.assert_array_equal(times, [10, 20, 30])
    np.testing.assert_allclose(order, [1, 2, 3])
    # n = 5: F, S, F, S, F. Order numbers 1, then (3*1 + 6)/4 = 2.25, then (1*2.25 + 6)/2 = 4.125.
    times, order = r.adjusted_ranks([1, 2, 3, 4, 5], [1, 0, 1, 0, 1])
    np.testing.assert_allclose(order, [1, 2.25, 4.125])
    # At a tie the failure is ranked before the censored unit (else its order would be 4/3).
    _, order = r.adjusted_ranks([2, 2, 3], [0, 1, 1])
    np.testing.assert_allclose(order, [1, 2.5])


@pytest.mark.parametrize("times,failed", [([1, 2], [1]), ([1, 2], [1, 2]), ([0, 1], None), ([-1, 2], [1, 1])])
def test_life_data_rejects_bad_input(times, failed):
    with pytest.raises(ValueError):
        r.adjusted_ranks(times, failed)


@pytest.mark.parametrize("regress", ["y", "x"])
def test_rank_regression_recovers_points_on_the_line(regress):
    eta, beta, n = 1500.0, 2.2, 10
    f = r.median_ranks(n)
    times = eta * (-np.log1p(-f)) ** (1 / beta)
    fit = r.weibull_rank_regression(times, regress=regress)
    assert np.isclose(fit.eta, eta) and np.isclose(fit.beta, beta)
    assert fit.method == f"rank regression on {regress}"
    assert np.isclose(fit.b_life(0.1), r.weibull_reliability(0, 1, 1) * eta * (-math.log(0.9)) ** (1 / beta))


def test_rank_regression_with_censoring_and_errors():
    rng = np.random.default_rng(1951)
    life = 1000 * rng.weibull(3.0, 400)
    censor = rng.uniform(0, 1600, 400)
    fit = r.weibull_rank_regression(np.minimum(life, censor), (life <= censor).astype(int))
    assert abs(fit.beta - 3) < 0.4 and abs(fit.eta - 1000) < 60
    assert fit.units == 400 and fit.failure_times.size == np.count_nonzero(life <= censor)
    with pytest.raises(ValueError):
        r.weibull_rank_regression([1, 2], regress="z")
    with pytest.raises(ValueError):
        r.weibull_rank_regression([1, 2, 3], [1, 0, 0])
    with pytest.raises(ValueError):
        r.weibull_rank_regression([2, 2, 3], [1, 1, 0])


@pytest.mark.parametrize("beta", [0.6, 3.5])
def test_weibull_mle_matches_scipy_on_complete_data(beta):
    times = 200 * np.random.default_rng(7).weibull(beta, 30)
    fit = r.weibull_mle(times)
    shape, _, scale = stats.weibull_min.fit(times, floc=0)
    assert np.isclose(fit.beta, shape, rtol=1e-5) and np.isclose(fit.eta, scale, rtol=1e-5)


def test_weibull_mle_maximises_the_censored_likelihood():
    rng = np.random.default_rng(1965)
    life = 500 * rng.weibull(1.8, 60)
    t = np.minimum(life, 450.0)
    d = life <= 450.0

    def negative_log_likelihood(p):
        eta, beta = np.exp(p)
        z = (t / eta) ** beta
        return -(np.sum(np.log(beta / eta) + (beta - 1) * np.log(t[d] / eta)) - np.sum(z))

    best = optimize.minimize(negative_log_likelihood, [np.log(400), 0.0], method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-12})
    fit = r.weibull_mle(t, d.astype(int))
    assert np.isclose(fit.eta, np.exp(best.x[0]), rtol=1e-4) and np.isclose(fit.beta, np.exp(best.x[1]), rtol=1e-4)
    assert fit.summary()["failures"] == np.count_nonzero(d)


def test_weibull_mle_rejects_degenerate_data():
    with pytest.raises(ValueError):
        r.weibull_mle([1, 2], [0, 0])
    with pytest.raises(ValueError):
        r.weibull_mle([3, 5, 5], [0, 1, 1])


def test_weibull_fit_result_summary_plot_and_b_life():
    fit = r.weibull_mle([16, 34, 53, 75, 93, 120], ranks="benard")
    summary = fit.summary()
    assert summary["n"] == 6 and summary["failures"] == 6
    assert np.isclose(r.weibull_cdf(summary["b10_life"], fit.eta, fit.beta), 0.1)
    with pytest.raises(ValueError):
        fit.b_life(1.0)
    fig, ax = fit.plot()
    assert ax.get_title() == "Weibull probability plot"
    plt.close(fig)


# Accelerated life testing -----------------------------------------------------


def test_arrhenius_factor_identities():
    ea = 0.7 * eV
    assert r.arrhenius_acceleration_factor(ea, 330.0, 330.0) == 1
    assert r.arrhenius_acceleration_factor(0.0, 300.0, 400.0) == 1
    chained = r.arrhenius_acceleration_factor(ea, 300.0, 350.0) * r.arrhenius_acceleration_factor(ea, 350.0, 400.0)
    assert np.isclose(r.arrhenius_acceleration_factor(ea, 300.0, 400.0), chained)
    assert np.isclose(r.arrhenius_acceleration_factor(ea, 400.0, 300.0), 1 / r.arrhenius_acceleration_factor(ea, 300.0, 400.0))
    with pytest.raises(ValueError):
        r.arrhenius_acceleration_factor(-ea, 300.0, 400.0)


def test_arrhenius_fit_recovers_the_model_and_predicts_use_life():
    ea, constant = 0.85 * eV, 2e-7
    temperatures = np.array([373.15, 398.15, 423.15])
    lives = constant * np.exp(ea / (k * temperatures))
    fit_ea, fit_constant = r.arrhenius_fit(temperatures, lives)
    assert np.isclose(fit_ea, ea) and np.isclose(fit_constant, constant)
    use_life = lives[0] * r.arrhenius_acceleration_factor(fit_ea, 328.15, 373.15)
    assert np.isclose(use_life, constant * np.exp(ea / (k * 328.15)))
    for bad in (([300.0, 300.0], [1.0, 2.0]), ([300.0, 310.0, 320.0], [1.0, 2.0]), ([300.0, 310.0], [1.0, -2.0])):
        with pytest.raises(ValueError):
            r.arrhenius_fit(*bad)


def test_eyring_and_inverse_power_law():
    ea = 0.5 * eV
    assert np.isclose(r.eyring_acceleration_factor(ea, 300.0, 360.0), 1.2 * r.arrhenius_acceleration_factor(ea, 300.0, 360.0))
    # Eyring's life L = (C/T) exp(E/kT) gives the factor directly.
    life = lambda temp: np.exp(ea / (k * temp)) / temp  # noqa: E731
    assert np.isclose(r.eyring_acceleration_factor(ea, 300.0, 360.0), life(300.0) / life(360.0))
    np.testing.assert_allclose(r.inverse_power_law_acceleration_factor(100.0, [100.0, 200.0, 400.0], 2.0), [1, 4, 16])
    with pytest.raises(ValueError):
        r.inverse_power_law_acceleration_factor(100.0, 200.0, -1.0)


# Block diagrams and standby ---------------------------------------------------


def _connected(components, up, source, sink):
    """Independent check: is sink reachable from source over working components (flood fill)?"""
    reached, frontier = {source}, [source]
    while frontier:
        node = frontier.pop()
        for name, (a, b) in components.items():
            if name in up and node in (a, b):
                other = b if node == a else a
                if other not in reached:
                    reached.add(other)
                    frontier.append(other)
    return sink in reached


def test_bridge_reliability_matches_its_polynomial_and_state_enumeration():
    p = 0.9
    exact = 2 * p**2 + 2 * p**3 - 5 * p**4 + 2 * p**5
    assert np.isclose(r.network_reliability(BRIDGE, dict.fromkeys(BRIDGE, p), "s", "t"), exact)
    rel = {"A": 0.95, "B": 0.7, "C": 0.6, "D": 0.85, "E": 0.9}
    brute = 0.0
    for states in itertools.product([0, 1], repeat=5):
        up = {name for name, s in zip(BRIDGE, states) if s}
        weight = math.prod(rel[n] if n in up else 1 - rel[n] for n in BRIDGE)
        brute += weight * _connected(BRIDGE, up, "s", "t")
    assert np.isclose(r.network_reliability(BRIDGE, rel, "s", "t"), brute)
    lower, upper = r.esary_proschan_bounds(BRIDGE, rel, "s", "t")
    assert lower <= brute <= upper


def test_cut_sets_disconnect_and_are_minimal():
    cuts = r.minimal_cut_sets(BRIDGE, "s", "t")
    assert len(cuts) == 4
    for cut in cuts:
        assert not _connected(BRIDGE, set(BRIDGE) - set(cut), "s", "t")
        for name in cut:
            assert _connected(BRIDGE, set(BRIDGE) - set(cut) | {name}, "s", "t")


def test_series_and_parallel_are_special_networks():
    series = {"a": ("in", 1), "b": (1, 2), "c": (2, "out")}
    parallel = {"a": ("in", "out"), "b": ("in", "out"), "c": ("in", "out")}
    rel = {"a": 0.9, "b": 0.8, "c": 0.7}
    assert np.isclose(r.network_reliability(series, rel, "in", "out"), r.series_reliability(0.9, 0.8, 0.7))
    assert np.isclose(r.network_reliability(parallel, rel, "in", "out"), r.parallel_reliability(0.9, 0.8, 0.7))
    assert r.minimal_cut_sets(series, "in", "out") == [("a",), ("b",), ("c",)]
    assert r.minimal_path_sets(parallel, "in", "out") == [("a",), ("b",), ("c",)]
    # For a series system the lower bound is exact; for a parallel one the upper bound is.
    assert np.isclose(r.esary_proschan_bounds(series, rel, "in", "out")[0], r.series_reliability(0.9, 0.8, 0.7))
    assert np.isclose(r.esary_proschan_bounds(parallel, rel, "in", "out")[1], r.parallel_reliability(0.9, 0.8, 0.7))


def test_disconnected_network_and_errors():
    split = {"a": ("in", 1), "b": (2, "out")}
    assert r.minimal_path_sets(split, "in", "out") == []
    assert r.minimal_cut_sets(split, "in", "out") == []
    assert r.network_reliability(split, {"a": 1, "b": 1}, "in", "out") == 0
    assert r.esary_proschan_bounds(split, {"a": 1, "b": 1}, "in", "out") == (0.0, 0.0)
    for components in ({}, [("a", (1, 2))], {"a": (1, 2, 3)}, {1: (1, 2)}):
        with pytest.raises(ValueError):
            r.minimal_path_sets(components, 1, 2)
    with pytest.raises(ValueError):
        r.minimal_path_sets(BRIDGE, "s", "s")
    for rel in ([0.9] * 5, {"A": 0.9}, dict.fromkeys(BRIDGE, 1.2)):
        with pytest.raises(ValueError):
            r.network_reliability(BRIDGE, rel, "s", "t")
    wide = {f"c{i:02d}": ("s", "t") for i in range(21)}
    with pytest.raises(ValueError):
        r.minimal_cut_sets(wide, "s", "t")
    with pytest.raises(ValueError):
        r.network_reliability(wide, dict.fromkeys(wide, 0.5), "s", "t")


def test_cold_standby_is_a_gamma_life_and_beats_active_parallel():
    lam, t = 2e-4, np.linspace(0, 30000, 31)
    for n in (1, 2, 4):
        np.testing.assert_allclose(r.standby_reliability(t, lam, n), r.gamma_reliability(t, n, 1 / lam), atol=1e-12)
        integral = quad(lambda u: r.standby_reliability(u / lam, lam, n), 0, np.inf)[0] / lam
        assert np.isclose(r.standby_mttf(lam, n), integral, rtol=1e-7)
    active = 1 - (1 - r.exponential_reliability(t, lam)) ** 2
    assert np.all(r.standby_reliability(t, lam, 2) >= active - 1e-15)


def test_standby_with_an_imperfect_switch_matches_simulation():
    lam, p, mission = 1e-3, 0.8, 1500.0
    rng = np.random.default_rng(1)
    first = rng.exponential(1 / lam, 200_000)
    second = rng.exponential(1 / lam, 200_000)
    switched = rng.random(200_000) < p
    survive = (first > mission) | (switched & (first + second > mission))
    assert np.isclose(r.standby_reliability(mission, lam, 2, p), survive.mean(), atol=4e-3)
    assert np.isclose(r.standby_mttf(lam, 3, p), (1 + p + p**2) / lam)
    with pytest.raises(ValueError):
        r.standby_reliability(1.0, 1e-3, 2, 1.5)


# FMEA ---------------------------------------------------------------------------


def test_risk_priority_number_and_ranking():
    assert r.risk_priority_number(10, 10, 10) == 1000
    np.testing.assert_array_equal(r.risk_priority_number(np.array([1, 2]), 3, np.array([4, 5])), [12, 30])
    result = r.fmea_ranking({"b": (5, 4, 3), "a": (3, 5, 4), "c": (10, 1, 2), "d": (6, 2, 5)})
    assert result.modes == ("d", "b", "a", "c")  # d and b tie at 60, broken by severity
    np.testing.assert_array_equal(result.rpn, [60, 60, 60, 20])
    np.testing.assert_array_equal(result.rpn, result.severity * result.occurrence * result.detection)
    assert result.summary() == {"modes": 4, "top_mode": "d", "top_rpn": 60, "total_rpn": 200}
    fig, ax = result.plot()
    assert [t.get_text() for t in ax.get_xticklabels()] == list(result.modes)
    plt.close(fig)


@pytest.mark.parametrize("ranks", [(0, 5, 5), (11, 5, 5), (2.5, 5, 5), (True, 5, 5)])
def test_risk_priority_number_rejects_bad_rankings(ranks):
    with pytest.raises(ValueError):
        r.risk_priority_number(*ranks)


@pytest.mark.parametrize("modes", [{}, [("a", (1, 1, 1))], {"a": (1, 1)}, {1: (1, 1, 1)}])
def test_fmea_ranking_rejects_bad_tables(modes):
    with pytest.raises(ValueError):
        r.fmea_ranking(modes)


def test_failure_mode_criticality_sums_to_expected_losses():
    modes = np.array([0.5, 0.3, 0.2])
    criticality = r.failure_mode_criticality(modes, 1.0, 1e-5, 100.0)
    assert np.isclose(criticality.sum(), 1e-5 * 100)
    with pytest.raises(ValueError):
        r.failure_mode_criticality(1.2, 0.5, 1e-5, 1.0)


# Monte Carlo next to Hasofer-Lind ---------------------------------------------


def test_crude_monte_carlo_agrees_with_the_exact_linear_result():
    mean, std = [500.0, 300.0], [50.0, 60.0]
    exact = r.hasofer_lind_index(lambda x: x[0] - x[1], mean, std).failure_probability
    result = r.monte_carlo_failure_probability(lambda x: x[0] - x[1], mean, std, samples=400_000, seed=1949)
    assert abs(result.failure_probability - exact) < 4 * result.std_error
    assert np.isclose(result.std_error, np.sqrt(exact * (1 - exact) / 400_000), rtol=0.05)
    low, high = result.confidence_interval(0.999)
    assert low < exact < high
    assert np.isclose(result.beta, r.reliability_index_from_probability(result.failure_probability))
    assert not result.importance_sampling and result.summary()["samples"] == 400_000


def test_importance_sampling_at_the_design_point_is_far_cheaper():
    mean, std = [500.0, 200.0], [40.0, 50.0]

    def g(x):
        return x[0] - x[1]

    hl = r.hasofer_lind_index(g, mean, std)
    crude = r.monte_carlo_failure_probability(g, mean, std, samples=20_000, seed=3)
    weighted = r.monte_carlo_failure_probability(g, mean, std, samples=20_000, seed=3, center=hl.design_point)
    assert crude.failures < 20  # p_f ~ 1e-5: crude sampling sees almost nothing
    assert weighted.importance_sampling and weighted.cov < 0.03
    assert np.isclose(weighted.failure_probability, hl.failure_probability, rtol=0.1)


def test_monte_carlo_nonlinear_limit_state_and_no_failures():
    # g = R*A - S with lognormal-free normal inputs: compare with a large crude run.
    result = r.monte_carlo_failure_probability(lambda x: x[0] * x[1] - x[2], [40.0, 5.0, 120.0], [4.0, 0.5, 30.0], samples=200_000, seed=5)
    rng = np.random.default_rng(6)
    x = rng.normal([40.0, 5.0, 120.0], [4.0, 0.5, 30.0], size=(1_000_000, 3))
    reference = np.mean(x[:, 0] * x[:, 1] < x[:, 2])
    assert abs(result.failure_probability - reference) < 4 * result.std_error
    safe = r.monte_carlo_failure_probability(lambda x: x[0] + 100, [0.0], [1.0], samples=1000, seed=0)
    assert safe.failure_probability == 0 and safe.cov == math.inf
    assert safe.confidence_interval() == (0.0, 0.0)


@pytest.mark.parametrize(
    "args,kwargs",
    [
        (("g", [1.0], [1.0]), {}),
        ((lambda x: x[0], [1.0, 2.0], [1.0]), {}),
        ((lambda x: x[0], [1.0], [1.0]), {"center": [1.0, 2.0]}),
        ((lambda x: 1.0, [1.0], [1.0]), {}),
        ((lambda x: x[0], [1.0], [1.0]), {"samples": 1}),
    ],
)
def test_monte_carlo_rejects_bad_input(args, kwargs):
    with pytest.raises(ValueError):
        r.monte_carlo_failure_probability(*args, **kwargs)


# Bayesian updating ------------------------------------------------------------


def test_conjugate_update_matches_a_numerical_posterior():
    a, b, failures, total = 2.0, 5e5, 4, 2e6
    posterior = r.bayesian_failure_rate(a, b, failures, total)
    assert (posterior.shape, posterior.rate) == (6.0, 2.5e6)
    lam = np.linspace(1e-9, 1e-5, 20001)
    density = stats.gamma.pdf(lam, a, scale=1 / b) * stats.poisson.pmf(failures, lam * total)
    density /= np.trapezoid(density, lam) if hasattr(np, "trapezoid") else np.trapz(density, lam)
    np.testing.assert_allclose(density, stats.gamma.pdf(lam, 6.0, scale=1 / 2.5e6), rtol=1e-3, atol=1e-3 * density.max())
    # The posterior mean is a weighted average of the prior mean and r/T.
    weight = b / (b + total)
    assert np.isclose(posterior.mean, weight * posterior.prior_mean + (1 - weight) * failures / total)


def test_posterior_bounds_summary_and_plot():
    posterior = r.bayesian_failure_rate(1.5, 1e6, 0, 3e6)
    assert np.isclose(posterior.upper_bound(0.9), stats.gamma.ppf(0.9, 1.5, scale=1 / 4e6))
    low, high = posterior.credible_interval(0.8)
    assert np.isclose(stats.gamma.cdf(high, 1.5, scale=1 / 4e6) - stats.gamma.cdf(low, 1.5, scale=1 / 4e6), 0.8)
    assert posterior.summary()["upper_90"] == posterior.upper_bound()
    with pytest.raises(ValueError):
        posterior.credible_interval(1.0)
    fig, ax = posterior.plot()
    assert [line.get_label() for line in ax.get_lines()] == ["prior", "posterior"]
    plt.close(fig)
    for bad in ((0.0, 1e6, 0, 1.0), (1.0, 1e6, -1, 1.0), (1.0, 1e6, 1.5, 1.0), (1.0, 1e6, 1, -1.0)):
        with pytest.raises(ValueError):
            r.bayesian_failure_rate(*bad)


def test_gamma_prior_from_bounds_round_trip():
    for median, upper, confidence in ((1e-6, 3e-6, 0.95), (2e-5, 2.2e-5, 0.9), (1e-6, 1e-4, 0.99)):
        shape, rate = r.gamma_prior_from_bounds(median, upper, confidence)
        assert np.isclose(stats.gamma.median(shape, scale=1 / rate), median, rtol=1e-9)
        assert np.isclose(stats.gamma.ppf(confidence, shape, scale=1 / rate), upper, rtol=1e-9)
    for bad in ((1e-6, 1e-6, 0.9), (1e-6, 2e-6, 0.4), (1e-300, 1e300, 0.95)):
        with pytest.raises(ValueError):
            r.gamma_prior_from_bounds(*bad)
