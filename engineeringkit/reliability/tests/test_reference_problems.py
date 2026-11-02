"""Worked examples from the NIST/SEMATECH *e-Handbook of Statistical Methods*, ch. 8.

The examples are from chapter 8, "Assessing Product Reliability"
(https://www.itl.nist.gov/div898/handbook/apr/apr.htm), cited by section.
Section D names Ebeling and O'Connor and Kleyner for this domain, but no
legitimate copy of either was at hand, so these public, citable examples
stand in. Each test reproduces the answer the handbook prints, to the
precision it prints it. Where the handbook rounds a constant or an
intermediate result, the test says so and the reliability guide records the
difference. The handbook works in hours; the tests convert to seconds.
"""

import numpy as np
from scipy.constants import eV, k
from scipy.special import gamma

from engineeringkit import reliability as r

HOUR = 3600.0
# Section 8.4.1.3 and 8.2.2.1: ten failures among 20 units, the other ten
# still running when the test stopped at 500 h.
FAILURES = [54, 187, 216, 240, 244, 335, 361, 373, 375, 386]
TIMES = np.array(FAILURES + [500] * 10, dtype=float) * HOUR
FAILED = [1] * 10 + [0] * 10


def test_nist_8_1_5_1_arrhenius_acceleration_factor():
    # 25 C to 125 C: 133 for dH = 0.5 eV and 17,597 for dH = 1.0 eV. The
    # handbook takes 1/k = 11605 K/eV and 273.16 K for 0 C (as in sec.
    # 8.4.3); CODATA k and 273.15 K give 17,600 for the second.
    assert round(float(r.arrhenius_acceleration_factor(0.5 * eV, 298.15, 398.15))) == 133
    assert round(float(r.arrhenius_acceleration_factor(1.0 * eV, 298.15, 398.15))) == 17600
    handbook_energy = 1.0 * eV * (k / eV) * 11605
    assert round(float(r.arrhenius_acceleration_factor(handbook_energy, 298.16, 398.16))) == 17597


def test_nist_8_4_5_1_mtbf_chi_square_bounds():
    # 800 h, 2 failures, Type I (time) censoring: MTBF 400 h and the 90 %
    # two-sided interval (127, 2251) h. With no failures the lower bound is
    # 800/(-ln 0.05) = 267 h.
    lower, upper = r.mtbf_confidence_bounds(800 * HOUR, 2, 0.9, time_terminated=True)
    assert 800 / 2 == 400
    assert (round(lower / HOUR), round(upper / HOUR)) == (127, 2251)
    lower, upper = r.mtbf_confidence_bounds(800 * HOUR, 0, 0.9, time_terminated=True)
    assert round(lower / HOUR) == 267 and np.isinf(upper)


def test_nist_8_4_1_3_weibull_maximum_likelihood():
    # eta = 606.5280 h, beta = 1.7208, log-likelihood -75.135, AIC 154.27,
    # mean 540.737 h and standard deviation 323.806 h. The likelihood is
    # maximised at eta = 606.5276 h; the handbook's optimiser stopped
    # 4e-4 h away, where the log-likelihood is the same to 1e-10.
    fit = r.weibull_mle(TIMES, FAILED)
    assert round(fit.beta, 4) == 1.7208
    assert round(fit.eta / HOUR, 2) == 606.53
    assert abs(fit.eta / HOUR - 606.5280) < 1e-3

    t, f = TIMES / HOUR, np.array(FAILURES, dtype=float)
    eta = fit.eta / HOUR
    loglik = np.sum(np.log(fit.beta / eta) + (fit.beta - 1) * np.log(f / eta)) - np.sum((t / eta) ** fit.beta)
    assert round(float(loglik), 3) == -75.135
    assert round(float(2 * 2 - 2 * loglik), 2) == 154.27

    mean = r.weibull_mttf(eta, fit.beta)
    variance = eta**2 * (gamma(1 + 2 / fit.beta) - gamma(1 + 1 / fit.beta) ** 2)
    assert round(float(mean), 3) == 540.737
    assert round(float(np.sqrt(variance)), 3) == 323.806


def test_nist_8_2_2_1_weibull_probability_plot():
    # Plotted points ln(1/(1 - F)) at Benard's median ranks of 20 units, then
    # the regression of log10 of them on log10 t: slope 1.46. The handbook
    # takes alpha = 10^(4.114/1.46) = 657 h from the rounded slope and an
    # intercept fitted to the 3-decimal points; unrounded, the fit is 667 h.
    fit = r.weibull_rank_regression(TIMES, FAILED, ranks="benard")
    printed = [0.035, 0.087, 0.142, 0.200, 0.262, 0.328, 0.398, 0.474, 0.556, 0.645]
    np.testing.assert_array_equal(np.round(-np.log1p(-fit.median_ranks), 3), printed)
    assert round(fit.beta, 2) == 1.46
    assert round(10 ** (4.114 / 1.46)) == 657
    assert round(fit.eta / HOUR) == 667


def test_nist_8_2_5_gamma_prior_from_50_and_95_percent_mtbf():
    # MTBF_50 = 600 h and MTBF_05 = 250 h, so the rate is below 1/600 with
    # probability 0.5 and below 1/250 with 0.95: a = 2.863, b = 1522.46 h.
    # The quantile ratio fixes a = 2.86306 and then b = 1522.51 h; the
    # handbook's solver stopped slightly short (b to 3e-5 of it).
    shape, rate = r.gamma_prior_from_bounds(1 / (600 * HOUR), 1 / (250 * HOUR), 0.95)
    assert round(shape, 3) == 2.863
    assert abs(rate / HOUR - 1522.46) / 1522.46 < 5e-5


def test_nist_8_4_6_gamma_posterior_mtbf_bounds():
    # Prior a = 2, b = 1400 h; 1909 h on test with 2 failures gives a' = 4,
    # b' = 3309 h. MTBF bounds 495 h (90 %), 600 h (80 %), 901 h (50 %) and
    # 1897 h (10 %); (495, 1897) is the 90 % credibility interval and
    # b'/a' = 827 h.
    posterior = r.bayesian_failure_rate(2.0, 1400 * HOUR, 2, 1909 * HOUR)
    assert (posterior.shape, posterior.rate / HOUR) == (4.0, 3309.0)
    assert round(1 / posterior.upper_bound(0.9) / HOUR) == 495
    assert round(1 / posterior.upper_bound(0.8) / HOUR) == 600
    assert round(1 / posterior.upper_bound(0.5) / HOUR) == 901
    assert round(1 / posterior.upper_bound(0.1) / HOUR) == 1897
    high_rate, low_rate = posterior.credible_interval(0.8)
    assert (round(1 / low_rate / HOUR), round(1 / high_rate / HOUR)) == (495, 1897)
    assert round(1 / posterior.mean / HOUR) == 827


def test_nist_8_2_1_5_kaplan_meier_product_limit():
    # 20 units: failures at 10, 32, 56, 98, 122 and 181 h, removals at 50,
    # 100, 125 and 150 h and ten survivors at 200 h. The handbook prints the
    # products, e.g. R(181) = 19/20 x 18/19 x 16/17 x 15/16 x 13/14 x 10/11.
    times = [10, 32, 56, 98, 122, 181, 50, 100, 125, 150] + [200] * 10
    result = r.kaplan_meier(np.array(times, dtype=float) * HOUR, [1] * 6 + [0] * 14)
    np.testing.assert_array_equal(result.at_risk, [20, 19, 17, 16, 14, 11])
    factors = [19 / 20, 18 / 19, 16 / 17, 15 / 16, 13 / 14, 10 / 11]
    np.testing.assert_allclose(result.survival, np.cumprod(factors), rtol=0, atol=1e-12)


def test_nist_8_4_3_lognormal_projection_to_use_conditions():
    # T50 = 507,383 h at use and sigma = 0.74: F(100,000 h) = 0.014.
    assert round(float(r.lognormal_cdf(1e5 * HOUR, np.log(507383 * HOUR), 0.74)), 3) == 0.014
