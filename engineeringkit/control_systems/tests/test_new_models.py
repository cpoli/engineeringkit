"""Analytical checks for the Padé delay, Smith's predictor, Jury's test and margins over every crossover."""

from math import factorial

import numpy as np
import pytest

from engineeringkit import control_systems as c


def _response(sys, w):
    return np.polyval(sys.num, 1j * w) / np.polyval(sys.den, 1j * w)


@pytest.mark.parametrize("order", [1, 2, 3, 6])
def test_pade_delay_matches_taylor_series_and_is_all_pass(order):
    theta = 0.7
    D = c.pade_delay(theta, order)
    # D(s) agrees with the Taylor series of exp(-s theta) through s^(2n): den * series - num = O(s^(2n+1)).
    taylor = np.array([(-theta) ** k / factorial(k) for k in range(2 * order + 1)])
    residual = np.convolve(D.den[::-1], taylor)[: 2 * order + 1] - np.pad(D.num[::-1], (0, order))
    np.testing.assert_allclose(residual / np.abs(D.den).max(), 0, atol=1e-13)
    w = np.logspace(-2, 1, 50)
    np.testing.assert_allclose(np.abs(_response(D, w)), 1, rtol=1e-12)
    assert np.sum(c.zeros(D).real > 0) == order


def test_pade_delay_phase_converges_with_order():
    w = np.array([0.5, 2.0, 4.0])
    errors = [np.abs(np.unwrap(np.angle(_response(c.pade_delay(1.0, n), w))) + w).max() for n in (2, 4, 8)]
    assert errors[0] > errors[1] > errors[2] and errors[2] < 1e-6


@pytest.mark.parametrize("bad", [dict(delay=0), dict(delay=-1), dict(delay=1, order=0), dict(delay=1, order=1.5)])
def test_pade_delay_rejects_invalid_input(bad):
    with pytest.raises(ValueError):
        c.pade_delay(**bad)


def test_smith_predictor_removes_delay_from_characteristic_equation():
    G = c.transfer_function([1], [1, 3, 2])
    C = c.transfer_function([3, 2], [1, 0])
    D = c.pade_delay(1.5, 4)
    Csp = c.smith_predictor(C, G, 1.5, order=4)
    loop = c.transfer_function(np.polymul(np.polymul(Csp.num, G.num), D.num), np.polymul(np.polymul(Csp.den, G.den), D.den))
    closed = c.poles(c.feedback(loop))
    for p in c.poles(c.feedback(c.transfer_function(np.polymul(C.num, G.num), np.polymul(C.den, G.den)))):
        assert np.min(np.abs(closed - p)) < 1e-6
    # Reference to output is C G D / (1 + C G) at every frequency.
    w = np.logspace(-2, 1, 30)
    cg = _response(C, w) * _response(G, w)
    np.testing.assert_allclose(_response(c.feedback(loop), w), cg * _response(D, w) / (1 + cg), rtol=1e-6)


def test_jury_stable_agrees_with_roots():
    rng = np.random.default_rng(7)
    for _ in range(2000):
        coeffs = rng.normal(size=rng.integers(2, 8)) * rng.choice([0.3, 1, 3])
        assert c.jury_stable(coeffs) == bool(np.all(np.abs(np.roots(coeffs)) < 1))


def test_jury_stable_boundary_and_input_checks():
    assert c.jury_stable([2, 1]) and not c.jury_stable([1, 1])  # root at -1/2, root at -1
    assert not c.jury_stable([1, 0, 1])  # roots on the unit circle
    assert c.jury_stable([0, 1, 0.5])  # leading zeros are stripped
    for bad in ([1], [1, np.nan], [0, 0]):
        with pytest.raises(ValueError):
            c.jury_stable(bad)


def test_margins_of_a_conditionally_stable_loop_use_the_nearest_crossover():
    # 5 (s+1)^2 / (s^3 (0.1 s + 1)^2): phase crosses -180 deg twice.
    num = 5 * np.polymul([1, 1], [1, 1])
    den = np.polymul([1, 0, 0, 0], np.polymul([0.1, 1], [0.1, 1]))
    m = c.stability_margins(c.transfer_function(num, den))
    assert m["gain_margins"].size == 2 and m["gain_margins"][0] < 1 < m["gain_margins"][1]
    assert m["gain_margin"] == m["gain_margins"][1] and m["phase_cross_frequency"] == m["phase_cross_frequencies"][1]
    for wp, gm in zip(m["phase_cross_frequencies"], m["gain_margins"]):
        edge = c.poles(c.feedback(c.transfer_function(gm * num, den)))
        assert np.min(np.abs(edge - 1j * wp)) < 1e-3


def test_phase_margin_is_the_crossover_of_smallest_magnitude():
    # LQR loop of an open-loop unstable plant crosses the unit circle on both sides of -1.
    A, B = np.array([[1.0, 1.0], [0.0, 1.0]]), np.array([[0.0], [1.0]])
    ((k1, k2),), _, _ = c.lqr(A, B, [[1, 1], [1, 1]], [[1.0]])
    loop = c.transfer_function([k2, k1 - k2], [1, -2, 1])  # K (sI - A)^-1 B
    m = c.stability_margins(loop)
    assert m["phase_margins_deg"].size == 2
    assert np.all(np.abs(m["phase_margins_deg"]) >= 60)
    assert abs(m["phase_margin_deg"]) == np.abs(m["phase_margins_deg"]).min()
    assert np.all(m["gain_margins"] <= 0.5)
