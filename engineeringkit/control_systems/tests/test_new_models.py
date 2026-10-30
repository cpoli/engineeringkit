"""Analytical checks for the Padé delay, Smith's predictor, Jury's test, margins over every crossover, compensators,
error constants, sensitivity, the Nichols chart, Tustin's transform, PID anti-windup and nonlinear loops."""

from math import factorial

import numpy as np
import pytest
from scipy import signal

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


# --- Lead, lag and lead-lag compensators ---------------------------------------------------------------


@pytest.mark.parametrize("phi", [10.0, 30.0, 55.0, 80.0])
def test_lead_compensator_peaks_at_its_frequency(phi):
    C = c.lead_compensator(phi, 3.0)
    w = np.logspace(-2, 3, 20001)
    phase = np.degrees(np.angle(_response(C, w)))
    assert np.isclose(w[np.argmax(phase)], 3.0, rtol=1e-3) and np.isclose(phase.max(), phi, atol=1e-6)
    alpha = (1 - np.sin(np.radians(phi))) / (1 + np.sin(np.radians(phi)))
    assert np.isclose(abs(_response(C, 3.0)), 1 / np.sqrt(alpha)) and c.dc_gain(C) == 1


def test_lag_and_lead_lag_compensators_have_unit_dc_gain():
    C = c.lag_compensator(8.0, 0.2)
    assert c.dc_gain(C) == 1 and np.isclose(C.num[0] / C.den[0], 1 / 8)
    np.testing.assert_allclose(c.zeros(C), [-0.2]) and np.testing.assert_allclose(c.poles(C), [-0.025])
    CL = c.lead_lag_compensator(40.0, 2.0, 8.0, 0.2)
    w = np.array([0.05, 2.0, 50.0])
    np.testing.assert_allclose(_response(CL, w), _response(c.lead_compensator(40.0, 2.0), w) * _response(C, w))
    for bad in (dict(beta=1.0, zero_frequency=1), dict(beta=2, zero_frequency=0)):
        with pytest.raises(ValueError):
            c.lag_compensator(**bad)
    for phi in (0, 90, -5):
        with pytest.raises(ValueError):
            c.lead_compensator(phi, 1.0)


def test_design_lead_meets_the_margin_and_keeps_the_error_constant():
    L = c.transfer_function([40], [1, 2, 0])  # Kv = 20, PM about 18 deg
    r = c.design_lead(L, 50.0)
    assert c.stability_margins(L)["phase_margin_deg"] < 20
    assert abs(r.phase_margin_deg - 50) < 2 and r.gain_cross_frequency > c.stability_margins(L)["gain_cross_frequency"]
    assert c.error_constants(r.compensated_loop)["Kv"] == pytest.approx(20)
    s = r.summary()
    assert s["high_frequency_gain"] > 1 and s["phase_margin_deg"] == r.phase_margin_deg and r.kind == "lead"


def test_design_lag_lowers_the_crossover_to_where_the_phase_suffices():
    L = c.transfer_function([5], [0.5, 1.5, 1, 0])  # Ogata's lag example: Kv = 5, unstable at this gain
    assert c.stability_margins(L)["phase_margin_deg"] < 0
    r = c.design_lag(L, 40.0)
    wc = r.gain_cross_frequency
    assert np.isclose(np.degrees(np.angle(_response(L, wc))), -180 + 45, atol=1.0)
    assert abs(r.phase_margin_deg - 40) < 1.5 and r.gain_margin > 1
    assert c.error_constants(r.compensated_loop)["Kv"] == pytest.approx(5)
    assert r.summary()["high_frequency_gain"] == pytest.approx(1 / abs(_response(L, wc)), rel=0.02)


def test_design_lead_lag_places_the_crossover():
    L = c.transfer_function([20], [1, 3, 2, 0])
    r = c.design_lead_lag(L, 45.0, 1.0)
    assert r.gain_cross_frequency == pytest.approx(1.0, rel=0.01) and abs(r.phase_margin_deg - 45) < 1.5
    assert c.dc_gain(r.compensator) == 1 and r.kind == "lead-lag"


def test_compensator_design_errors_and_plot():
    import matplotlib.pyplot as plt

    L = c.transfer_function([40], [1, 2, 0])
    with pytest.raises(ValueError, match="already"):
        c.design_lead(L, 10.0)
    with pytest.raises(ValueError, match="one lead stage"):
        c.design_lead(c.transfer_function([1000], [1, 1, 0, 0]), 60.0)
    with pytest.raises(ValueError, match="no gain crossover"):
        c.design_lead(c.transfer_function([0.1], [1, 1]), 40.0)
    with pytest.raises(ValueError, match="never reaches"):
        c.design_lag(c.transfer_function([10], [1, 0, 0]), 40.0)
    with pytest.raises(ValueError, match="below one"):
        c.design_lag(c.transfer_function([0.01], [1, 3, 2, 0]), 40.0)
    with pytest.raises(ValueError, match="lead"):
        c.design_lead_lag(L, 45.0, 0.01)
    with pytest.raises(ValueError, match="design_lead"):
        c.design_lead_lag(c.transfer_function([0.5], [1, 3, 2, 0]), 45.0, 1.0)
    for bad in (dict(zero_ratio=1.0), dict(phase_margin_deg=95.0), dict(extra_deg=-1.0)):
        kwargs = dict(phase_margin_deg=40.0) | bad
        with pytest.raises(ValueError):
            c.design_lag(L, **kwargs)
    with pytest.raises(ValueError):
        c.design_lead_lag(L, 45.0, 1.0, zero_ratio=0.5)
    fig, axes = c.design_lead(L, 50.0).plot(np.logspace(-1, 2, 50))
    assert len(axes[0].lines) == 2 and axes[0].get_legend() is not None
    plt.close(fig)


# --- System type and error constants -------------------------------------------------------------------


def test_error_constants_by_system_type():
    e0 = c.error_constants(c.transfer_function([6], [1, 5, 6]))  # Kp = 1
    assert (e0["type"], e0["Kp"], e0["Kv"], e0["Ka"]) == (0, 1.0, 0.0, 0.0)
    assert e0["step_error"] == 0.5 and e0["ramp_error"] == np.inf and e0["parabola_error"] == np.inf
    e1 = c.error_constants(c.transfer_function([10, 10], [1, 5, 0]))
    assert (e1["type"], e1["Kp"], e1["Kv"], e1["step_error"], e1["ramp_error"]) == (1, np.inf, 2.0, 0.0, 0.5)
    e2 = c.error_constants(c.transfer_function([4, 2], [1, 2, 0, 0]))
    assert (e2["type"], e2["Ka"], e2["ramp_error"], e2["parabola_error"]) == (2, 1.0, 0.0, 1.0)
    assert c.system_type(c.transfer_function([1, 0], [1, 1])) == -1
    assert c.error_constants(c.transfer_function([1, 0], [1, 1]))["Kp"] == 0.0
    with pytest.raises(ValueError, match="not stable"):
        c.error_constants(c.transfer_function([10], [1, 3, 2, 0]))


def test_ramp_error_matches_simulation():
    L = c.transfer_function([10, 10], [1, 5, 0])
    t = np.linspace(0, 20, 4001)
    _, y, _ = signal.lsim(c.feedback(L), t, t)
    assert np.isclose(t[-1] - y[-1], c.error_constants(L)["ramp_error"], atol=1e-6)


# --- Sensitivity ---------------------------------------------------------------------------------------


def test_sensitivity_and_complementary_sensitivity_sum_to_one():
    L = c.transfer_function([2, 1], [1, 3, 2, 0])
    w = np.logspace(-2, 2, 40)
    S, T = _response(c.sensitivity(L), w), _response(c.complementary_sensitivity(L), w)
    np.testing.assert_allclose(S + T, 1, atol=1e-12)
    np.testing.assert_allclose(T, _response(c.feedback(L), w), rtol=1e-12)


def test_sensitivity_peak_bounds_the_margins():
    L = c.transfer_function([4], [1, 3, 2, 0])
    p = c.sensitivity_peak(L)
    m = c.stability_margins(L)
    w = np.logspace(-4, 4, 20000)
    assert np.isclose(p["Ms"], 1 / np.min(np.abs(1 + c.nyquist(L, w))))
    assert m["gain_margin"] >= p["gain_margin_bound"] and m["phase_margin_deg"] >= p["phase_margin_bound_deg"]
    assert p["Mt"] > 1 and 0 < p["Ms_frequency"] < 10
    far = c.sensitivity_peak(c.transfer_function([10, 0], [1, 1]), np.logspace(-2, 2, 50))  # |1 + L| > 1 everywhere
    assert far["Ms"] < 1 and far["gain_margin_bound"] == np.inf and far["phase_margin_bound_deg"] > 60


# --- Nichols chart -------------------------------------------------------------------------------------


def test_nichols_closed_loop_matches_direct_evaluation():
    L = 10 ** (np.array([-6.0, 0.0, 3.0]) / 20) * np.exp(1j * np.radians([-200.0, -120.0, -150.0]))
    m, p = c.nichols_closed_loop(20 * np.log10(np.abs(L)), np.degrees(np.angle(L)))
    T = L / (1 + L)
    np.testing.assert_allclose(m, 20 * np.log10(np.abs(T))) and np.testing.assert_allclose(p, np.degrees(np.angle(T)))
    assert c.nichols_closed_loop(0.0, -180.0)[0] > 300  # the critical point
    with pytest.raises(ValueError):
        c.nichols_closed_loop(np.nan, 0.0)


def test_nichols_resonant_peak_of_a_second_order_loop():
    # 1/(s(s+1)) closes to 1/(s^2 + s + 1): zeta = 1/2, Mr = 1/(2 zeta sqrt(1 - zeta^2)), wr = wn sqrt(1 - 2 zeta^2).
    n = c.nichols(c.transfer_function([1], [1, 1, 0]), np.linspace(0.5, 1.0, 50001))
    assert np.isclose(n["resonant_peak_db"], 20 * np.log10(2 / np.sqrt(3)), atol=1e-9)
    assert np.isclose(n["resonant_frequency"], np.sqrt(0.5), atol=1e-4)
    assert np.all(n["phase_deg"] < -90) and n["closed_phase_deg"].shape == n["w"].shape


# --- Tustin's bilinear transform -----------------------------------------------------------------------


def test_tustin_maps_poles_and_matches_at_the_prewarp_frequency():
    G = c.transfer_function([2, 3], [1, 3, 2])
    dt = 0.2
    Gd = c.discretize_tustin(G, dt)
    np.testing.assert_allclose(np.sort(np.roots(Gd.den)), np.sort((2 / dt + c.poles(G)) / (2 / dt - c.poles(G))))
    assert Gd.dt == dt and np.isclose(np.polyval(Gd.num, 1) / np.polyval(Gd.den, 1), c.dc_gain(G))
    w0 = 6.0
    Gp = c.discretize_tustin(G, dt, prewarp_frequency=w0)
    z = np.exp(1j * w0 * dt)
    assert np.isclose(np.polyval(Gp.num, z) / np.polyval(Gp.den, z), _response(G, w0))
    zd = np.exp(1j * w0 * dt)
    assert not np.isclose(np.polyval(Gd.num, zd) / np.polyval(Gd.den, zd), _response(G, w0), rtol=1e-3)


def test_tustin_input_checks():
    with pytest.raises(ValueError, match="proper"):
        c.discretize_tustin(c.transfer_function([1, 0, 0], [1, 1]), 0.1)
    with pytest.raises(ValueError, match="Nyquist"):
        c.discretize_tustin(c.transfer_function([1], [1, 1]), 0.1, prewarp_frequency=40.0)
    with pytest.raises(ValueError):
        c.discretize_tustin(c.transfer_function([1], [1, 1]), 0.0)


# --- PID anti-windup -----------------------------------------------------------------------------------


def test_saturated_pid_without_active_limits_matches_the_linear_loop():
    plant = c.transfer_function([1], [1, 3, 2])
    pid = c.PID(3.0, 2.0)
    t = np.linspace(0, 10, 2001)
    linear = pid.simulate(plant, t)
    for mode in (None, "clamping", "back_calculation"):
        r = pid.simulate_saturated(plant, t, -100, 100, anti_windup=mode)
        np.testing.assert_allclose(r.y, linear.y, atol=1e-8)


def test_anti_windup_cuts_the_overshoot_of_a_saturated_integrator():
    plant = c.transfer_function([1], [1, 0])  # an integrating plant winds up badly
    pid = c.PID(2.0, 1.0)
    t = np.linspace(0, 30, 6001)
    runs = {mode: pid.simulate_saturated(plant, t, -0.4, 0.4, anti_windup=mode) for mode in (None, "clamping", "back_calculation")}
    for r in runs.values():
        assert np.all(np.abs(r.u) <= 0.4 + 1e-12) and abs(r.y[-1] - 1) < 2e-3
        assert r.summary()["input_max"] == pytest.approx(0.4)
    assert runs[None].summary()["peak"] > 1.5 * runs["clamping"].summary()["peak"] - 0.5
    assert runs[None].summary()["peak"] - 1 > 3 * (runs["back_calculation"].summary()["peak"] - 1)
    assert runs[None].summary()["peak"] - 1 > 3 * (runs["clamping"].summary()["peak"] - 1)


def test_saturated_pid_derivative_on_measurement_and_input_checks():
    plant = c.transfer_function([1], [1, 3, 2])
    t = np.linspace(0, 8, 1601)
    r = c.PID(4.0, 2.0, 0.5).simulate_saturated(plant, t, -50, 50, "back_calculation")
    assert abs(r.y[-1] - 1) < 1e-2 and r.u[0] == pytest.approx(4.0)  # no derivative kick
    no_integral = c.PID(4.0).simulate_saturated(plant, t, -1, 1, "back_calculation")
    assert no_integral.y[-1] < 1
    explicit = c.PID(0.0, 1.0).simulate_saturated(plant, t, -1, 1, "back_calculation", tracking_time=1.0)
    assert np.all(np.isfinite(explicit.y))
    with pytest.raises(ValueError, match="tracking_time"):
        c.PID(0.0, 1.0).simulate_saturated(plant, t, -1, 1, "back_calculation")
    with pytest.raises(ValueError, match="anti_windup"):
        c.PID(1.0, 1.0).simulate_saturated(plant, t, -1, 1, "reset")
    with pytest.raises(ValueError, match="below"):
        c.PID(1.0, 1.0).simulate_saturated(plant, t, 1, 1)
    with pytest.raises(ValueError, match="strictly proper"):
        c.PID(1.0, 1.0).simulate_saturated(c.transfer_function([1, 0], [1, 1]), t, -1, 1)
    with pytest.raises(ValueError, match="relative degree"):
        c.PID(1.0, 1.0, 1.0).simulate_saturated(c.transfer_function([1], [1, 1]), t, -1, 1)


# --- Describing functions, limit cycles and the phase plane --------------------------------------------


def test_describing_functions_against_numerical_harmonics():
    A = np.array([0.3, 1.0, 2.5, 7.0])
    relay = c.describing_function(lambda x: 2.0 * np.sign(x), A, samples=2**16)
    np.testing.assert_allclose(relay, c.relay_describing_function(A, 2.0), rtol=1e-4)
    sat = c.describing_function(lambda x: 3.0 * np.clip(x, -1.0, 1.0), A)
    np.testing.assert_allclose(sat, c.saturation_describing_function(A, 1.0, 3.0), rtol=1e-6)
    dead = c.describing_function(lambda x: 2.0 * np.sign(x) * np.maximum(np.abs(x) - 0.5, 0), A)
    np.testing.assert_allclose(dead, c.dead_zone_describing_function(A, 0.5, 2.0), atol=1e-6)
    np.testing.assert_allclose(c.describing_function(lambda x: x**3, A), 0.75 * A**2, rtol=1e-12)
    assert isinstance(c.describing_function(np.sin, 1.0), complex)


def test_relay_with_hysteresis_lags_by_arcsin():
    A, M, h = 2.0, 1.5, 0.6
    theta = 2 * np.pi * (np.arange(2**16) + 0.5) / 2**16
    # Up from the rising crossing of +h to the falling crossing of -h.
    y = np.where((theta > np.arcsin(h / A)) & (theta < np.pi + np.arcsin(h / A)), M, -M)
    numeric = 2 * np.mean(y * (np.sin(theta) + 1j * np.cos(theta))) / A
    assert np.isclose(c.relay_describing_function(A, M, h), numeric, rtol=1e-4)
    assert np.angle(c.relay_describing_function(A, M, h)) == pytest.approx(-np.arcsin(h / A))
    with pytest.raises(ValueError, match="hysteresis"):
        c.relay_describing_function(0.5, 1.0, 0.6)


def test_limit_cycle_of_a_relay_and_of_a_saturation():
    (relay,) = c.limit_cycle(c.transfer_function([1], [1, 3, 3, 1]), lambda A: c.relay_describing_function(A, 1.0), (0.01, 10))
    assert relay["amplitude"] == pytest.approx(1 / (2 * np.pi), rel=1e-4) and relay["stable"]
    assert relay["frequency"] == pytest.approx(np.sqrt(3), rel=1e-3) and relay["period"] == pytest.approx(2 * np.pi / np.sqrt(3), rel=1e-3)
    L = c.transfer_function([10], [1, 3, 2, 0])  # gain margin 0.6
    (sat,) = c.limit_cycle(L, lambda A: c.saturation_describing_function(A, 1.0), (0.1, 100))
    assert c.saturation_describing_function(sat["amplitude"], 1.0) == pytest.approx(0.6) and sat["stable"]
    # A dead zone makes the small-amplitude gain zero: the cycle it predicts is unstable.
    (dz,) = c.limit_cycle(L, lambda A: c.dead_zone_describing_function(A, 1.0), (0.1, 100))
    assert not dz["stable"] and c.dead_zone_describing_function(dz["amplitude"], 1.0) == pytest.approx(0.6)
    assert c.limit_cycle(c.transfer_function([1], [1, 1]), lambda A: c.relay_describing_function(A, 1.0), (0.1, 1)) == []


def test_limit_cycle_input_checks():
    L = c.transfer_function([1], [1, 3, 3, 1])
    with pytest.raises(ValueError, match="real"):
        c.limit_cycle(L, lambda A: c.relay_describing_function(A, 1.0, 0.01), (0.1, 1))
    with pytest.raises(ValueError, match="increasing"):
        c.limit_cycle(L, lambda A: A, (1.0, 0.1))
    with pytest.raises(ValueError):
        c.describing_function(np.sign, 1.0, samples=4)


def test_phase_portrait_of_a_relay_servo_conserves_its_energy():
    M = 2.0
    r = c.phase_portrait(lambda x, v: -M * np.sign(x), 1.0, 0.0, np.linspace(0, 10, 20001))
    energy = M * np.abs(r.x) + r.v**2 / 2
    np.testing.assert_allclose(energy, M, rtol=2e-3)
    s = r.summary()
    assert s["max_abs_x"] == pytest.approx(1.0, rel=1e-3) and set(s) == {"final_x", "final_v", "max_abs_x"}
    with pytest.raises(ValueError):
        c.phase_portrait(lambda x, v: -x, np.nan, 0.0, [0, 1])


def test_phase_portrait_plot_draws_field_and_trajectory():
    import matplotlib.pyplot as plt

    r = c.phase_portrait(lambda x, v: -x - 0.5 * v, 1.0, 0.0, np.linspace(0, 10, 501))
    fig, ax = r.plot()
    assert len(ax.collections) == 1 and ax.get_xlabel() == "x"
    _, ax2 = r.plot(field=False)
    assert not ax2.collections
    plt.close("all")
