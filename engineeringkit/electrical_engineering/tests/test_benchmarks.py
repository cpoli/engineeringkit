import numpy as np
import pytest
from scipy.integrate import trapezoid

from engineeringkit import electrical_engineering as e


def test_rlc_step_is_current_and_charge_conserves():
    R, L, C = 2.0, 1.0, 0.25
    result = e.SeriesRLC(R, L, C).step_response(np.linspace(0, 20, 20001))
    wd = np.sqrt(1 / (L * C) - (R / (2 * L)) ** 2)
    expected = np.exp(-R * result.t / (2 * L)) * np.sin(wd * result.t) / (L * wd)
    np.testing.assert_allclose(result.response, expected, atol=2e-12)
    # A 1 V step deposits Q=C*V, not a nonzero steady current.
    assert np.isclose(trapezoid(result.response, result.t), C, rtol=1e-6)
    assert abs(result.response[-1]) < 1e-8


def test_rc_charge_current_energy_balance():
    t = np.linspace(0, 20, 100001)
    R, C, V = 2.0, 0.5, 3.0
    voltage = e.rc_charging_voltage(V, R, C, t)
    current = (V - voltage) / R
    source = trapezoid(V * current, t)
    dissipated = trapezoid(R * current**2, t)
    stored = 0.5 * C * voltage[-1] ** 2
    assert np.isclose(source, dissipated + stored, rtol=1e-7)


def test_rlc_resonance_and_lossless_line_identity():
    omega = 2 * np.pi * e.resonance_frequency(0.1, 1e-4)
    assert np.isclose(e.series_rlc_impedance(10, 0.1, 1e-4, omega), 10)
    L, C = 2e-7, 8e-11
    assert np.isclose(e.transmission_line_z0(L, C) * e.transmission_line_velocity(L, C), 1 / C)


@pytest.mark.parametrize("R,L,C", [(-1, 1, 1), (1, 0, 1), (1, 1, 0), (np.nan, 1, 1)])
def test_rlc_invalid_components(R, L, C):
    with pytest.raises(ValueError):
        e.SeriesRLC(R, L, C)


@pytest.mark.parametrize("R", [0.5, 2.0, 5.0])
def test_rlc_step_voltage_against_scipy_and_series_rlc_current(R):
    from scipy import signal

    L, C = 1.0, 1.0
    t = np.linspace(0, 15, 3001)
    _, v_ref = signal.step(signal.TransferFunction([1], [L * C, R * C, 1]), T=t)
    v = e.rlc_step_voltage(1.0, R, L, C, t)
    np.testing.assert_allclose(v, v_ref, atol=1e-6)
    # Loop current i = C dv/dt is what SeriesRLC.step_response returns.
    current = e.SeriesRLC(R, L, C).step_response(t).response
    np.testing.assert_allclose((C * np.gradient(v, t))[1:-1], current[1:-1], atol=5e-5)


def test_diode_operating_point_against_lambert_w():
    from scipy.special import lambertw

    Vs, R, Is, n = 3.0, 220.0, 2e-9, 1.8
    nvt = n * e.thermal_voltage(300.0)
    # Closed form of the load line crossing the Shockley curve.
    i_exact = nvt / R * lambertw(Is * R / nvt * np.exp((Vs + Is * R) / nvt)).real - Is
    v, i = e.diode_operating_point(Vs, R, Is, ideality=n)
    assert i == pytest.approx(i_exact, rel=1e-10)


def test_newton_raphson_matches_gauss_seidel_on_a_meshed_network():
    lines = [
        (0, 1, 0.02 + 0.06j, 0.06j),
        (0, 2, 0.08 + 0.24j, 0.05j),
        (1, 2, 0.06 + 0.18j, 0.04j),
        (1, 3, 0.06 + 0.18j, 0.04j),
        (2, 3, 0.01 + 0.03j, 0.02j),
    ]
    Y = e.admittance_matrix(4, lines)
    p, q = [0.2, 0.45, 0.4], [0.1, 0.15, 0.05]
    gs = e.gauss_seidel_load_flow(Y, p, q, tol=1e-12, max_iter=5000)
    nr = e.newton_raphson_load_flow(Y, p, q, tol=1e-12)
    assert gs.converged and nr.converged and nr.iterations < 8 < gs.iterations
    np.testing.assert_allclose(nr.voltages, gs.voltages, atol=1e-10)
    # Power balance: the bus injections sum to the line losses, which are positive.
    assert nr.losses.real > 0 and nr.power[1:].real == pytest.approx(-np.array(p))


def test_swing_curve_conserves_energy_and_brackets_critical_clearing_time():
    H, Pm, P1, P3, f = 4.0, 0.8, 1.6, 1.6, 60.0
    ws = 2 * np.pi * f
    tcr = e.critical_clearing_time(H, Pm, P1, frequency=f, p_max_post=P3)
    assert e.swing_response(H, Pm, P1, 0.0, P3, 0.97 * tcr, 2.0, frequency=f).stable
    assert not e.swing_response(H, Pm, P1, 0.0, P3, 1.03 * tcr, 2.0, frequency=f).stable
    # After clearing, kinetic plus potential energy is constant without damping.
    r = e.swing_response(H, Pm, P1, 0.0, P3, 0.8 * tcr, 2.0, frequency=f)
    post = r.t > r.clearing_time
    energy = H / ws * r.speed_deviation[post] ** 2 - Pm * r.delta[post] - P3 * np.cos(r.delta[post])
    assert np.ptp(energy) < 1e-6
    # During a bolted fault the angle is a parabola.
    fault = r.t <= r.clearing_time
    np.testing.assert_allclose(r.delta[fault], r.delta[0] + ws * Pm * r.t[fault] ** 2 / (4 * H), atol=1e-8)
    # Damping lowers the first swing.
    damped = e.swing_response(H, Pm, P1, 0.0, P3, 0.8 * tcr, 2.0, frequency=f, damping=0.05)
    assert damped.delta.max() < r.delta.max()
    for kwargs in [{"clearing_time": 0.0}, {"clearing_time": 2.0}, {"n_points": 2}, {"damping": -1.0}]:
        args = dict(inertia_constant=H, p_mech=Pm, p_max_pre=P1, p_max_fault=0.0, p_max_post=P3, clearing_time=0.1, t_end=2.0) | kwargs
        with pytest.raises(ValueError):
            e.swing_response(**args)
