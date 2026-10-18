"""Analytical checks for nodal analysis, transmission lines, symmetrical components and load flow."""

import numpy as np
import pytest

from engineeringkit import electrical_engineering as e


def test_nodal_analysis_voltage_divider_and_bridge():
    r = e.nodal_analysis(2, [(1, 2, 1e3), (2, 0, 3e3)], voltage_sources=[(1, 0, 8.0)])
    np.testing.assert_allclose(r.voltages, [0, 8, 6])
    np.testing.assert_allclose(r.resistor_currents, [2e-3, 2e-3])
    np.testing.assert_allclose(r.source_currents, [2e-3])
    assert r.kcl_residual < 1e-15
    assert r.summary() == {"nodes": 3, "max_voltage": 8.0, "kcl_residual": r.kcl_residual}

    # Unbalanced Wheatstone bridge: galvanometer current from the closed form.
    R1, R2, R3, R4, Rg, V = 100.0, 200.0, 150.0, 250.0, 50.0, 10.0
    bridge = [(1, 2, R1), (2, 0, R2), (1, 3, R3), (3, 0, R4), (2, 3, Rg)]
    r = e.nodal_analysis(3, bridge, voltage_sources=[(1, 0, V)])
    vth = V * (R2 / (R1 + R2) - R4 / (R3 + R4))
    rth = R1 * R2 / (R1 + R2) + R3 * R4 / (R3 + R4)
    assert r.resistor_currents[4] == pytest.approx(vth / (rth + Rg))


def test_nodal_analysis_current_sources_and_superposition():
    resistors = [(1, 0, 2.0), (1, 2, 4.0), (2, 0, 4.0)]
    # Each source acting alone, the other zeroed (a 0 V source is a short).
    a = e.nodal_analysis(2, resistors, current_sources=[(0, 1, 3.0)], voltage_sources=[(2, 0, 0.0)])
    b = e.nodal_analysis(2, resistors, voltage_sources=[(2, 0, 2.0)])
    both = e.nodal_analysis(2, resistors, current_sources=[(0, 1, 3.0)], voltage_sources=[(2, 0, 2.0)])
    np.testing.assert_allclose(both.voltages, a.voltages + b.voltages)
    np.testing.assert_allclose(both.source_currents, a.source_currents + b.source_currents)
    # 3 A into 2 || 4 ohm with node 2 shorted.
    assert a.voltages[1] == pytest.approx(3 * 4 / 3)
    # A current source between two non-ground nodes.
    r = e.nodal_analysis(2, [(1, 0, 1.0), (2, 0, 1.0)], current_sources=[(1, 2, 1.0)])
    np.testing.assert_allclose(r.voltages, [0, -1, 1])
    assert r.kcl_residual < 1e-15


@pytest.mark.parametrize(
    "args",
    [
        (0, [(1, 0, 1.0)]),
        (1, []),
        (1, [(1, 0)]),
        (1, [(2, 0, 1.0)]),
        (1, [(1, 1, 1.0)]),
        (1, [(1, 0, 0.0)]),
        (1, [(1, 0, np.nan)]),
        (2, [(1, 0, 1.0)]),
        (1, [(1, 0, 1.0)], (), [(1, 0, 1.0), (1, 0, 2.0)]),
    ],
)
def test_nodal_analysis_rejects_bad_networks(args):
    with pytest.raises(ValueError):
        e.nodal_analysis(*args)


def test_lossless_and_distortionless_lines():
    L, C = 250e-9, 100e-12
    w = np.array([1e3, 1e6, 1e9])
    gamma = e.propagation_constant(0, L, 0, C, w)
    np.testing.assert_allclose(gamma.real, 0, atol=1e-18)
    np.testing.assert_allclose(w / gamma.imag, e.transmission_line_velocity(L, C))

    R, G = 0.1, 0.1 * C / L
    assert e.distortionless_condition(R, L, G, C)
    assert not e.distortionless_condition(R, L, G / 10, C)
    gamma = e.propagation_constant(R, L, G, C, w)
    np.testing.assert_allclose(gamma.real, np.sqrt(R * G))
    np.testing.assert_allclose(w / gamma.imag, 1 / np.sqrt(L * C))


def test_reflection_vswr_and_input_impedance():
    z0 = 50.0
    assert e.reflection_coefficient(z0, z0) == 0
    assert e.reflection_coefficient(0, z0) == -1
    assert e.standing_wave_ratio(e.reflection_coefficient(150.0, z0)) == pytest.approx(3)
    # Quarter-wave transformer: Z_in = Z0^2 / Z_L, and a half-wave line repeats the load.
    zl = 30 + 40j
    assert e.input_impedance(zl, z0, np.pi / 2) == pytest.approx(z0**2 / zl)
    assert e.input_impedance(zl, z0, np.pi) == pytest.approx(zl)
    # Against the tangent formula.
    bl = np.linspace(0, 3, 7)
    expected = z0 * (zl + 1j * z0 * np.tan(bl)) / (z0 + 1j * zl * np.tan(bl))
    np.testing.assert_allclose(e.input_impedance(zl, z0, bl), expected)
    # A shorted stub is a pure reactance j Z0 tan(beta l).
    np.testing.assert_allclose(e.input_impedance(0, z0, 0.3), 1j * z0 * np.tan(0.3))
    for fn, args in [
        (e.reflection_coefficient, (-1 + 0j, z0)),
        (e.reflection_coefficient, (np.inf, z0)),
        (e.reflection_coefficient, (50, 0)),
        (e.standing_wave_ratio, (1.0,)),
        (e.input_impedance, (50, z0, -1)),
        (e.propagation_constant, (-1, 1, 0, 1, 1)),
        (e.distortionless_condition, (1, 0, 1, 1)),
    ]:
        with pytest.raises(ValueError):
            fn(*args)


def test_symmetrical_components_round_trip_and_balanced_set():
    a = np.exp(2j * np.pi / 3)
    v0, v1, v2 = e.symmetrical_components(1, a**2, a)
    assert abs(v0) < 1e-15 and v1 == pytest.approx(1) and abs(v2) < 1e-15
    va, vb, vc = np.array([1.0, 0.0, 0.0]), np.array([0.5 * a**2, 0.0, 1.0]), np.array([0.2 * a, 0.0, 1.0])
    seq = e.symmetrical_components(va, vb, vc)
    back = e.phase_components(*seq)
    for x, y in zip(back, (va, vb, vc)):
        np.testing.assert_allclose(x, y, atol=1e-15)
    # Zero sequence is the average of the phases.
    np.testing.assert_allclose(seq[0], (va + vb + vc) / 3)
    assert e.phase_components(1, 0, 0) == pytest.approx((1, 1, 1))
    with pytest.raises(ValueError):
        e.symmetrical_components(1, np.nan, 0)


def test_admittance_matrix_and_validation():
    Y = e.admittance_matrix(3, [(0, 1, 0.1j, 0.02j), (1, 2, 0.05 + 0.1j)])
    assert np.allclose(Y, Y.T)
    assert Y[0, 1] == pytest.approx(10j)
    assert Y[0, 0] == pytest.approx(-10j + 0.01j)
    np.testing.assert_allclose(Y.sum(axis=1), [0.01j, 0.01j, 0])
    for args in [(1, [(0, 1, 1j)]), (2, []), (2, [(0, 1)]), (2, [(0, 2, 1j)]), (2, [(1, 1, 1j)]), (2, [(0, 1, 0)]), (2, [(0, 1, np.nan)])]:
        with pytest.raises(ValueError):
            e.admittance_matrix(*args)


def test_two_bus_load_flow_matches_closed_form():
    # Receiving-end voltage of a lossless line feeding P + jQ: solve V^4 + ... exactly.
    x, P, Q = 0.1, 1.0, 0.5
    Y = e.admittance_matrix(2, [(0, 1, 1j * x)])
    r = e.gauss_seidel_load_flow(Y, [P], [Q])
    assert r.converged
    V2 = abs(r.voltages[1])
    # |V1|^2 = (V2 + Q x / V2)^2 + (P x / V2)^2 with V1 = 1.
    assert (V2 + Q * x / V2) ** 2 + (P * x / V2) ** 2 == pytest.approx(1, abs=1e-7)
    assert r.power[1] == pytest.approx(-(P + 1j * Q), abs=1e-6)
    assert r.losses.real == pytest.approx(0, abs=1e-6)
    assert r.losses.imag == pytest.approx(x * (P**2 + Q**2) / V2**2, rel=1e-5)
    summary = r.summary()
    assert summary["converged"] and summary["slack_power"] == r.power[0]


def test_load_flow_acceleration_nonconvergence_and_validation():
    # On a long radial feeder over-relaxation cuts the sweep count several-fold.
    feeder = e.admittance_matrix(8, [(k, k + 1, 0.01 + 0.03j) for k in range(7)])
    slow = e.gauss_seidel_load_flow(feeder, [0.1] * 7, [0.05] * 7)
    fast = e.gauss_seidel_load_flow(feeder, [0.1] * 7, [0.05] * 7, acceleration=1.6)
    assert slow.converged and fast.converged and fast.iterations < slow.iterations / 3
    np.testing.assert_allclose(fast.voltages, slow.voltages, atol=1e-6)

    Y = e.admittance_matrix(3, [(0, 1, 0.02 + 0.06j), (0, 2, 0.08 + 0.24j), (1, 2, 0.06 + 0.18j)])
    plain = e.gauss_seidel_load_flow(Y, [0.4, 0.45], [0.2, 0.15])
    assert plain.converged
    assert plain.losses.real > 0
    # Power balance: the slack supplies the loads plus the losses.
    assert plain.power[0].real == pytest.approx(0.85 + plain.losses.real)
    short = e.gauss_seidel_load_flow(Y, [0.4, 0.45], [0.2, 0.15], max_iter=2)
    assert not short.converged and short.iterations == 2 and short.history.size == 2
    for kwargs in [
        {"Y": np.ones(3)},
        {"Y": np.zeros((3, 3))},
        {"p_load": [0.4]},
        {"v_slack": 0},
        {"tol": 0},
        {"max_iter": 0},
        {"acceleration": 2.0},
    ]:
        args = {"Y": Y, "p_load": [0.4, 0.45], "q_load": [0.2, 0.15]} | kwargs
        with pytest.raises(ValueError):
            e.gauss_seidel_load_flow(**args)


def test_series_rlc_properties_and_phasors():
    c = e.SeriesRLC(10.0, 10e-3, 1e-6)
    assert c.quality_factor == pytest.approx(10)
    assert c.impedance(2 * np.pi * c.resonance_frequency) == pytest.approx(10)
    assert e.SeriesRLC(0, 1, 1).quality_factor == np.inf
    assert e.polar(e.phasor(2.0, 0.5)) == pytest.approx((2.0, 0.5))
    mag, ph = e.polar([1j, -1])
    np.testing.assert_allclose(mag, [1, 1])
    np.testing.assert_allclose(ph, [np.pi / 2, np.pi])
    with pytest.raises(ValueError):
        e.phasor(-1, 0)
    with pytest.raises(ValueError):
        e.power_factor(11, 10)
