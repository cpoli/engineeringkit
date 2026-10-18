"""Analytical checks for the belt, bearing, linkage, robot, absorber, Holzer and Jeffcott models."""

import numpy as np
import pytest
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.linalg import eigh

from engineeringkit import mechanical_engineering as m
from engineeringkit.mechanical_engineering.utils import bracketed_roots


def test_belt_tension_ratio_is_exponential_in_wrap_and_wedge():
    assert m.belt_tension_ratio(0.0, 5.0) == 1
    assert np.isclose(m.belt_tension_ratio(0.25, 2 * np.pi), np.exp(np.pi / 2))
    np.testing.assert_allclose(m.belt_tension_ratio([0.1, 0.2], np.pi), np.exp([0.1 * np.pi, 0.2 * np.pi]))
    assert np.isclose(m.belt_tension_ratio(0.2, 1.0, groove_angle=np.pi / 3), np.exp(0.4))
    for bad in ({"mu": -0.1, "wrap_angle": 1}, {"mu": 0.1, "wrap_angle": 1, "groove_angle": np.pi}):
        with pytest.raises(ValueError):
            m.belt_tension_ratio(**bad)


def test_petroff_torque_matches_couette_shear():
    viscosity, r, L, c, w = 0.05, 0.02, 0.04, 2e-5, 100.0
    shear = viscosity * w * r / c
    assert np.isclose(m.petroff_friction_torque(viscosity, r, L, c, w), shear * 2 * np.pi * r * L * r)
    assert m.petroff_friction_torque(viscosity, r, L, c, -w) < 0
    with pytest.raises(ValueError):
        m.petroff_friction_torque(viscosity, r, L, 0.0, w)


def test_sommerfeld_number_and_eccentricity_match_reynolds_equation():
    viscosity, omega, r, L, c, eps = 0.04, 250.0, 0.025, 1.0, 30e-6, 0.6
    theta = np.linspace(0, 2 * np.pi, 20001)
    h = c * (1 + eps * np.cos(theta))
    h_star = trapezoid(h**-2, theta) / trapezoid(h**-3, theta)
    p = cumulative_trapezoid(6 * viscosity * omega * r**2 * (h - h_star) / h**3, theta, initial=0)
    load = np.hypot(trapezoid(p * np.cos(theta), theta), trapezoid(p * np.sin(theta), theta)) * r * L
    S = m.sommerfeld_number(viscosity, omega, load, r, L, c)
    assert np.isclose(S, (2 + eps**2) * np.sqrt(1 - eps**2) / (12 * np.pi**2 * eps), rtol=1e-6)
    assert np.isclose(m.long_bearing_eccentricity(S), eps, rtol=1e-6)
    eccentricities = m.long_bearing_eccentricity(np.array([[1e-3, 0.1], [1.0, 100.0]]))
    assert eccentricities.shape == (2, 2)
    assert np.all(np.diff(eccentricities.ravel()) < 0) and np.all((eccentricities > 0) & (eccentricities < 1))
    with pytest.raises(ValueError):
        m.long_bearing_eccentricity(0.0)


def test_screw_rejects_impossible_angle_for_any_array_element():
    np.testing.assert_allclose(2 * np.pi * m.screw_torque_raise([10.0, 20.0], 0.02, 0.01, 0.0), [0.1, 0.2])
    with pytest.raises(ValueError):
        m.screw_torque_raise(100, 0.01, 1.0, [0.0, 100.0])


def test_fourbar_output_angle_closes_the_loop_on_both_circuits():
    a, b, c, d = 1.0, 3.0, 2.0, 3.0
    theta2 = np.linspace(0, 2 * np.pi, 37)
    for crossed in (False, True):
        theta4 = m.fourbar_output_angle(a, b, c, d, theta2, crossed=crossed)
        coupler = np.hypot(d + c * np.cos(theta4) - a * np.cos(theta2), c * np.sin(theta4) - a * np.sin(theta2))
        np.testing.assert_allclose(coupler, b)
    assert np.sin(m.fourbar_output_angle(a, b, c, d, np.pi / 2)) > 0
    assert isinstance(m.fourbar_output_angle(a, b, c, d, 0.3), float)
    with pytest.raises(ValueError):
        m.fourbar_output_angle(3.0, 5.0, 9.0, 4.0, theta2)


def test_freudenstein_synthesis_round_trip():
    lengths = np.array([1.0, 3.0, 2.0, 3.0])
    t2 = np.radians([20.0, 45.0, 80.0, 120.0])
    t4 = m.fourbar_output_angle(*lengths, t2)
    K = m.freudenstein_coefficients(t2, t4)
    np.testing.assert_allclose(K, [3.0, 1.5, 1.25])
    np.testing.assert_allclose(m.freudenstein_link_lengths(K, ground=3.0), lengths)
    with pytest.raises(ValueError):
        m.freudenstein_coefficients(t2[:2], t4[:2])
    with pytest.raises(ValueError):
        m.freudenstein_coefficients(t2, t4[:3])
    with pytest.raises(ValueError):
        m.freudenstein_coefficients([0.0, 0.0, 0.0], [0.1, 0.1, 0.1])
    for bad in ([3.0, 1.5], [3.0, 1.5, 1.25, 1.0], [-3.0, 1.5, 1.25], [1.0, 1.0, 5.0]):
        with pytest.raises(ValueError):
            m.freudenstein_link_lengths(bad)


def test_dh_forward_kinematics_of_planar_and_spatial_arms():
    T = m.dh_transform(0.5, np.pi / 2, 0.2, np.pi / 2)
    np.testing.assert_allclose(T[:3, 3], [0, 0.5, 0.2], atol=1e-15)
    np.testing.assert_allclose(T[:3, :3] @ T[:3, :3].T, np.eye(3), atol=1e-15)
    frames = m.forward_kinematics([[1.0, 0, 0, np.pi / 2], [1.0, 0, 0, -np.pi / 2]])
    assert frames.shape == (3, 4, 4)
    np.testing.assert_allclose(frames[-1][:3, 3], [1, 1, 0], atol=1e-15)
    q1, q2 = 0.4, 0.7
    hand = m.forward_kinematics([[0, np.pi / 2, 0.3, q1], [0.5, 0, 0, q2]])[-1][:3, 3]
    np.testing.assert_allclose(hand, [0.5 * np.cos(q2) * np.cos(q1), 0.5 * np.cos(q2) * np.sin(q1), 0.3 + 0.5 * np.sin(q2)])
    for bad in ([], [[1, 2, 3]], [1, 2, 3, 4]):
        with pytest.raises(ValueError):
            m.forward_kinematics(bad)
    with pytest.raises(ValueError):
        m.dh_transform(np.nan, 0, 0, 0)


def test_tuned_mass_damper_frahm_zero_and_den_hartog_fixed_points():
    frahm = m.tuned_mass_damper(10.0, 1e4, 0.2, frequency_ratio=1.0, damping_ratio=0.0)
    assert frahm.absorber_mass == 2 and np.isclose(frahm.absorber_stiffness, 2e3) and frahm.absorber_damping == 0
    assert np.isclose(frahm.primary_frequency, np.sqrt(1e3))
    assert np.isclose(frahm.amplification(1.0), 0, atol=1e-12)
    X1, X2 = frahm.response(1.0)
    assert np.isclose(abs(X2), 1 / 0.2)  # the absorber spring force cancels the applied force
    assert np.isclose(frahm.amplification(0.0), 1)
    mu = 0.2
    fixed = np.sqrt((1 + np.array([-1, 1]) * np.sqrt(mu / (2 + mu))) / (1 + mu))
    for zeta in (0.05, 0.3):
        np.testing.assert_allclose(m.tuned_mass_damper(10.0, 1e4, mu, damping_ratio=zeta).amplification(fixed), np.sqrt(1 + 2 / mu))
    optimum = m.tuned_mass_damper(10.0, 1e4, mu)
    assert np.isclose(optimum.frequency_ratio, 1 / 1.2)
    assert np.isclose(optimum.amplification(np.linspace(0.5, 1.5, 4001)).max(), np.sqrt(1 + 2 / mu), rtol=0.01)
    summary = optimum.summary()
    assert summary["mass_ratio"] == mu and np.isclose(summary["fixed_point_amplification"], np.sqrt(11))
    for bad in ({"mass_ratio": 0}, {"mass_ratio": 0.1, "frequency_ratio": 0}, {"mass_ratio": 0.1, "damping_ratio": -1}):
        with pytest.raises(ValueError):
            m.tuned_mass_damper(10.0, 1e4, **bad)
    with pytest.raises(ValueError):
        optimum.response(-1.0)


def test_holzer_matches_matrix_eigenvalues():
    J, K = [2.0, 1.0, 3.0, 0.5], [4e3, 2e3, 6e3]
    result = m.holzer_frequencies(J, K, 300.0)
    stiffness = np.zeros((4, 4))
    for i, k in enumerate(K):
        stiffness[i : i + 2, i : i + 2] += k * np.array([[1, -1], [-1, 1]])
    np.testing.assert_allclose(result.frequencies, np.sqrt(eigh(stiffness, np.diag(J), eigvals_only=True)[1:]), rtol=1e-10)
    np.testing.assert_allclose(np.max(np.abs(result.mode_shapes), axis=1), 1)
    np.testing.assert_allclose(m.holzer_residual(result.frequencies, J, K), 0, atol=1e-6)
    assert isinstance(m.holzer_residual(10.0, J, K), float)
    assert result.summary()["modes"] == 3
    empty = m.holzer_frequencies(J, K, 1.0)
    assert empty.frequencies.size == 0 and empty.mode_shapes.shape == (0, 4)
    for args in (([1.0], []), ([1.0, 1.0], [1.0, 1.0]), ([1.0, -1.0], [1.0])):
        with pytest.raises(ValueError):
            m.holzer_residual(1.0, *args)
    with pytest.raises(ValueError):
        m.holzer_frequencies(J, K, 300.0, n_grid=5)


def test_jeffcott_whirl_limits_and_bearing_force():
    mass, k, e, c = 10.0, 1e6, 1e-4, 200.0
    wc = np.sqrt(k / mass)
    result = m.jeffcott_response(mass, k, e, np.array([0.0, wc, 100 * wc]), c)
    zeta = c / (2 * np.sqrt(k * mass))
    assert result.amplitude[0] == 0 and np.isclose(result.amplitude[1], e / (2 * zeta))
    assert np.isclose(result.amplitude[2], e, rtol=1e-3) and np.isclose(result.phase[1], np.pi / 2)
    assert result.phase[2] > 0.99 * np.pi
    w = np.linspace(0, 3 * wc, 31)
    r = m.jeffcott_response(mass, k, e, w, c)
    bearing = np.abs(k + 1j * c * w) * r.amplitude
    np.testing.assert_allclose(bearing, m.rotating_unbalance_force(mass, e, w) * m.transmissibility(w / wc, zeta))
    undamped = m.jeffcott_response(mass, k, e, wc)
    assert np.isinf(undamped.amplitude) and np.isclose(undamped.summary()["critical_speed_rad_s"], wc)
    with pytest.raises(ValueError):
        m.jeffcott_response(mass, k, e, -1.0)


def test_bracketed_roots_finds_every_sign_change_and_exact_zeros():
    roots = bracketed_roots(np.sin, np.linspace(0.5, 10, 50))
    np.testing.assert_allclose(roots, np.pi * np.arange(1, 4))
    np.testing.assert_allclose(bracketed_roots(lambda x: x - 1, np.array([0.0, 1.0, 2.0])), [1.0])


def test_vectorised_scalar_helpers():
    np.testing.assert_allclose(m.damped_frequency(4.0, 1.0, np.array([0.0, 4.0, 8.0])), [2.0, 0.0, 0.0])
    np.testing.assert_allclose(m.slider_crank_displacement(1.0, 3.0, np.array([0.0, np.pi])), [4.0, 2.0])
    assert m.fourbar_grashof(np.array([1.0, 2.0, 2.0, 3.0])) is True
    with pytest.raises(ValueError):
        m.work_constant_force(np.inf, 1.0)
    with pytest.raises(ValueError):
        m.power_torque(1.0, np.nan)
