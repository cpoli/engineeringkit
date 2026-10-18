"""Analytical checks for the elasticity, fatigue, fracture, creep and plasticity models added in the rebuild."""

import numpy as np
import pytest
from scipy.integrate import quad

from engineeringkit import materials as m
from engineeringkit.materials.utils import power_law_fit


def test_lame_constants_reproduce_every_isotropic_modulus():
    E, nu = np.array([70e9, 200e9]), np.array([0.33, 0.3])
    lam, mu = m.lame_parameters(E, nu)
    np.testing.assert_allclose(mu * (3 * lam + 2 * mu) / (lam + mu), E)
    np.testing.assert_allclose(lam / (2 * (lam + mu)), nu)
    np.testing.assert_allclose(lam + 2 * mu / 3, m.bulk_modulus(E, nu))
    with pytest.raises(ValueError):
        m.shear_modulus(200e9, 0.5)


def test_duhamel_constrained_stress_cancels_the_free_strain():
    E, alpha, dT = 200e9, 12e-6, np.array([-40.0, 0.0, 80.0])
    stress = m.constrained_thermal_stress(E, alpha, dT)
    np.testing.assert_allclose(m.thermal_strain(alpha, dT) + stress / E, 0, atol=1e-18)
    assert stress[-1] < 0
    with pytest.raises(ValueError):
        m.thermal_strain(alpha, np.nan)


def test_tresca_and_von_mises_agree_in_uniaxial_tension_and_differ_in_shear():
    principal = np.array([[250.0, 0.0, 0.0], [125.0, 0.0, -125.0]])
    np.testing.assert_allclose(m.tresca_equivalent(principal), [250, 250])
    assert m.von_mises_plane_stress(250.0, 0.0, 0.0) == 250
    assert np.isclose(m.tresca_equivalent([0.0, 100.0, -100.0]) / m.von_mises_plane_stress(0, 0, 100.0), 2 / np.sqrt(3))
    with pytest.raises(ValueError):
        m.tresca_equivalent(100.0)
    with pytest.raises(ValueError):
        m.von_mises_plane_stress(np.inf, 0, 0)


def test_basquin_and_coffin_manson_invert_their_power_laws():
    N = m.basquin_cycles(np.array([300e6, 400e6]), 900e6, -0.09)
    np.testing.assert_allclose(900e6 * (2 * N) ** -0.09, [300e6, 400e6])
    N = m.coffin_manson_cycles(np.array([0.002, 0.01]), 0.6, -0.6)
    np.testing.assert_allclose(0.6 * (2 * N) ** -0.6, [0.002, 0.01])
    for bad in (0.0, [-0.1, 0.1]):
        with pytest.raises(ValueError):
            m.basquin_cycles(300e6, 900e6, bad)
    with pytest.raises(ValueError):
        m.goodman_safety_factor([0.0, 10.0], [0.0, 10.0], 100.0, 200.0)


def test_strain_life_curve_splits_at_the_transition_life():
    E, sf, b, ef, c = 200e9, 1000e6, -0.09, 0.6, -0.6
    Nt = 0.5 * (ef * E / sf) ** (1 / (b - c))
    total, elastic, plastic = m.strain_life_amplitude(Nt, E, sf, b, ef, c)
    assert np.isclose(elastic, plastic) and np.isclose(total, 2 * elastic)
    assert np.isclose(m.coffin_manson_cycles(plastic, ef, c), Nt)
    assert np.isclose(m.basquin_cycles(elastic * E, sf, b), Nt)


def test_griffith_stress_is_where_irwin_k_reaches_sqrt_2_e_gamma():
    E, gamma, a = 70e9, 1.0, 1e-6
    sigma_f = m.griffith_fracture_stress(E, gamma, a)
    assert np.isclose(m.stress_intensity(1.0, sigma_f, a), np.sqrt(2 * E * gamma))
    assert np.isclose(np.pi * sigma_f**2 * a / E, 2 * gamma)


def test_crack_growth_life_matches_closed_form_and_quadrature():
    C, n, ds, a0, af, Y = 1e-11 * 1e-18, 3.0, 200e6, 1e-3, 15e-3, 1.12
    closed = 2 / ((n - 2) * C * (Y * ds * np.sqrt(np.pi)) ** n) * (a0 ** (1 - n / 2) - af ** (1 - n / 2))
    result = m.crack_growth_life(C, n, ds, a0, af, Y=Y)
    assert np.isclose(result.life, closed, rtol=1e-5)
    assert np.isclose(m.crack_growth_life(C, n, ds, a0, af, Y=lambda a: Y).life, result.life)
    width = 0.05

    def secant(a):
        return np.sqrt(1 / np.cos(np.pi * a / (2 * width)))

    corrected = m.crack_growth_life(C, n, ds, a0, af, Y=secant)
    exact = quad(lambda a: 1 / (C * (secant(a) * ds * np.sqrt(np.pi * a)) ** n), a0, af)[0]
    assert np.isclose(corrected.life, exact, rtol=1e-5) and corrected.life < m.crack_growth_life(C, n, ds, a0, af).life
    summary = result.summary()
    assert summary["cycles_to_double"] < summary["life_cycles"]
    assert np.isclose(result.cycles_to(af), result.life)
    with pytest.raises(ValueError):
        result.cycles_to(2 * af)
    with pytest.raises(ValueError):
        m.crack_growth_life(C, n, ds, af, a0)
    with pytest.raises(ValueError):
        m.crack_growth_life(C, n, ds, a0, af, Y=lambda a: -a)
    short = m.crack_growth_life(C, n, ds, a0, 1.5 * a0)
    assert short.summary()["cycles_to_double"] == short.life


def test_larson_miller_trades_temperature_for_time():
    P = m.larson_miller_parameter(900.0, 3600.0 * 1e3)
    assert np.isclose(P, 900 * 23)
    assert np.isclose(m.larson_miller_rupture_time(P, 900.0), 3600.0 * 1e3)
    hotter = m.larson_miller_rupture_time(P, 950.0)
    assert np.isclose(m.larson_miller_parameter(950.0, hotter), P) and hotter < 3600.0 * 1e3
    rate = m.arrhenius_creep_rate(1e-30, 300e3, 8.314, np.array([800.0, 900.0]), 100e6, n=5)
    assert np.isclose(rate[1] / rate[0], np.exp(300e3 / 8.314 * (1 / 800 - 1 / 900)))


def test_dislocations_close_the_gap_between_ideal_and_real_strength():
    G, nu, b = 26e9, 0.35, 2.86e-10
    assert np.isclose(m.ideal_shear_strength(G, b, b), G / (2 * np.pi))
    a = np.array([0.8, 1.0, 1.2]) * b
    tau = m.peierls_stress(G, nu, b, a)
    np.testing.assert_allclose(tau, 2 * G / (1 - nu) * np.exp(-2 * np.pi * a / ((1 - nu) * b)))
    assert np.all(np.diff(tau) < 0) and m.ideal_shear_strength(G, b, b) / tau[1] > 100


def test_hall_petch_is_linear_in_inverse_root_grain_size():
    d = np.array([1e-6, 1e-5, 1e-4])
    sy = m.hall_petch(70e6, 0.74e6, d)
    slope = np.polyfit(d**-0.5, sy, 1)[0]
    assert np.isclose(slope, 0.74e6)


def test_ramberg_osgood_offset_yield_is_recovered_from_the_curve():
    E, sy, n = 70e9, 300e6, 15
    stress = np.linspace(0, 1.1 * sy, 2001)
    strain = m.ramberg_osgood_strain(stress, E, sy, n)
    assert np.isclose(m.ramberg_osgood_strain(sy, E, sy, n) - sy / E, 0.002)
    assert m.ramberg_osgood_strain(-sy, E, sy, n) == -m.ramberg_osgood_strain(sy, E, sy, n)
    curve = m.StressStrainCurve(strain, stress)
    assert np.isclose(curve.offset_yield_strength(), sy, rtol=2e-3)
    assert curve.summary()["ultimate_strength_Pa"] == stress[-1]
    with pytest.raises(ValueError):
        curve.offset_yield_strength(offset=0.5)
    with pytest.raises(ValueError):
        m.ramberg_osgood_strain(sy, E, sy, 0.5)


def test_power_law_fit_recovers_basquin_constants():
    sa = np.geomspace(200e6, 600e6, 7)
    N = m.basquin_cycles(sa, 900e6, -0.09)
    A, p = power_law_fit(N, sa)
    assert np.isclose(p, -0.09) and np.isclose(A, 900e6 * 2**-0.09)
    for bad in (([1.0], [1.0]), ([1.0, 2.0], [1.0]), ([2.0, 2.0], [1.0, 3.0])):
        with pytest.raises(ValueError):
            power_law_fit(*bad)
    with pytest.raises(ValueError):
        power_law_fit([1.0, -2.0], [1.0, 2.0])
