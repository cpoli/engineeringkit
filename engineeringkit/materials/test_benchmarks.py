import numpy as np
import pytest

from engineeringkit import materials as m


def test_isotropic_elastic_identity_and_auxetic_range():
    E, nu = 200e9, 0.3
    G, K = m.shear_modulus(E, nu), m.bulk_modulus(E, nu)
    assert np.isclose(9 * K * G / (3 * K + G), E)
    assert m.shear_modulus(E, -0.2) > 0


def test_toughness_elastic_triangle_and_slope():
    strain = np.linspace(0, 0.001, 11)
    stress = 200e9 * strain
    curve = m.StressStrainCurve(strain, stress)
    assert np.isclose(curve.young_modulus(), 200e9)
    assert np.isclose(curve.toughness(), 0.5 * stress[-1] * strain[-1])
    assert np.isclose(curve.toughness(), m.resilience_modulus(stress[-1], 200e9))


def test_fracture_inverse_and_composite_bounds():
    K, sigma, Y = 50e6, 100e6, 1.12
    length = m.critical_crack_length(K, sigma, Y)
    assert np.isclose(m.stress_intensity(Y, sigma, length), K)
    for fraction in [0, 0.25, 0.5, 1]:
        voigt = m.rule_of_mixtures_E_longitudinal(200e9, 3e9, fraction)
        reuss = m.rule_of_mixtures_E_transverse(200e9, 3e9, fraction)
        assert 3e9 <= reuss <= voigt <= 200e9


def test_yield_invariance_and_fatigue_accumulation():
    assert m.von_mises_plane_stress(100, 100, 0) == 100
    assert np.isclose(m.von_mises_plane_stress(0, 0, 100), np.sqrt(3) * 100)
    assert m.tresca_equivalent([-50, 0, 50]) == 100
    assert m.miner_damage([100, 200], [1000, 400]) == 0.6
    assert m.basquin_cycles(200, 1000, -0.1) > m.basquin_cycles(300, 1000, -0.1)


def test_curve_and_damage_shapes():
    with pytest.raises(ValueError):
        m.StressStrainCurve([0, 0.1, 0.05], [0, 1, 2])
    with pytest.raises(ValueError):
        m.miner_damage([1, 2], [3])
