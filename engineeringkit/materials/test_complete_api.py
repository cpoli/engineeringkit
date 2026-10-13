import numpy as np
import pytest

from engineeringkit import materials as m


def test_stress_strain_and_thermal_constraint():
    assert m.normal_stress(100, 2) == 50
    assert m.normal_strain(0.01, 1) == 0.01
    np.testing.assert_allclose(m.hooke_stress(100, [0.01, 0.02]), [1, 2])
    lam, G = m.lame_parameters(100, 0.25)
    assert lam == 40 and G == 40
    assert m.thermal_strain(1e-5, 100) == 0.001
    assert m.constrained_thermal_stress(100, 1e-5, 100) == -0.1
    assert m.factor_of_safety(100, 20) == 5
    assert np.isclose(m.griffith_fracture_stress(np.pi, 2, 1), 2)
    assert m.paris_crack_growth(2, 3, 2) == 16
    assert m.goodman_safety_factor(10, 20, 100, 200) == 5
    with pytest.raises(ValueError):
        m.goodman_safety_factor(0, 0, 100, 200)
    assert np.isclose(m.arrhenius_creep_rate(1, 1, 1, 1, 2), 2 / np.e)
    assert m.hall_petch(100, 2, 4) == 101
    assert m.density_rule_of_mixtures(100, 200, 0.5) == 150


def test_curve_plot_and_invalid_fit_window():
    import matplotlib.pyplot as plt

    curve = m.StressStrainCurve([0, 0.1, 0.2], [0, 1, 2])
    with pytest.raises(ValueError):
        curve.young_modulus(4)
    fig, ax = curve.plot()
    assert ax.get_ylabel() == "stress [Pa]"
    plt.close(fig)
