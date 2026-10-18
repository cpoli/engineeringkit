import numpy as np
import pytest

from engineeringkit.mechanical_engineering import simulate_sdof


def test_undamped_energy_conservation():
    r = simulate_sdof(2, 0, 8, x0=0.1, t=np.linspace(0, 10, 501))
    np.testing.assert_allclose(r.displacement, 0.1 * np.cos(2 * r.t), atol=1e-8)
    energy = r.velocity**2 + 4 * r.displacement**2
    np.testing.assert_allclose(energy, 0.04, rtol=1e-6)
    assert r.history.shape == (501, 2)


def test_critical_damping_analytic_solution():
    r = simulate_sdof(1, 4, 4, x0=1, t=np.linspace(0, 5, 101))
    np.testing.assert_allclose(r.displacement, (1 + 2 * r.t) * np.exp(-2 * r.t), atol=1e-8)


@pytest.mark.parametrize("properties", [(0, 1, 1), (1, -1, 1), (1, 1, -1), (np.nan, 1, 1)])
def test_invalid_passive_oscillator(properties):
    with pytest.raises(ValueError):
        simulate_sdof(*properties)
