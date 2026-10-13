import numpy as np
import pytest

from engineeringkit import heat_transfer as h


def test_wall_energy_and_interface_continuity():
    result = h.composite_wall(400, 300, [(0.1, 1), (0.2, 2)], A=2, h_hot=10, h_cold=20)
    resistances = np.array([0.05, 0.05, 0.05, 0.025])
    assert np.isclose(result.heat_rate, 100 / resistances.sum())
    np.testing.assert_allclose(-np.diff(result.interface_temperatures), result.heat_rate * resistances)
    assert np.isclose(result.interface_temperatures[-1], 300)


def test_cylinder_thin_wall_limit():
    radius, thickness, k, L = 1.0, 1e-6, 20.0, 2.0
    planar = h.plane_wall_resistance(thickness, k, 2 * np.pi * radius * L)
    assert np.isclose(h.cylinder_resistance(radius, radius + thickness, k, L), planar, rtol=1e-6)


def test_lumped_cooling_initial_and_asymptotic():
    assert h.lumped_temperature(400, 300, 10, 2, 1000, 1000, 0.01, 0) == 400
    assert np.isclose(h.lumped_temperature(400, 300, 10, 2, 1000, 1000, 0.01, 1e6), 300)
    assert h.radiation_heat_rate(1, 1, 300, 300) == 0
    assert h.radiation_heat_rate(0, 1, 400, 300) == 0


def test_heat_exchanger_singular_and_small_ntu_limits():
    assert h.log_mean_temperature_difference(20, 20) == 20
    assert np.isclose(h.log_mean_temperature_difference(20 + 1e-8, 20), 20 + 5e-9, rtol=1e-12)
    for Cr in [0, 0.5, 1 - 1e-12, 1]:
        assert h.exchanger_effectiveness_counterflow(0, Cr) == 0
        assert np.isclose(h.exchanger_effectiveness_counterflow(1e-12, Cr), 1e-12, rtol=1e-10)
    assert np.isclose(h.exchanger_effectiveness_counterflow(2, 1), 2 / 3)
    assert h.exchanger_effectiveness_counterflow(2, 0.5) > h.exchanger_effectiveness_parallel(2, 0.5)


@pytest.mark.parametrize("kwargs", [{"h_hot": 0}, {"h_cold": -1}, {"A": 0}])
def test_wall_invalid_boundary_data(kwargs):
    with pytest.raises(ValueError):
        h.composite_wall(400, 300, [(0.1, 1)], **kwargs)
