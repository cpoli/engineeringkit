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
