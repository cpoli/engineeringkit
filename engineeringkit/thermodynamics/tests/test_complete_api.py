import numpy as np
import pytest

from engineeringkit import thermodynamics as t


def test_first_law_properties_exergy_and_cops():
    assert t.first_law_heat(10, 20) == 30
    assert t.isentropic_temperature_ratio(1) == 1
    with pytest.raises(ValueError):
        t.isentropic_temperature_ratio(2, 1)
    assert t.specific_internal_energy_ideal_gas(300, 700) == 210000
    assert t.specific_enthalpy_ideal_gas(300, 1000) == 300000
    assert t.exergy_heat(100, 600, 300) == 50
    with pytest.raises(ValueError):
        t.heat_engine_work(1, 2)
    assert t.thermal_efficiency(100, 50) == 0.5
    assert np.isclose(t.isentropic_pressure_volume(100, 1, 2, 1.4), 100 / 2**1.4)
    assert t.refrigeration_cop(100, 50) == 2
    assert t.heat_pump_cop(150, 50) == 3
    cycle = t.OttoCycle(8)
    assert cycle.summary()["efficiency"] == t.otto_efficiency(8)


def test_gamma_and_ratios_vectorize():
    gamma = np.array([1.3, 1.4, 5 / 3])
    for f, args in [
        (t.isentropic_pressure_volume, (100, 1, 2)),
        (t.isentropic_temperature_ratio, (3,)),
        (t.otto_efficiency, (8,)),
        (t.diesel_efficiency, (8, 2)),
        (t.brayton_efficiency, (6,)),
    ]:
        assert np.allclose(f(*args, gamma), [f(*args, g) for g in gamma])
    cutoff = np.array([1.0, 1.5, 2.0])
    assert np.allclose(t.diesel_efficiency(8, cutoff), [t.diesel_efficiency(8, c) for c in cutoff])
    assert t.diesel_efficiency(8, cutoff)[0] == t.otto_efficiency(8)
    assert np.allclose(t.brayton_efficiency([1, 6]), [0, t.brayton_efficiency(6)])
    with pytest.raises(ValueError):
        t.otto_efficiency(8, [1.4, 1.0])
    with pytest.raises(ValueError):
        t.diesel_efficiency([8, 2], 3)
    with pytest.raises(ValueError):
        t.brayton_efficiency([0.5, 6])
