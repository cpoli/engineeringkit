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
