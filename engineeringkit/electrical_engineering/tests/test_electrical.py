import numpy as np

from engineeringkit.electrical_engineering import *


def test_ohm():
    assert ohms_law(voltage=10, resistance=5) == 2


def test_equivalents():
    assert series_resistance(2, 3) == 5 and np.isclose(parallel_resistance(2, 2), 1)
    assert np.isclose(capacitors_series(2, 2), 1) and capacitors_parallel(2, 3) == 5


def test_resonance():
    assert np.isclose(resonance_frequency(1, 1), 1 / (2 * np.pi))


def test_max_power():
    assert maximum_power_transfer(10, 5) == (5, 5)


def test_rc_limit():
    assert np.isclose(rc_charging_voltage(5, 1000, 1e-3, 100)[()], 5)
