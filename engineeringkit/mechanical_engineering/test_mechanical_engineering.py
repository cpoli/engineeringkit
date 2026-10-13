import numpy as np

import engineeringkit.mechanical_engineering as m


def test_vibration_and_gears():
    assert np.isclose(m.natural_frequency(100, 4), 5)
    assert m.gear_ratio(40, 20) == 2


def test_sdof_initial():
    r = m.simulate_sdof(1, 0, 4, x0=1, t=np.linspace(0, 1, 20))
    assert np.isclose(r.displacement[0], 1)
