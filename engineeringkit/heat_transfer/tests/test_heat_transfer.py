import numpy as np

import engineeringkit.heat_transfer as m


def test_conduction_and_lumped():
    assert np.isclose(m.fourier_heat_rate(10, 2, 20, 0.5), 800)
    assert np.isclose(m.lumped_temperature(400, 300, 10, 1, 1000, 1000, 0.01, 0), 400)


def test_composite_wall_energy():
    r = m.composite_wall(400, 300, [(0.1, 10), (0.2, 1)], A=1)
    assert np.isclose(r.heat_rate, 100 / (0.01 + 0.2))
    assert np.isclose(r.interface_temperatures[-1], 300)
