import numpy as np

import engineeringkit.materials as m


def test_elastic_and_fracture():
    assert np.isclose(m.shear_modulus(210e9, 0.3), 210e9 / 2.6)
    a = m.critical_crack_length(50e6, 100e6)
    assert np.isclose(m.stress_intensity(1, 100e6, a), 50e6)


def test_miner():
    assert np.isclose(m.miner_damage([10, 20], [100, 200]), 0.2)
