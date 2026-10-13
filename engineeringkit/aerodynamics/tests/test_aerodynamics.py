import numpy as np

import engineeringkit.aerodynamics as m


def test_thin_airfoil_and_shock():
    assert np.isclose(m.thin_airfoil_cl(0.1), 2 * np.pi * 0.1)
    assert m.normal_shock_mach2(2) < 1


def test_drag_polar():
    assert m.drag_polar(0, 0.02, 8) == 0.02
