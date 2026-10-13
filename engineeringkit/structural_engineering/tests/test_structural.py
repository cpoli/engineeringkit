import numpy as np

from engineeringkit.structural_engineering import *


def test_hooke_and_section():
    assert hooke_stress(200e9, 1e-3) == 2e8
    assert np.isclose(rectangular_section_modulus(0.2, 0.3), 0.003)


def test_beam_midpoint():
    b = SimplySupportedBeam(4, 200e9, 8e-6)
    r = b.solve_point_load(1000, 2)
    assert np.isclose(max(r.moment), 1000)
    assert np.isclose(min(r.deflection), -1000 * 4**3 / (48 * 200e9 * 8e-6), rtol=2e-5)


def test_energy_and_reciprocity():
    assert np.isclose(strain_energy_axial(1000, 2, 0.01, 200e9), 0.0005)
    assert maxwell_betti_reciprocity(1e-3, 1e-3)


def test_triangle_truss_equilibrium():
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]])
    tr = PlanarTruss(nodes, ((0, 1), (0, 2), (1, 2)), 200e9, 1e-3)
    r = tr.solve([[0, 0], [0, 0], [0, -1000]], [0, 1, 3])
    assert np.isclose(r.reactions[1] + r.reactions[3], 1000)
    assert np.all(np.isfinite(r.member_forces))
