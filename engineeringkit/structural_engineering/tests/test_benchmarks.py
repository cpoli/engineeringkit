import numpy as np
import pytest
from scipy.integrate import quad

from engineeringkit import structural_engineering as s


def test_castigliano_recovers_axial_displacement():
    P, L, A, E = 1000.0, 2.0, 0.01, 200e9
    step = 1.0
    dU = (s.strain_energy_axial(P + step, L, A, E) - s.strain_energy_axial(P - step, L, A, E)) / (2 * step)
    assert np.isclose(dU, s.axial_deformation(P, L, A, E), rtol=1e-10)


def test_offcenter_beam_equilibrium_and_curvature():
    beam = s.SimplySupportedBeam(6, 200e9, 1e-5)
    r = beam.solve_point_load(12000, position=2, n=6001)
    ra, rb = r.reactions
    assert np.isclose(ra + rb, 12000)
    assert np.isclose(rb * 6, 12000 * 2)
    np.testing.assert_allclose(r.deflection[[0, -1]], 0, atol=1e-15)
    curvature = np.gradient(np.gradient(r.deflection, r.x), r.x)
    np.testing.assert_allclose(curvature[3:-3], r.moment[3:-3] / (200e9 * 1e-5), atol=2e-6)


def test_uniform_beam_energy_and_midspan():
    L, E, I, w = 4.0, 200e9, 2e-5, 3000.0
    result = s.SimplySupportedBeam(L, E, I).solve_uniform_load(w, n=1001)
    assert np.isclose(result.deflection[500], -5 * w * L**4 / (384 * E * I))
    energy = quad(lambda x: (w * x * (L - x) / 2) ** 2 / (2 * E * I), 0, L)[0]
    assert np.isclose(energy, w * w * L**5 / (240 * E * I))


def test_mohr_stress_invariants():
    c, r, maximum, minimum = s.mohr_circle(100, -20, 35)
    assert np.isclose(maximum + minimum, 80)
    assert np.isclose(maximum * minimum, 100 * (-20) - 35**2)
    assert np.isclose(r, (maximum - minimum) / 2)


def test_single_bar_equilibrium_and_strain():
    truss = s.PlanarTruss(np.array([[0.0, 0.0], [2.0, 0.0]]), ((0, 1),), 200e9, 0.01)
    r = truss.solve([0, 0, 1000, 0], [0, 1, 3])
    assert np.isclose(r.displacements[2], 1e-6)
    assert np.isclose(r.member_forces[0], 1000)
    assert np.isclose(r.reactions[0], -1000)
    assert np.isclose(r.reactions.sum() + 1000, 0)


@pytest.mark.parametrize("fixed", [[-1, 0], [0.5, 1], [0, 4]])
def test_truss_rejects_invalid_dofs(fixed):
    truss = s.PlanarTruss(np.array([[0.0, 0.0], [2.0, 0.0]]), ((0, 1),), 200e9, 0.01)
    with pytest.raises(ValueError):
        truss.solve([0, 0, 1000, 0], fixed)


def test_unrestrained_truss_has_clear_error():
    truss = s.PlanarTruss(np.array([[0.0, 0.0], [2.0, 0.0]]), ((0, 1),), 200e9, 0.01)
    with pytest.raises(ValueError, match="unstable"):
        truss.solve([0, 0, 1000, 0], [0])
