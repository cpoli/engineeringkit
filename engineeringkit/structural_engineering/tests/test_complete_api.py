import numpy as np
import pytest

from engineeringkit import structural_engineering as s


def test_elastic_section_and_torsion_values():
    assert s.thermal_strain(1e-5, 100) == 0.001
    assert s.poisson_lateral_strain(0.01, 0.3) == -0.003
    assert np.isclose(s.euler_buckling_load(1, 1, 1), np.pi**2)
    assert s.radius_of_gyration(4, 1) == 2
    assert s.slenderness_ratio(10, 1, 2) == 5
    assert s.torsion_shear_stress(10, 2, 4) == 5
    assert s.torsion_angle(10, 2, 5, 4) == 1
    assert s.rectangular_second_moment(12, 1) == 1
    assert s.parallel_axis(1, 2, 3) == 19
    assert s.bending_stress(10, 2, 4) == 5
    assert s.castigliano_deflection(0.01) == 0.01
    assert s.strain_energy_bending_constant_moment(10, 2, 5, 4) == 5


def test_beam_plot_summary_and_load_position():
    import matplotlib.pyplot as plt

    beam = s.SimplySupportedBeam(4, 200e9, 1e-5)
    r = beam.solve_point_load(100)
    assert r.summary()["reactions"] == (50, 50)
    fig, axs = r.plot()
    assert len(axs) == 3
    plt.close(fig)
    with pytest.raises(ValueError):
        beam.solve_point_load(100, position=5)
    edge = beam.solve_point_load(100, position=0)
    np.testing.assert_allclose(edge.deflection, 0)


def test_truss_plot_summary_and_invalid_meshes():
    import matplotlib.pyplot as plt

    xy = np.array([[0.0, 0.0], [2.0, 0.0]])
    truss = s.PlanarTruss(xy, ((0, 1),), 200e9, 0.01)
    r = truss.solve([0, 0, 1000, 0], [0, 1, 3])
    assert r.summary()["max_abs_member_force"] == pytest.approx(1000)
    fig, ax = r.plot(xy, ((0, 1),), scale=100)
    assert len(ax.lines) == 2
    plt.close(fig)
    for nodes, members, E, loads in [
        ([0, 1], ((0, 1),), 200e9, [0, 0, 0, 0]),
        (xy, ((0, 1),), 0, [0, 0, 0, 0]),
        (xy, (), 200e9, [0, 0, 0, 0]),
        (xy, ((0, 1),), 200e9, [0, 0]),
        (xy, ((0, 2),), 200e9, [0, 0, 0, 0]),
        (xy, ((0, 0),), 200e9, [0, 0, 0, 0]),
    ]:
        with pytest.raises(ValueError):
            s.PlanarTruss(nodes, members, E, 0.01).solve(loads, [0, 1, 3])
