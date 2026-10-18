"""Analytical checks for the continuous-beam, influence-line, energy, truss, plasticity, buckling and frame models."""

import warnings

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

from engineeringkit import structural_engineering as s
from engineeringkit.structural_engineering.utils import fixed_dof_vector, mesh_geometry, solve_reduced


def test_three_moment_equation_matches_textbook_continuous_beams():
    w, L = 1000.0, 4.0
    np.testing.assert_allclose(s.three_moment_equation([L], w), [0, 0])
    np.testing.assert_allclose(s.three_moment_equation([L, L, L], w), [0, -w * L**2 / 10, -w * L**2 / 10, 0])
    # one loaded span of two: the support moment is -wL^2/16
    np.testing.assert_allclose(s.three_moment_equation([L, L], [w, 0.0]), [0, -w * L**2 / 16, 0])
    # a stiffer second span draws more moment: compare with the plane-frame solution
    moments = s.three_moment_equation([L, 2 * L], w, I=[1.0, 3.0])
    frame = s.PlanarFrame([[0, 0], [L, 0], [3 * L, 0]], ((0, 1), (1, 2)), 1e9, 1e3, [1e-3, 3e-3])
    result = frame.solve(np.zeros(9), [0, 1, 4, 7], member_loads=w)
    assert np.isclose(moments[1], -result.member_end_forces[1, 2], rtol=1e-6)  # counterclockwise end moment = -sagging
    with pytest.raises(ValueError):
        s.three_moment_equation([L, -1.0])
    with pytest.raises(ValueError, match="real"):
        s.three_moment_equation([L + 1j, L])


def test_influence_lines_equal_point_load_responses():
    L, at = 10.0, 3.0
    beam = s.SimplySupportedBeam(L, 200e9, 1e-4)
    k = 300  # station at x = 3 m on 1001 stations
    for xi in (1.0, 3.0, 7.5):
        r = beam.solve_point_load(1.0, position=xi, n=1001)
        assert np.isclose(s.influence_line(L, xi, "moment", at=at), r.moment[k])
        assert np.isclose(s.influence_line(L, xi, "shear", at=at), r.shear[k])
        assert np.isclose(s.influence_line(L, xi, "left_reaction"), r.reactions[0])
        assert np.isclose(s.influence_line(L, xi, "right_reaction"), r.reactions[1])
    for kwargs in ({"quantity": "torque"}, {"quantity": "moment"}, {"quantity": "shear", "at": 11.0}):
        with pytest.raises(ValueError):
            s.influence_line(L, 1.0, **kwargs)
    with pytest.raises(ValueError):
        s.influence_line(L, 12.0, "left_reaction")


def test_castigliano_integrates_cantilever_and_simply_supported_deflections():
    L, E, I, P = 3.0, 200e9, 4e-6, 1000.0
    x = np.linspace(0, L, 101)
    tip = s.castigliano_deflection(x=x, moment=-P * (L - x), dmoment_dP=-(L - x), E=E, I=I)
    assert np.isclose(tip, P * L**3 / (3 * E * I))
    # midspan of a simply supported beam under a uniform load, with a dummy load
    w = 2000.0
    M = w * x * (L - x) / 2
    dM = np.where(x <= L / 2, x / 2, (L - x) / 2)
    assert np.isclose(s.castigliano_deflection(x=x, moment=M, dmoment_dP=dM, E=E, I=I), 5 * w * L**4 / (384 * E * I))
    with pytest.raises(ValueError, match="not both"):
        s.castigliano_deflection(0.1, x=x)
    with pytest.raises(ValueError, match="required"):
        s.castigliano_deflection(x=x, moment=M)
    with pytest.raises(ValueError, match="increasing"):
        s.castigliano_deflection(x=x[::-1], moment=M, dmoment_dP=dM, E=E, I=I)


def test_maxwell_count_classifies_trusses():
    assert s.truss_determinacy(3, 3, 3) == 0
    assert s.truss_determinacy(4, 4, 3) == -1
    assert s.truss_determinacy(4, 6, 3) == 1
    with pytest.raises(ValueError):
        s.truss_determinacy(1, 1, 0)


def _warren():
    nodes = np.array([[0, 0], [4, 0], [8, 0], [2, 3], [6, 3]], float)
    members = ((0, 1), (1, 2), (0, 3), (3, 1), (1, 4), (4, 2), (3, 4))
    loads = np.zeros(10)
    loads[3], loads[6] = -10e3, 4e3
    return nodes, members, loads, [0, 1, 5]


def test_method_of_joints_matches_direct_stiffness():
    nodes, members, loads, fixed = _warren()
    joints = s.method_of_joints(nodes, members, loads, fixed)
    stiff = s.PlanarTruss(nodes, members, 200e9, 1e-3).solve(loads, fixed)
    np.testing.assert_allclose(joints.member_forces, stiff.member_forces, rtol=1e-9, atol=1e-6)
    np.testing.assert_allclose(joints.reactions, stiff.reactions, rtol=1e-9, atol=1e-6)
    summary = joints.summary()
    assert summary["tension_members"] + summary["compression_members"] == len(members)
    assert sorted(summary["order"]) == list(range(5))


def test_method_of_joints_rejects_indeterminate_and_stalled_trusses():
    nodes, members, loads, fixed = _warren()
    with pytest.raises(ValueError, match="count"):
        s.method_of_joints(nodes, members, loads, [0, 1, 4, 5])
    with pytest.raises(ValueError, match="two components"):
        s.method_of_joints(nodes, members, loads[:4], fixed)
    # two collinear bars pinned at both ends: every joint has too many or dependent unknowns
    with pytest.raises(ValueError, match="stalls"):
        s.method_of_joints([[0, 0], [1, 0], [2, 0]], ((0, 1), (1, 2)), [0, 0, 0, -1, 0, 0], [0, 1, 4, 5])


def test_truss_conditioning_diagnostics():
    nodes, members, loads, fixed = _warren()
    result = s.PlanarTruss(nodes, members, 200e9, 1e-3).solve(loads, fixed)
    assert 1 < result.summary()["condition_number"] < 1e3
    # Maxwell's count is negative: refused before solving
    with pytest.raises(ValueError, match="Maxwell"):
        s.PlanarTruss(nodes, members[:-1], 200e9, 1e-3).solve(loads, fixed)
    # count zero but a mechanism: an unbraced square on two pins sways
    square = s.PlanarTruss([[0, 0], [1, 0], [1, 1], [0, 1]], ((0, 1), (1, 2), (2, 3), (3, 0)), 200e9, 1e-3)
    with pytest.raises(ValueError, match="singular"):
        square.solve(np.zeros(8), [0, 1, 2, 3])
    # a nearly flat two-bar truss is invertible but ill-conditioned
    flat = s.PlanarTruss([[0, 0], [1, 1e-7], [2, 0]], ((0, 1), (1, 2)), 200e9, 1e-3)
    with pytest.warns(RuntimeWarning, match="ill-conditioned"):
        flat.solve([0, 0, 0, -1.0, 0, 0], [0, 1, 4, 5])


def test_stiffness_utilities_edge_cases():
    u, cond = solve_reduced(np.eye(2), np.zeros(2), np.array([0, 1]), 1e12)
    assert cond == 1.0 and not u.any()
    assert fixed_dof_vector([], 4).size == 0
    with pytest.raises(ValueError):
        fixed_dof_vector([[0, 1]], 4)
    with pytest.raises(ValueError, match="positive length"):
        mesh_geometry([[0, 0], [0, 0]], ((0, 1),))


def test_plastic_section_modulus_and_shape_factor():
    b, h = 0.1, 0.3
    first_moment = 2 * quad(lambda y: b * y, 0, h / 2)[0]
    assert np.isclose(s.plastic_section_modulus(b, h), first_moment)
    assert np.isclose(s.plastic_section_modulus(b, h) / s.rectangular_section_modulus(b, h), 1.5)


def _upper_bound(work_ratio, L):
    return minimize_scalar(work_ratio, bounds=(1e-6 * L, (1 - 1e-6) * L), method="bounded", options={"xatol": 1e-10}).fun


def test_plastic_collapse_loads_minimise_the_mechanism_work():
    Mp, L = 100e3, 6.0
    assert np.isclose(s.plastic_collapse_load(Mp, L, "simple", "uniform"), 8 * Mp / L**2)
    assert np.isclose(s.plastic_collapse_load(Mp, L, "fixed", "uniform"), 16 * Mp / L**2)
    # propped cantilever, uniform load: hinges at the wall and at a from it
    propped = _upper_bound(lambda a: 2 * Mp * (2 / a + 1 / (L - a)) / L, L)
    assert np.isclose(s.plastic_collapse_load(Mp, L, "propped", "uniform"), propped, rtol=1e-8)
    a = 2.0
    assert np.isclose(s.plastic_collapse_load(Mp, L, "simple", position=a), Mp * L / (a * (L - a)))
    assert np.isclose(s.plastic_collapse_load(Mp, L, "propped", position=a), Mp * (2 / a + 1 / (L - a)))
    assert np.isclose(s.plastic_collapse_load(Mp, L), 8 * Mp / L)
    for kwargs in ({"supports": "pinned"}, {"load": "moment"}, {"load": "uniform", "position": 1.0}, {"position": L}):
        with pytest.raises(ValueError):
            s.plastic_collapse_load(Mp, L, **kwargs)


def test_tangent_modulus_load_solves_the_inelastic_buckling_stress():
    E, I, A, L = 200e9, 8e-6, 4e-3, 2.0
    assert s.tangent_modulus_load(E, I, L) == s.euler_buckling_load(E, I, L)
    assert np.isclose(s.tangent_modulus_load(lambda sigma: E, I, L, area=A), s.euler_buckling_load(E, I, L))
    # E_t = E (1 - sigma/sigma_y) gives sigma = c E / (1 + c E / sigma_y), with c = pi^2 r^2 / L^2
    sy = 250e6
    c = np.pi**2 * I / (A * L**2)
    load = s.tangent_modulus_load(lambda sigma: E * max(0.0, 1 - sigma / sy), I, L, area=A)
    assert np.isclose(load, A * c * E / (1 + c * E / sy), rtol=1e-10)
    for Et, kwargs in [
        (lambda sigma: E, {}),
        (lambda sigma: E, {"area": [A, A]}),
        (lambda sigma: E * (1 + sigma / sy), {"area": A}),
        (lambda sigma: 0.0, {"area": A}),
        (lambda sigma: E - 2 * sigma * E, {"area": A}),
    ]:
        with pytest.raises(ValueError):
            s.tangent_modulus_load(Et, I, L, **kwargs)


def test_stress_transformation_rotates_to_principal_axes():
    sx, sy, txy = 80e6, -20e6, 30e6
    _, radius, s1, s2 = s.mohr_circle(sx, sy, txy)
    theta_p = 0.5 * np.arctan2(2 * txy, sx - sy)
    a, b, t = s.stress_transformation(sx, sy, txy, theta_p)
    np.testing.assert_allclose([a, b, t], [s1, s2, 0], atol=1e-6)
    a, b, t = s.stress_transformation(sx, sy, txy, theta_p + np.pi / 4)
    assert np.isclose(abs(t), radius)
    assert np.isclose(a + b, sx + sy)


def test_moment_distribution_matches_the_three_moment_equation():
    w, L, EI = 10e3, 5.0, 1e7
    result = s.moment_distribution([(0, 1, L, EI, w), (1, 2, L, EI, w), (2, 3, L, EI, w)])
    support = s.three_moment_equation([L, L, L], w)
    # sagging moment at the j end of a member is -end_moments[k, 1]
    np.testing.assert_allclose(-result.end_moments[:, 1], support[1:], atol=1e-6)
    np.testing.assert_allclose(result.end_moments[:, 0], support[:-1], atol=1e-6)
    assert result.unbalance[-1] <= 1e-10 * w * L**2 / 12
    # a fixed-fixed member needs no balancing
    fixed = s.moment_distribution([(0, 1, L, EI, w)], fixed_joints=(0, 1))
    assert fixed.iterations == 0
    np.testing.assert_allclose(fixed.end_moments, [[-w * L**2 / 12, w * L**2 / 12]])
    assert fixed.summary()["final_unbalance"] == 0


def test_moment_distribution_errors():
    member = (0, 1, 5.0, 1e7, 1e3)
    for members, kwargs in [
        ([], {}),
        ([(0, 1, 5.0)], {}),
        ([(0, 0, 5.0, 1e7, 1e3)], {}),
        ([member], {"fixed_joints": (2,)}),
    ]:
        with pytest.raises(ValueError):
            s.moment_distribution(members, **kwargs)
    with pytest.raises(RuntimeError, match="converge"):
        s.moment_distribution([member, (1, 2, 5.0, 1e7, 0.0)], fixed_joints=(0,), tol=1e-15, max_iterations=1)


def test_planar_frame_cantilever_and_fixed_beam():
    E, A, I, L, P = 200e9, 1e-2, 1e-4, 3.0, 5e3
    frame = s.PlanarFrame([[0, 0], [L, 0]], ((0, 1),), E, A, I)
    r = frame.solve([0, 0, 0, 2 * P, -P, 0], [0, 1, 2])
    np.testing.assert_allclose(r.displacements[3:], [2 * P * L / (E * A), -P * L**3 / (3 * E * I), -P * L**2 / (2 * E * I)])
    np.testing.assert_allclose(r.reactions[:3], [-2 * P, P, P * L])
    assert np.isclose(r.axial_forces()[0], 2 * P)
    w = 4e3
    fixed = frame.solve(np.zeros(6), [0, 1, 2, 3, 4, 5], member_loads=w)
    np.testing.assert_allclose(fixed.member_end_forces[0], [0, w * L / 2, w * L**2 / 12, 0, w * L / 2, -w * L**2 / 12])
    assert fixed.condition_number == 1.0


def test_planar_frame_portal_sway_and_errors():
    E, A, I, h, L, H = 200e9, 1e-2, 1e-4, 4.0, 6.0, 10e3
    nodes = [[0, 0], [0, h], [L, h], [L, 0]]
    frame = s.PlanarFrame(nodes, ((0, 1), (1, 2), (2, 3)), E, A, I)
    r = frame.solve([0, 0, 0, H, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 2, 9, 10, 11])
    assert np.isclose(r.reactions[0] + r.reactions[9], -H)
    # with an axially rigid, infinitely stiff beam the sway would be H h^3/(24 E I); the real beam is softer
    assert r.displacements[3] > H * h**3 / (24 * E * I)
    summary = r.summary()
    assert summary["max_abs_translation"] >= r.displacements[3]
    with pytest.raises(ValueError, match="three components"):
        frame.solve(np.zeros(6), [0, 1, 2])
    with pytest.raises(ValueError, match="singular"):
        frame.solve(np.zeros(12), [0, 1])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(RuntimeWarning, match="ill-conditioned"):
            frame.solve(np.zeros(12), [0, 1, 2, 9, 10, 11], cond_limit=1.0)
