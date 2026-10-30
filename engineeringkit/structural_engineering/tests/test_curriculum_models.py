"""Analytical checks for the ROADMAP D11 models: Macaulay beams and superposition, built-up and
transformed sections, unsymmetrical bending, Jourawski's shear and the shear centre, Bredt's tubes,
pressure vessels, three-dimensional stress and rosettes, column design curves, the unit-load method,
slope-deflection, cables and three-hinged arches."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from scipy.integrate import quad

from engineeringkit import structural_engineering as s

E, I = 200e9, 8e-6


# --- Macaulay beams and superposition -------------------------------------------------------------


def test_beam_reproduces_the_simply_supported_solutions():
    L = 6.0
    beam = s.Beam(L, E, I, ((0.0, "pin"), (L, "roller")))
    simple = s.SimplySupportedBeam(L, E, I)
    for ours, theirs in [
        (beam.solve(point_loads=[(2.0, 12e3)], n=61), simple.solve_point_load(12e3, 2.0, n=61)),
        (beam.solve(distributed_loads=[(0.0, L, 3e3)], n=61), simple.solve_uniform_load(3e3, n=61)),
    ]:
        for field in ("shear", "moment", "deflection"):
            np.testing.assert_allclose(getattr(ours, field), getattr(theirs, field), atol=1e-9 * np.max(np.abs(getattr(theirs, field))))
        np.testing.assert_allclose(ours.reactions, theirs.reactions)
        assert ours.support_moments == ()


def test_beam_cantilevers_and_indeterminate_beams_match_handbook_formulas():
    L, P, w, C = 3.0, 1e3, 2e3, 5e3
    cantilever = s.Beam(L, E, I, ((0.0, "fixed"),))
    tip = cantilever.solve(point_loads=[(L, P)])
    assert np.isclose(tip.deflection[-1], -P * L**3 / (3 * E * I))
    assert np.isclose(tip.slope[-1], -P * L**2 / (2 * E * I))
    assert tip.reactions == pytest.approx((P,)) and tip.support_moments == pytest.approx((-P * L,))
    assert tip.shear[-1] == pytest.approx(P)  # just left of the load, as in SimplySupportedBeam
    # a clockwise couple at the free end hogs the whole cantilever: the moment steps up by C to zero at the tip
    couple = cantilever.solve(moments=[(L, C)])
    np.testing.assert_allclose(couple.moment[:-1], -C)
    assert np.isclose(couple.deflection[-1], -C * L**2 / (2 * E * I))
    # propped cantilever and fixed-fixed beam under a uniform load
    propped = s.Beam(L, E, I, ((0.0, "fixed"), (L, "roller"))).solve(distributed_loads=[(0.0, L, w)])
    assert propped.reactions == pytest.approx((5 * w * L / 8, 3 * w * L / 8))
    assert propped.support_moments == pytest.approx((-w * L**2 / 8,))
    fixed = s.Beam(L, E, I, ((0.0, "fixed"), (L, "fixed"))).solve(distributed_loads=[(0.0, L, w)], n=401)
    assert fixed.support_moments == pytest.approx((-w * L**2 / 12, w * L**2 / 12))
    assert np.isclose(fixed.deflection[200], -w * L**4 / (384 * E * I))


def test_beam_with_triangular_load_and_overhang():
    L, w0 = 6.0, 4e3
    # load rising from zero at the left: reactions wL/6 and wL/3, peak deflection 0.00652 w L^4 / EI
    r = s.Beam(L, E, I, ((0.0, "pin"), (L, "roller"))).solve(distributed_loads=[(0.0, L, 0.0, w0)], n=20001)
    assert r.reactions == pytest.approx((w0 * L / 6, w0 * L / 3))
    assert np.isclose(r.deflection.min(), -0.0065216 * w0 * L**4 / (E * I), rtol=1e-4)
    assert np.isclose(r.x[np.argmin(r.deflection)], 0.5193 * L, rtol=1e-3)
    # overhang: supports at 0 and 4 m, load at the 6 m tip lifts the left support
    over = s.Beam(L, E, I, ((0.0, "pin"), (4.0, "roller"))).solve(point_loads=[(L, 1e3)], n=601)
    assert over.reactions == pytest.approx((-500.0, 1500.0))
    assert np.isclose(over.moment[np.argmin(np.abs(over.x - 4.0))], -2e3)
    # partial uniform load: equilibrium and curvature M/EI
    part = s.Beam(L, E, I, ((1.0, "pin"), (5.0, "roller"))).solve(distributed_loads=[(2.0, 6.0, w0)], n=6001)
    assert np.isclose(sum(part.reactions), 4 * w0)
    curvature = np.gradient(np.gradient(part.deflection, part.x), part.x)
    smooth = np.min(np.abs(part.x[:, None] - np.array([0.0, 1.0, 2.0, 5.0, 6.0])), axis=1) > 0.005
    np.testing.assert_allclose(curvature[smooth], part.moment[smooth] / (E * I), atol=1e-9)


def test_beam_matches_the_three_moment_equation_on_two_spans():
    L, w = 4.0, 1e3
    r = s.Beam(2 * L, E, I, ((0.0, "pin"), (L, "roller"), (2 * L, "roller"))).solve(distributed_loads=[(0.0, 2 * L, w)], n=401)
    assert np.isclose(r.moment[200], s.three_moment_equation([L, L], w)[1])
    assert r.reactions == pytest.approx((3 * w * L / 8, 10 * w * L / 8, 3 * w * L / 8))


def test_superposition_builds_the_propped_cantilever():
    L, w = 5.0, 2e3
    cantilever = s.Beam(L, E, I, ((0.0, "fixed"),))
    R = 3 * w * L / 8
    total = s.superpose_beam_results(cantilever.solve(distributed_loads=[(0.0, L, w)]), cantilever.solve(point_loads=[(L, -R)]))
    assert np.isclose(total.deflection[-1], 0, atol=1e-15)
    assert total.support_moments == pytest.approx((-w * L**2 / 8,))
    mixed = s.superpose_beam_results(
        s.SimplySupportedBeam(L, E, I).solve_uniform_load(w), s.Beam(L, E, I, ((0, "pin"), (L, "roller"))).solve()
    )
    assert mixed.slope is None and mixed.reactions == pytest.approx((w * L / 2, w * L / 2))


def test_beam_errors():
    L = 4.0
    for supports in [(), ((0.0, "hinge"),), ((0.0,),), ((5.0, "pin"),)]:
        with pytest.raises(ValueError):
            s.Beam(L, E, I, supports)
    beam = s.Beam(L, E, I, ((0.0, "pin"), (L, "roller")))
    for kwargs in [
        {"point_loads": [(5.0, 1.0)]},
        {"moments": [(-1.0, 1.0)]},
        {"distributed_loads": [(0.0, 1.0)]},
        {"distributed_loads": [(2.0, 1.0, 1.0)]},
        {"n": 1},
    ]:
        with pytest.raises(ValueError):
            beam.solve(**kwargs)
    with pytest.raises(ValueError, match="mechanism"):
        s.Beam(L, E, I, ((1.0, "pin"),)).solve(point_loads=[(2.0, 1.0)])
    with pytest.raises(ValueError, match="mechanism"):
        s.Beam(L, E, I, ((1.0, "pin"), (1.0, "roller"))).solve()
    a, b = beam.solve(n=11), beam.solve(n=21)
    with pytest.raises(ValueError, match="stations"):
        s.superpose_beam_results(a, b)
    with pytest.raises(ValueError, match="supports"):
        s.superpose_beam_results(a, s.Beam(L, E, I, ((0.0, "fixed"),)).solve(n=11))
    with pytest.raises(ValueError):
        s.superpose_beam_results()


# --- sections, unsymmetrical bending and shear ------------------------------------------------------


def test_standard_sections_match_hand_formulas():
    d, bf, tf, tw = 0.3, 0.15, 0.012, 0.008
    I_sec = s.i_section(d, bf, tf, tw)
    assert np.isclose(I_sec.Iz, (bf * d**3 - (bf - tw) * (d - 2 * tf) ** 3) / 12)
    assert np.isclose(I_sec.Iy, 2 * tf * bf**3 / 12 + (d - 2 * tf) * tw**3 / 12)
    assert np.isclose(I_sec.y_centroid, d / 2) and I_sec.Iyz == 0
    box = s.box_section(0.2, 0.1, 0.01)
    assert np.isclose(box.Iz, (0.1 * 0.2**3 - 0.08 * 0.18**3) / 12)
    assert np.isclose(box.width_at(0.1), 0.02) and np.isclose(box.width_at(0.005), 0.1)
    tee = s.t_section(0.2, 0.15, 0.02, 0.01)
    A1, A2 = 0.01 * 0.18, 0.15 * 0.02
    assert np.isclose(tee.y_centroid, (A1 * 0.09 + A2 * 0.19) / (A1 + A2))
    channel = s.channel_section(0.2, 0.08, 0.01, 0.008)
    assert np.isclose(channel.z_centroid, (0.2 * 0.008 * 0.004 + 2 * 0.072 * 0.01 * 0.044) / channel.area)
    for args, factory in [((0.02, 0.15, 0.012, 0.008), s.i_section), ((0.3, 0.005, 0.012, 0.008), s.i_section)]:
        with pytest.raises(ValueError):
            factory(*args)
    for bad in [(0.02, 0.15, 0.02, 0.01), (0.2, 0.005, 0.02, 0.01)]:
        with pytest.raises(ValueError):
            s.t_section(*bad)
    with pytest.raises(ValueError):
        s.channel_section(0.2, 0.008, 0.01, 0.008)
    with pytest.raises(ValueError):
        s.angle_section(0.01, 0.1, 0.01)
    with pytest.raises(ValueError):
        s.box_section(0.2, 0.02, 0.01)


def test_composite_section_transformed_and_errors():
    # flitch beam: a steel plate (n = 20) between two timber joists bends like the sum of their EI
    timber, plate, n = (0.05, 0.25), (0.01, 0.25), 20.0
    flitch = s.composite_section([(*timber, -0.035, 0.125), (*plate, 0.0, 0.125, n), (*timber, 0.035, 0.125)])
    EI_sum = 10e9 * 2 * s.rectangular_second_moment(*timber) + 200e9 * s.rectangular_second_moment(*plate)
    assert np.isclose(10e9 * flitch.Iz, EI_sum)
    assert np.isclose(flitch.width_at(0.1), 0.11)
    for rectangles in [[], [(1, 1, 0)], [(1, 1, 0, 0, 0)], [(1, 1, 0, 0, -1)], [(-1, 1, 0, 0)]]:
        with pytest.raises(ValueError):
            s.composite_section(rectangles)


def test_principal_axes_of_an_angle_maximise_the_rotated_moment():
    angle = s.angle_section(0.15, 0.1, 0.01)
    assert angle.Iyz < 0
    I1, I2, theta = angle.principal_moments()
    assert np.isclose(I1 + I2, angle.Iz + angle.Iy)
    phi = np.linspace(-np.pi / 2, np.pi / 2, 20001)
    rotated = (angle.Iz + angle.Iy) / 2 + (angle.Iz - angle.Iy) / 2 * np.cos(2 * phi) - angle.Iyz * np.sin(2 * phi)
    assert np.isclose(rotated.max(), I1) and np.isclose(rotated.min(), I2)
    assert np.isclose(phi[np.argmax(rotated)], theta, atol=1e-3)
    fig, ax = angle.plot()
    assert len(ax.patches) == 2
    plt.close(fig)
    assert set(angle.summary()) == {"area", "y_centroid", "z_centroid", "Iz", "Iy", "Iyz"}


def _grid(section):
    """Gauss-Legendre points (z, y, dA * weight) over the rectangles, exact for the quadratic integrands."""
    u, g = np.polynomial.legendre.leggauss(3)
    points = []
    for b, h, z, y, weight in section.rectangles:
        Z, Y = np.meshgrid(z + b * u / 2, y + h * u / 2)
        points.append((Z.ravel(), Y.ravel(), (weight * b * h / 4 * np.outer(g, g)).ravel()))
    return [np.concatenate(column) for column in zip(*points)]


def test_unsymmetrical_bending_stress_resultants_and_neutral_axis():
    angle = s.angle_section(0.15, 0.1, 0.01)
    My, Mz = 2e3, -5e3
    Z, Y, dA = _grid(angle)
    y, z = Y - angle.y_centroid, Z - angle.z_centroid
    sigma = s.unsymmetrical_bending_stress(My, Mz, y, z, angle.Iy, angle.Iz, angle.Iyz)
    assert np.isclose(np.sum(sigma * dA), 0, atol=1e-9 * np.max(np.abs(sigma)) * angle.area)
    assert np.isclose(-np.sum(sigma * y * dA), Mz, rtol=1e-12)
    assert np.isclose(np.sum(sigma * z * dA), My, rtol=1e-12)
    alpha = s.neutral_axis_angle(My, Mz, angle.Iy, angle.Iz, angle.Iyz)
    assert np.isclose(s.unsymmetrical_bending_stress(My, Mz, np.sin(alpha), np.cos(alpha), angle.Iy, angle.Iz, angle.Iyz), 0, atol=1e-3)
    # principal axes: the flexure formula, and tan(alpha) = (Iz/Iy) tan(theta)
    assert np.isclose(s.unsymmetrical_bending_stress(0, 10, 2, 0, 1, 4), -s.bending_stress(10, 2, 4))
    theta = np.radians([30.0, 150.0, -100.0])
    alpha = s.neutral_axis_angle(np.sin(theta), np.cos(theta), 1.0, 4.0)
    np.testing.assert_allclose(np.tan(alpha), 4 * np.tan(theta))
    assert np.all((alpha > -np.pi / 2) & (alpha <= np.pi / 2))
    with pytest.raises(ValueError, match="both be zero"):
        s.neutral_axis_angle(0, 0, 1, 1)
    with pytest.raises(ValueError, match="positive"):
        s.unsymmetrical_bending_stress(1, 1, 0, 0, 1, 1, 1)


def test_jourawski_shear_stress_integrates_to_the_shear_force():
    b, h, V = 0.1, 0.3, 30e3
    rect = s.composite_section([(b, h, 0, h / 2)])
    assert np.isclose(rect.first_moment(h / 2), b * h**2 / 8)
    assert np.isclose(rect.first_moment(0.0), 0) and rect.first_moment(h) == 0
    peak = s.transverse_shear_stress(V, rect.first_moment(h / 2), rect.Iz, b)
    assert np.isclose(peak, 1.5 * V / (b * h))
    I_sec = s.i_section(0.3, 0.15, 0.012, 0.008)

    def force_per_height(y):
        return s.transverse_shear_stress(V, I_sec.first_moment(y), I_sec.Iz, I_sec.width_at(y)) * I_sec.width_at(y)

    total = quad(force_per_height, 0, 0.3, points=[0.012, 0.288], limit=200)[0]
    assert np.isclose(total, V, rtol=1e-6)
    q = s.shear_flow(V, I_sec.first_moment(0.288), I_sec.Iz)
    assert np.isclose(q, s.transverse_shear_stress(V, I_sec.first_moment(0.288), I_sec.Iz, 0.15) * 0.15)


def test_channel_shear_centre_balances_the_flange_couple():
    h, b, tf, tw = 0.2, 0.08, 0.01, 0.006
    I_thin = tw * h**3 / 12 + 2 * b * tf * (h / 2) ** 2
    flange_force = b**2 * h * tf / (4 * I_thin)  # per unit shear
    assert np.isclose(s.channel_shear_center(h, b, tf, tw), flange_force * h)
    assert np.isclose(s.channel_shear_center(h, b, 0.01, 0.01), 3 * b**2 / (6 * b + h))
    assert s.channel_shear_center(h, 0.0, tf, tw) == 0
    with pytest.raises(ValueError):
        s.channel_shear_center(h, -b, tf, tw)


# --- torsion, vessels, three-dimensional stress and rosettes ---------------------------------------


def test_bredt_reduces_to_the_thin_circular_tube():
    T, r, t, L, G = 2e3, 0.05, 0.002, 1.5, 80e9
    J = 2 * np.pi * r**3 * t
    assert np.isclose(s.thin_walled_torsion_shear_stress(T, np.pi * r**2, t), s.torsion_shear_stress(T, r, J))
    assert np.isclose(s.thin_walled_torsion_angle(T, L, G, np.pi * r**2, 2 * np.pi * r / t), s.torsion_angle(T, L, G, J))
    # a box with walls of two thicknesses: the thin wall carries the higher stress, the shear flow is common
    Am = 0.19 * 0.09
    assert np.isclose(s.thin_walled_torsion_shear_stress(T, Am, 0.005) * 0.005, s.thin_walled_torsion_shear_stress(T, Am, 0.01) * 0.01)


def test_lame_stresses_satisfy_equilibrium_and_boundary_conditions():
    ri, ro, pi, po = 0.1, 0.2, 30e6, 5e6
    r = np.linspace(ri, ro, 2001)
    radial, hoop, axial = s.thick_wall_cylinder_stresses(r, ri, ro, pi, po)
    assert np.isclose(radial[0], -pi) and np.isclose(radial[-1], -po)
    np.testing.assert_allclose(np.gradient(radial, r)[2:-2], ((hoop - radial) / r)[2:-2], rtol=1e-5)
    assert np.isclose(axial[0] * np.pi * (ro**2 - ri**2), pi * np.pi * ri**2 - po * np.pi * ro**2)
    # thin-wall limit: the mean hoop stress is p r / t with r the inner radius
    hoop_thin, long_thin = s.thin_wall_vessel_stresses(1e6, 1.0, 0.005)
    rr = np.linspace(1.0, 1.005, 101)
    mean_hoop = np.trapezoid(s.thick_wall_cylinder_stresses(rr, 1.0, 1.005, 1e6)[1], rr) / 0.005
    assert np.isclose(mean_hoop, hoop_thin) and np.isclose(long_thin, hoop_thin / 2)
    assert s.thin_wall_vessel_stresses(1e6, 1.0, 0.005, "sphere") == pytest.approx((1e8, 1e8))
    with pytest.raises(ValueError):
        s.thin_wall_vessel_stresses(1e6, 1.0, 0.005, "torus")
    with pytest.raises(ValueError):
        s.thick_wall_cylinder_stresses(0.3, ri, ro, pi)
    with pytest.raises(ValueError):
        s.thick_wall_cylinder_stresses(0.1, ri, ri, pi)


def test_principal_stresses_3d_invariants_and_plane_stress():
    sx, sy, sz, txy, tyz, tzx = 50.0, -20.0, 30.0, 25.0, -10.0, 15.0
    s1, s2, s3 = s.principal_stresses_3d(sx, sy, sz, txy, tyz, tzx)
    assert s1 >= s2 >= s3
    assert np.isclose(s1 + s2 + s3, sx + sy + sz)
    det = np.linalg.det([[sx, txy, tzx], [txy, sy, tyz], [tzx, tyz, sz]])
    assert np.isclose(s1 * s2 * s3, det)
    # plane stress with both in-plane principal stresses positive: the absolute maximum shear is out of plane
    _, radius, p1, p2 = s.mohr_circle(80.0, 40.0, 0.0)
    assert s.principal_stresses_3d(80.0, 40.0, 0.0) == pytest.approx((p1, p2, 0.0))
    assert (p1 - 0.0) / 2 > radius
    arrays = s.principal_stresses_3d([80.0, 50.0], [40.0, -20.0], 0.0, [0.0, 25.0])
    assert arrays[0].shape == (2,) and np.isclose(arrays[0][1], s.mohr_circle(50.0, -20.0, 25.0)[2])


def _gauge(ex, ey, gxy, theta):
    return ex * np.cos(theta) ** 2 + ey * np.sin(theta) ** 2 + gxy * np.sin(theta) * np.cos(theta)


def test_strain_rosettes_recover_the_strain_state():
    state = (420e-6, -130e-6, 260e-6)
    for angles in [(0, np.pi / 4, np.pi / 2), (0, np.pi / 3, 2 * np.pi / 3)]:
        readings = [_gauge(*state, a) for a in angles]
        assert s.strain_rosette(readings, angles) == pytest.approx(state)
    series = np.array([[_gauge(k * 1e-4, 0, 0, a) for k in (1, 2)] for a in (0, np.pi / 4, np.pi / 2)])
    ex, _, _ = s.strain_rosette(series)
    np.testing.assert_allclose(ex, [1e-4, 2e-4])
    # Hooke's law round trip: plane stress back to strain
    Em, nu = 70e9, 0.33
    sxx, syy, txy = s.plane_stress_from_strains(*state, Em, nu)
    assert np.isclose((sxx - nu * syy) / Em, state[0]) and np.isclose((syy - nu * sxx) / Em, state[1])
    assert np.isclose(txy * 2 * (1 + nu) / Em, state[2])
    for strains, angles in [([1, 2], (0, 1, 2)), ([1, 2, 3], (0, 1)), ([1, 2, 3], (0, np.pi, np.pi / 2)), (np.ones((3, 2, 2)), (0, 1, 2))]:
        with pytest.raises(ValueError):
            s.strain_rosette(strains, angles)


# --- column design ----------------------------------------------------------------------------------


def test_effective_length_factors():
    K = s.effective_length_factor("fixed-pinned")
    kL = np.pi / K
    assert np.isclose(np.tan(kL), kL) and np.isclose(K, 0.6992, atol=1e-4)
    assert [s.effective_length_factor(e) for e in ("pinned-pinned", "fixed-fixed", "fixed-free", "fixed-guided")] == [1, 0.5, 2, 1]
    with pytest.raises(ValueError):
        s.effective_length_factor("free-free")


def test_secant_formula_limits_and_inverse():
    A, c, r, L = 4e-3, 0.05, 0.04, 3.0
    P_euler = s.euler_buckling_load(E, A * r**2, L)
    assert np.isclose(s.secant_formula_stress(1e5, A, 0.0, c, r, L, E), 1e5 / A)
    # short column: sec -> 1, the stress of eccentric compression P/A (1 + e c / r^2)
    assert np.isclose(s.secant_formula_stress(1e3, A, 0.01, c, r, 1e-3, E), 1e3 / A * (1 + 0.01 * c / r**2))
    near = s.secant_formula_stress(0.999999 * P_euler, A, 0.001, c, r, L, E)
    assert near > 100 * P_euler / A
    with pytest.raises(ValueError, match="Euler"):
        s.secant_formula_stress(P_euler * 1.01, A, 0.01, c, r, L, E)
    P = s.secant_formula_load(250e6, A, 0.01, c, r, L, E)
    assert np.isclose(s.secant_formula_stress(P, A, 0.01, c, r, L, E), 250e6)
    assert s.secant_formula_load(250e6, A, 0.0, c, r, L, E) == pytest.approx(min(250e6 * A, P_euler))
    assert s.secant_formula_load(250e6, A, 0.0, c, r, 30.0, E) == pytest.approx(s.euler_buckling_load(E, A * r**2, 30.0))
    with pytest.raises(ValueError, match="scalars"):
        s.secant_formula_load(250e6, [A, A], 0.01, c, r, L, E)


def test_johnson_parabola_meets_euler_tangentially():
    sy = 250e6
    Cc = np.sqrt(2 * np.pi**2 * E / sy)
    assert s.johnson_critical_stress(sy, E, 0.0) == sy
    h = 1e-6 * Cc
    below, above = s.johnson_critical_stress(sy, E, [Cc - h, Cc + h])
    assert np.isclose(below, sy / 2, rtol=1e-5) and np.isclose(above, sy / 2, rtol=1e-5)
    slope = (above - below) / (2 * h)
    assert np.isclose(slope, -2 * np.pi**2 * E / Cc**3, rtol=1e-4)
    assert np.isclose(s.johnson_critical_stress(sy, E, 200.0), np.pi**2 * E / 200.0**2)


# --- unit-load method and slope-deflection ----------------------------------------------------------


def test_unit_load_method_matches_the_stiffness_solution():
    nodes = np.array([[0, 0], [4, 0], [8, 0], [2, 3], [6, 3]], float)
    members = ((0, 1), (1, 2), (0, 3), (3, 1), (1, 4), (4, 2), (3, 4))
    fixed = [0, 1, 5]
    loads = np.zeros(10)
    loads[3] = -10e3
    A = 1e-3
    truss = s.PlanarTruss(nodes, members, E, A)
    real = truss.solve(loads, fixed)
    unit = np.zeros(10)
    unit[3] = -1.0
    virtual = truss.solve(unit, fixed)
    lengths = [np.hypot(*(nodes[j] - nodes[i])) for i, j in members]
    delta = s.unit_load_truss_deflection(real.member_forces, virtual.member_forces, lengths, A, E)
    assert np.isclose(delta, -real.displacements[3])
    with pytest.raises(ValueError, match="one entry per member"):
        s.unit_load_truss_deflection([1.0], [1.0, 2.0], [1.0], A, E)
    L, P = 3.0, 1e3
    x = np.linspace(0, L, 101)
    assert np.isclose(s.unit_load_beam_deflection(x, -P * (L - x), -(L - x), E, I), P * L**3 / (3 * E * I))


def test_slope_deflection_matches_moment_distribution_and_settlement():
    members = [(0, 1, 5.0, 1e7, 10e3), (1, 2, 4.0, 2e7, 0.0), (2, 3, 6.0, 1e7, 5e3)]
    ours = s.slope_deflection(members, fixed_joints=(0,))
    cross = s.moment_distribution(members, fixed_joints=(0,))
    np.testing.assert_allclose(ours.end_moments, cross.end_moments, atol=1e-6)
    assert ours.rotations[0] == 0
    assert set(ours.summary()) == {"max_abs_end_moment", "max_abs_rotation"}
    # fixed-fixed beam, right end settles by Delta: end moments -6 EI Delta / L^2 (clockwise)
    EI, L, delta = 1e7, 5.0, 0.01
    settled = s.slope_deflection([(0, 1, L, EI, 0.0)], fixed_joints=(0, 1), chord_rotations=[delta / L])
    np.testing.assert_allclose(settled.end_moments, [[-6 * EI * delta / L**2] * 2])
    # pinned far end: the moment drops to 3 EI Delta / L^2 at the fixed end
    pinned = s.slope_deflection([(0, 1, L, EI, 0.0)], fixed_joints=(0,), chord_rotations=[delta / L])
    assert np.isclose(pinned.end_moments[0, 0], -3 * EI * delta / L**2) and np.isclose(pinned.end_moments[0, 1], 0)
    with pytest.raises(ValueError, match="one entry per member"):
        s.slope_deflection(members, chord_rotations=[0.0])


# --- cables and arches ------------------------------------------------------------------------------


def test_parabolic_and_catenary_cables():
    w, L, f = 2e3, 100.0, 12.0
    H, T, length = s.parabolic_cable(w, L, f)
    assert np.isclose(H, w * L**2 / (8 * f)) and np.isclose(T, np.hypot(H, w * L / 2))
    arc = quad(lambda x: np.hypot(1, 4 * f * (L - 2 * x) / L**2), 0, L)[0]
    assert np.isclose(length, arc)
    Hc, Tc, sc = s.catenary_cable(w, L, f)
    a = Hc / w
    assert np.isclose(a * (np.cosh(L / (2 * a)) - 1), f)
    assert np.isclose(Tc, np.hypot(Hc, w * sc / 2))  # vertical support force carries half the cable weight
    # a shallow cable: the catenary approaches the parabola
    shallow = np.array(s.catenary_cable(w, L, 0.5)) / np.array(s.parabolic_cable(w, L, 0.5))
    np.testing.assert_allclose(shallow, 1, rtol=1e-4)
    with pytest.raises(ValueError, match="scalars"):
        s.catenary_cable([w, w], L, f)


def test_three_hinged_arch_is_funicular_for_a_uniform_load():
    L, f, w, P = 30.0, 6.0, 4e3, 50e3
    uniform = s.three_hinged_arch(L, f, w=w)
    assert np.isclose(uniform.thrust, w * L**2 / (8 * f))
    np.testing.assert_allclose(uniform.moment, 0, atol=1e-6 * w * L**2)
    np.testing.assert_allclose(-uniform.axial_force, np.hypot(uniform.thrust, w * (L / 2 - uniform.x)), rtol=1e-12)
    crown = s.three_hinged_arch(L, f, point_loads=[(L / 2, P)], n=301)
    assert np.isclose(crown.thrust, P * L / (4 * f))
    off = s.three_hinged_arch(L, f, point_loads=[(7.5, P)], n=301)
    np.testing.assert_allclose(off.moment[[0, 150, -1]], 0, atol=1e-6)
    assert off.reactions == pytest.approx((0.75 * P, 0.25 * P))
    summary = off.summary()
    assert summary["max_abs_moment"] > 0 and summary["max_compression"] > 0
    fig, axes = off.plot()
    assert len(axes) == 2
    plt.close(fig)
    with pytest.raises(ValueError):
        s.three_hinged_arch(L, f, point_loads=[(31.0, P)])
