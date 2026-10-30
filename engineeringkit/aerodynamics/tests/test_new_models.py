"""Analytical checks for circulation lift, NACA sections, Breguet range, the Mach cone, the ISA and the D1 curriculum models."""

import numpy as np
import pytest

from engineeringkit import aerodynamics as a


def test_kutta_joukowski_matches_coefficient_lift_for_thin_airfoil_circulation():
    # Thin-airfoil circulation for a flat plate: Gamma = pi c V alpha.
    rho, V, c, alpha = 1.225, 40.0, 1.5, 0.05
    gamma = np.pi * c * V * alpha
    lift_per_span = a.kutta_joukowski_lift(rho, V, gamma)
    assert lift_per_span == pytest.approx(a.dynamic_pressure(rho, V) * c * a.thin_airfoil_cl(alpha))


def test_symmetric_section_has_zero_lift_at_zero_incidence():
    assert a.naca4_zero_lift_angle(0.0, 0.4) == pytest.approx(0.0, abs=1e-15)


def test_naca2412_zero_lift_angle_matches_tabulated_thin_airfoil_value():
    # Anderson, Fundamentals of Aerodynamics, Example 4.6: -2.077 deg.
    assert np.degrees(a.naca4_zero_lift_angle(0.02, 0.4)) == pytest.approx(-2.077, abs=2e-3)


def test_naca4_geometry_recovers_designation():
    foil = a.naca4_airfoil("2412", n=401)
    summary = foil.summary()
    assert summary["max_thickness"] == pytest.approx(0.12, rel=1e-3)
    assert summary["max_thickness_at"] == pytest.approx(0.30, abs=0.01)
    assert summary["max_camber"] == pytest.approx(0.02, rel=1e-4)
    assert summary["max_camber_at"] == pytest.approx(0.40, abs=0.01)
    assert foil.upper[0] == pytest.approx([0, 0])


def test_naca4_closed_trailing_edge_and_symmetric_section():
    closed = a.naca4_airfoil("0012", closed_trailing_edge=True)
    assert closed.upper[-1, 1] == pytest.approx(0.0, abs=1e-5)
    assert np.allclose(closed.camber, 0)
    assert np.allclose(closed.upper[:, 1], -closed.lower[:, 1])


@pytest.mark.parametrize("designation", ["241", "24x2", 2412, "0000", "2012"])
def test_naca4_rejects_bad_designations(designation):
    with pytest.raises(ValueError):
        a.naca4_airfoil(designation)


def test_naca4_rejects_too_few_stations():
    with pytest.raises(ValueError):
        a.naca4_airfoil("0012", n=2)


def test_breguet_range_is_linear_in_lift_to_drag_and_zero_without_fuel():
    r1 = a.breguet_range(230, 1.6e-5, 10, 7e5, 5e5)
    assert a.breguet_range(230, 1.6e-5, 20, 7e5, 5e5) == pytest.approx(2 * r1)
    assert a.breguet_range(230, 1.6e-5, 10, 5e5, 5e5) == 0
    with pytest.raises(ValueError):
        a.breguet_range(230, 1.6e-5, 10, 4e5, 5e5)


def test_mach_angle_limits():
    assert a.mach_angle(1.0) == pytest.approx(np.pi / 2)
    assert np.sin(a.mach_angle(np.array([2.0, 4.0]))) == pytest.approx([0.5, 0.25])
    with pytest.raises(ValueError):
        a.mach_angle(0.9)


def test_isa_sea_level_and_tropopause():
    sea = a.isa_atmosphere(0.0)
    assert sea.density == pytest.approx(1.225, rel=1e-4)
    assert sea.speed_of_sound == pytest.approx(340.294, rel=1e-5)
    trop = a.isa_atmosphere(11000.0)
    assert trop.pressure == pytest.approx(22632.06, rel=1e-5)
    top = a.isa_atmosphere(20000.0)
    assert top.pressure == pytest.approx(5474.89, rel=1e-4)
    assert top.temperature == pytest.approx(216.65)


def test_isa_is_hydrostatic():
    h = np.linspace(0, 20000, 2001)
    state = a.isa_atmosphere(h)
    dpdh = np.gradient(state.pressure, h)
    assert dpdh[1:-1] == pytest.approx(-state.density[1:-1] * 9.80665, rel=1e-4)
    summary = state.summary()
    assert summary["max_altitude_m"] == 20000
    assert 0 < summary["min_density_ratio"] < summary["min_pressure_ratio"] * 2


def test_isa_rejects_out_of_range_altitude():
    with pytest.raises(ValueError):
        a.isa_atmosphere(25000.0)


def test_compressible_relations_accept_arrays():
    M = np.array([1.0, 2.0, 3.0])
    assert a.normal_shock_mach2(M)[0] == pytest.approx(1.0)
    assert a.normal_shock_pressure_ratio(M)[1] == pytest.approx(4.5)
    assert a.area_mach_ratio(M)[0] == pytest.approx(1.0)
    assert a.prandtl_glauert_cp(-1, np.array([0.0, 0.6]))[1] == pytest.approx(-1.25)
    assert a.induced_drag_coefficient(0.5, 8, np.array([0.5, 1.0]))[1] == pytest.approx(0.25 / (8 * np.pi))


def test_isentropic_factor_and_gamma_checks():
    assert a.isentropic_factor(1.0, 1.4) == pytest.approx(1.2)
    with pytest.raises(ValueError):
        a.isentropic_factor(1.0, 1.0)
    with pytest.raises(ValueError):
        a.thin_airfoil_cl(np.nan)


# Supersonic flow: Anderson, Modern Compressible Flow, 3rd ed., example 4.1 and appendices B and C.


def test_oblique_shock_matches_anderson_mach_two_twenty_degree_wedge():
    shock = a.oblique_shock(2.0, np.radians(20.0))
    assert np.isclose(np.degrees(shock.beta), 53.42, atol=0.01)
    assert np.isclose(shock.mach2, 1.21, atol=0.005)
    assert np.isclose(shock.pressure_ratio, a.normal_shock_pressure_ratio(2.0 * np.sin(shock.beta)))
    assert 0 < shock.stagnation_pressure_ratio < 1
    strong = a.oblique_shock(2.0, np.radians(20.0), strong=True)
    assert strong.beta > shock.beta and strong.mach2 < 1


def test_oblique_shock_branches_invert_the_theta_beta_mach_relation():
    theta = np.radians(np.linspace(0.5, 30.0, 40))
    for strong in (False, True):
        beta = a.oblique_shock_angle(4.0, theta, strong=strong)
        assert np.allclose(a.oblique_shock_deflection(4.0, beta), theta, atol=1e-10)
    theta_max, beta_max = a.oblique_shock_max_deflection(2.0)
    assert np.isclose(np.degrees(theta_max), 22.97, atol=0.01)
    assert np.isclose(a.oblique_shock_angle(2.0, theta_max), beta_max)


def test_mach_wave_and_normal_shock_are_zero_deflection_limits():
    M = np.array([1.5, 2.0, 5.0])
    assert np.allclose(a.oblique_shock_deflection(M, a.mach_angle(M)), 0, atol=1e-12)
    assert np.allclose(a.oblique_shock_deflection(M, np.pi / 2), 0, atol=1e-12)


@pytest.mark.parametrize(
    "fn,args",
    [
        (a.oblique_shock_deflection, (2.0, 0.1)),
        (a.oblique_shock_angle, (2.0, np.radians(25.0))),
        (a.oblique_shock_angle, (1.0, 0.1)),
        (a.oblique_shock_angle, (2.0, 0.0)),
        (a.inverse_prandtl_meyer, (np.radians(131.0),)),
        (a.prandtl_meyer_expansion, (2.0, -0.1)),
    ],
)
def test_supersonic_input_contracts(fn, args):
    with pytest.raises(ValueError):
        fn(*args)


def test_prandtl_meyer_tabulated_values_and_inverse():
    assert np.allclose(np.degrees(a.prandtl_meyer_angle([1.0, 1.5, 2.0, 3.0])), [0.0, 11.905, 26.380, 49.757], atol=1e-3)
    assert np.isclose(a.prandtl_meyer_expansion(1.5, np.radians(20.0)), 2.207, atol=1e-3)
    assert a.inverse_prandtl_meyer(0.0) == 1.0
    M = np.array([1.2, 2.0, 7.5, 40.0])
    assert np.allclose(a.inverse_prandtl_meyer(a.prandtl_meyer_angle(M)), M, rtol=1e-9)


# Duct flow: Anderson, Modern Compressible Flow, 3rd ed., appendices C and D.


def test_fanno_and_rayleigh_tables_at_mach_two():
    fanno = a.fanno_flow(2.0)
    assert np.isclose(fanno.friction_parameter / 4, 0.30500, atol=1e-5)
    assert np.allclose(
        [fanno.temperature_ratio, fanno.pressure_ratio, fanno.stagnation_pressure_ratio], [0.66667, 0.40825, 1.6875], atol=1e-5
    )
    rayleigh = a.rayleigh_flow(2.0)
    assert np.allclose(
        [rayleigh.pressure_ratio, rayleigh.temperature_ratio, rayleigh.stagnation_temperature_ratio, rayleigh.stagnation_pressure_ratio],
        [0.36364, 0.52893, 0.79339, 1.5031],
        atol=1e-4,
    )


def test_duct_flows_are_sonic_reference_at_mach_one():
    for state in (a.fanno_flow(1.0), a.rayleigh_flow(1.0)):
        for name in ("temperature_ratio", "pressure_ratio", "density_ratio", "stagnation_pressure_ratio", "stagnation_temperature_ratio"):
            assert np.isclose(getattr(state, name), 1.0)
        assert np.isclose(state.friction_parameter, 0.0)


def test_choking_length_and_heat():
    assert np.isclose(a.fanno_choking_length(0.5, 0.05, 0.02), a.fanno_flow(0.5).friction_parameter * 0.05 / 0.02)
    assert np.isclose(a.rayleigh_choking_heat(1.0, 300.0), 0.0, atol=1e-9)
    q = a.rayleigh_choking_heat(0.2, 300.0)
    cp = 1.4 * a.R_AIR / 0.4
    assert np.isclose(a.rayleigh_flow(0.2).stagnation_temperature_ratio * (300.0 + q / cp), 300.0)


# Wings and sections: Anderson, Fundamentals of Aerodynamics, 6th ed., sec. 5.3 and example 4.6.


def test_lifting_line_elliptic_wing_is_exact():
    b, c0, alpha = 10.0, 1.25, 0.1
    wing = a.lifting_line(b, lambda y: c0 * np.sqrt(1 - (2 * y / b) ** 2), alpha, n_terms=21)
    assert np.isclose(wing.aspect_ratio, b / (np.pi * c0 / 4))
    assert np.isclose(wing.span_efficiency, 1.0)
    assert np.isclose(wing.CL, a.finite_wing_cl(alpha, wing.aspect_ratio))
    assert np.isclose(wing.CDi, a.induced_drag_coefficient(wing.CL, wing.aspect_ratio))
    assert np.allclose(wing.cl, wing.cl[0])
    assert np.allclose(wing.induced_angle, wing.CL / (np.pi * wing.aspect_ratio))


def test_lifting_line_rectangular_wing_converges_below_elliptic():
    coarse, fine = (a.lifting_line(8.0, 1.0, 0.1, n_terms=n) for n in (20, 80))
    assert np.isclose(coarse.CL, fine.CL, rtol=2e-3)
    assert 0.9 < fine.span_efficiency < 1.0
    assert fine.CL < a.finite_wing_cl(0.1, 8.0)


def test_lifting_line_twist_and_zero_lift():
    antisymmetric = a.lifting_line(8.0, 1.0, lambda y: 0.01 * y)
    assert np.isclose(antisymmetric.CL, 0.0, atol=1e-12) and antisymmetric.CDi > 0
    assert np.isclose(antisymmetric.span_efficiency, 0.0)
    assert np.isnan(a.lifting_line(8.0, 1.0, 0.0).span_efficiency)
    with pytest.raises(ValueError):
        a.lifting_line(8.0, lambda y: 1 - abs(y), 0.1)
    with pytest.raises(ValueError):
        a.lifting_line(8.0, 1.0, 0.1, n_terms=0)


def test_naca2412_moment_and_aerodynamic_center():
    assert np.isclose(a.naca4_moment_coefficient(0.02, 0.4), -0.053, atol=5e-4)
    assert a.naca4_moment_coefficient(0.0, 0.4) == 0.0
    assert a.aerodynamic_center(2 * np.pi, 0.0) == 0.25
    assert np.isclose(a.aerodynamic_center(6.0, -0.06), 0.26)


def test_panel_method_circle_is_exact_and_sections_follow_thin_airfoil_theory():
    t = np.linspace(0, -2 * np.pi, 65)
    circle = a.panel_method(np.column_stack((np.cos(t), np.sin(t))), 0.0)
    assert np.allclose(circle.cp, 1 - 4 * np.sin(np.arctan2(circle.y, circle.x)) ** 2, atol=1e-12)
    assert np.isclose(circle.cl, 0.0, atol=1e-12)
    symmetric = a.panel_method(a.naca4_airfoil("0012", 81), np.radians(5.0))
    assert np.isclose(symmetric.cl, 2 * np.pi * np.radians(5.0) * (1 + 0.77 * 0.12), rtol=0.01)
    assert np.isclose(symmetric.cl, symmetric.cl_pressure, rtol=0.01)
    cambered = a.panel_method(a.naca4_airfoil("2412", 81), 0.0)
    assert np.isclose(cambered.cm_quarter_chord, a.naca4_moment_coefficient(0.02, 0.4), atol=0.003)
    assert cambered.cl > 0
    assert np.isclose(a.panel_method(a.naca4_airfoil("0012", 81), 0.0).cl, 0.0, atol=1e-12)


@pytest.mark.parametrize(
    "geometry",
    [
        np.zeros((3, 2)),
        np.zeros((5, 3)),
        np.column_stack((np.cos(np.linspace(0, 2 * np.pi, 9)), np.sin(np.linspace(0, 2 * np.pi, 9)))),
        [[1, 0], [0, -1], [0, -1], [-1, 0], [0, 1], [1, 0]],
    ],
)
def test_panel_method_rejects_bad_geometry(geometry):
    with pytest.raises(ValueError):
        a.panel_method(geometry, 0.0)


# Stability and performance: Anderson, Introduction to Flight, 8th ed., chs. 6 and 7.


def test_neutral_point_static_margin_and_tail_volume():
    V_H = a.tail_volume_ratio(3.0, 5.0, 16.0, 1.5)
    assert np.isclose(V_H, 0.625)
    h_n = a.neutral_point(0.24, 4.8, 4.0, V_H, a.downwash_gradient(4.8, 8.0))
    assert np.isclose(h_n, 0.24 + V_H * 4.0 / 4.8 * (1 - 2 * 4.8 / (np.pi * 8.0)))
    assert np.isclose(a.static_margin(h_n, 0.3), h_n - 0.3)
    assert a.pitching_moment_slope(0.3, 0.24, 4.8, 4.0, V_H) < 0 < a.pitching_moment_slope(0.9, 0.24, 4.8, 4.0, V_H)
    assert np.isclose(a.neutral_point(0.25, 5.0, 4.0, 0.5, eta=0.0 + 1e-12), 0.25, atol=1e-12)
    with pytest.raises(ValueError):
        a.neutral_point(0.25, 5.0, 4.0, 0.5, deps_dalpha=1.0)


def test_climb_turn_and_takeoff_relations():
    T, D, V, W = 3000.0, 1500.0, 60.0, 15000.0
    assert np.isclose(a.rate_of_climb(T, D, V, W), V * np.sin(a.climb_angle(T, D, W)))
    n = a.load_factor(np.radians(60.0))
    assert np.isclose(n, 2.0)
    assert np.isclose(a.turn_rate(V, n) * a.turn_radius(V, n), V)
    with pytest.raises(ValueError):
        a.climb_angle(20000.0, 0.0, W)
    with pytest.raises(ValueError):
        a.turn_radius(V, 1.0)
    with pytest.raises(ValueError):
        a.takeoff_ground_roll(W, 1.225, 16.0, 1.6, 400.0, 400.0)
    with pytest.raises(ValueError):
        a.breguet_range_propeller(0.8, 7.6e-8, 13.0, 1.0e4, 1.1e4)
    assert a.breguet_range_propeller(0.8, 7.6e-8, 13.0, 1.0e4, 1.0e4) == 0.0
