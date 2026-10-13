"""Analytical checks for circulation lift, NACA sections, Breguet range, the Mach cone and the ISA."""

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
