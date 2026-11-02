"""Worked examples from Anderson's textbooks, checked against the answers the books print (ROADMAP B4).

J. D. Anderson, *Fundamentals of Aerodynamics*, 3rd ed. (McGraw-Hill, 2001), and *Modern
Compressible Flow*, 3rd ed. (McGraw-Hill, 2003). Anderson reads ratios from his appendix tables
(four significant figures) or charts, so each tolerance is half a unit in the last printed digit,
widened only where the table rounding is larger; the differences of method are in the guide.
"""

import numpy as np
import pytest
from scipy.optimize import brentq

from engineeringkit import aerodynamics as a

R = 287.0
ATM = 101325.0


def test_fundamentals_example_4_4_thin_airfoil_lift_at_five_degrees():
    # Fundamentals, 3rd ed., Example 4.4: symmetric thin airfoil at 5 deg, c_l = 0.548.
    assert a.thin_airfoil_cl(np.radians(5.0)) == pytest.approx(0.548, abs=5e-4)


def test_fundamentals_example_4_8_turbulent_flat_plate_skin_friction():
    # Fundamentals, 3rd ed., Example 4.8: c = 1.5 m, Re_c = 3.1e6, fully turbulent:
    # C_f = 0.00372 per side and 0.00745 for both sides.
    cf = a.skin_friction_turbulent(3.1e6)
    assert cf == pytest.approx(0.00372, abs=5e-6)
    assert 2 * cf == pytest.approx(0.00745, abs=5e-6)


def test_fundamentals_example_5_2_induced_drag_at_a_larger_aspect_ratio():
    # Fundamentals, 3rd ed., Example 5.2: a wing of AR 6 with delta = tau = 0.055 has
    # C_D,i = 0.01 at 3.4 deg (zero-lift angle -2 deg). The same section on a wing of AR 10
    # at the same angle has C_D,i = 0.00741.
    delta, alpha, alpha_l0 = 0.055, np.radians(3.4), np.radians(-2.0)
    e = 1 / (1 + delta)
    CL6 = np.sqrt(0.01 * np.pi * e * 6)
    a6 = CL6 / (alpha - alpha_l0)
    a0 = a6 / (1 - a6 / (np.pi * e * 6))
    assert a.finite_wing_lift_slope(6, e, a0) == pytest.approx(a6)
    CL10 = a.finite_wing_lift_slope(10, e, a0) * (alpha - alpha_l0)
    assert a.induced_drag_coefficient(CL10, 10, e) == pytest.approx(0.00741, abs=5e-6)


def test_fundamentals_example_8_5_stagnation_and_sonic_reference_state():
    # Fundamentals, 3rd ed., Example 8.5: M = 3.5, p = 0.3 atm, T = 180 K gives
    # p0 = 22.88 atm, T0 = 621 K, T* = 517.5 K, a* = 456 m/s and M* = 2.064.
    M, p, T = 3.5, 0.3 * ATM, 180.0
    assert a.stagnation_pressure(p, M) / ATM == pytest.approx(22.88, abs=5e-3)
    T0 = a.stagnation_temperature(T, M)
    assert T0 == pytest.approx(621.0, abs=0.5)
    T_star = T0 * 2 / 2.4
    assert T_star == pytest.approx(517.5, abs=0.05)
    a_star = a.speed_of_sound(1.4, R, T_star)
    assert a_star == pytest.approx(456.0, abs=0.5)
    assert M * a.speed_of_sound(1.4, R, T) / a_star == pytest.approx(2.064, abs=5e-4)


def test_fundamentals_example_9_2_oblique_shock_on_a_twenty_degree_wedge():
    # Fundamentals, 3rd ed., Example 9.2: M1 = 2, theta = 20 deg, p1 = 1 atm, T1 = 288 K.
    # Anderson reads beta = 53.4 deg from the chart, rounds M_n1 = 2 sin(beta) = 1.606 to the
    # table entry 1.6 and prints p2 = 2.82 atm and T2 = 399.7 K. The exact M_n1 gives
    # 2.843 atm and 401.0 K: a difference of method, not of theory.
    shock = a.oblique_shock(2.0, np.radians(20.0))
    assert np.degrees(shock.beta) == pytest.approx(53.4, abs=0.05)
    assert shock.mach2 == pytest.approx(1.21, abs=5e-3)
    assert 2 * np.sin(shock.beta) == pytest.approx(1.6, abs=0.01)
    table = a.normal_shock_pressure_ratio(1.6)
    assert table == pytest.approx(2.82, abs=5e-3)
    assert table * (2 + 0.4 * 1.6**2) / (2.4 * 1.6**2) * 288 == pytest.approx(399.7, abs=0.05)
    assert shock.pressure_ratio == pytest.approx(2.82, rel=0.01)


def test_fundamentals_example_9_8_prandtl_meyer_expansion_around_a_corner():
    # Fundamentals, 3rd ed., Example 9.8: M1 = 1.5, p1 = 1 atm, T1 = 288 K, 15 deg corner.
    # nu1 = 11.91 deg, so nu2 = 26.91 deg; Anderson takes the nearest table entry M2 = 2.0 and
    # prints p2 = 0.4692 atm, T2 = 232 K, a forward Mach line at 41.81 deg and a rearward one at
    # 15 deg. The root of nu(M2) = 26.91 deg is 2.019, which gives 0.455 atm and 230.0 K.
    assert np.degrees(a.prandtl_meyer_angle(1.5)) == pytest.approx(11.91, abs=5e-3)
    M2 = a.prandtl_meyer_expansion(1.5, np.radians(15.0))
    assert M2 == pytest.approx(2.0, abs=0.025)
    p2 = a.stagnation_pressure(1.0, 1.5) / a.stagnation_pressure(1.0, 2.0)
    T2 = a.stagnation_temperature(288.0, 1.5) / a.stagnation_temperature(1.0, 2.0)
    assert p2 == pytest.approx(0.4692, abs=5e-4)
    assert T2 == pytest.approx(232.0, abs=0.5)
    assert np.degrees(a.mach_angle(1.5)) == pytest.approx(41.81, abs=5e-3)
    assert np.degrees(a.mach_angle(2.0)) - 15.0 == pytest.approx(15.0, abs=5e-3)


def test_compressible_flow_example_3_4_normal_shock_in_a_wind_tunnel():
    # Modern Compressible Flow, 3rd ed., Example 3.4: M1 = 3, p1 = 0.5 atm, T1 = 200 K gives
    # M2 = 0.4752, p2 = 5.165 atm, T2 = 535.8 K and u2 = 220.5 m/s. p2 uses the table
    # ratio 10.33 for the exact 31/3, hence the wider tolerance on it.
    M2 = a.normal_shock_mach2(3.0)
    ratio = a.normal_shock_pressure_ratio(3.0)
    T2 = 200.0 * ratio * (2 + 0.4 * 9) / (2.4 * 9)
    assert M2 == pytest.approx(0.4752, abs=5e-5)
    assert 0.5 * ratio == pytest.approx(5.165, abs=2e-3)
    assert T2 == pytest.approx(535.8, abs=0.05)
    assert M2 * a.speed_of_sound(1.4, R, T2) == pytest.approx(220.5, abs=0.05)


def test_compressible_flow_example_3_10_heat_that_chokes_a_supersonic_duct():
    # Modern Compressible Flow, 3rd ed., Example 3.10: M1 = 3 and T01 = 840 K choke after
    # q* = 4.46e5 J/kg (cp = 1004.5). The table ratio T0/T0* = 0.6540 gives 4.464e5; the exact
    # 0.65398 gives 4.465e5, which rounds the other way.
    q = a.rayleigh_choking_heat(3.0, 840.0)
    assert q == pytest.approx(4.465e5, abs=50)
    assert 1004.5 * 840.0 * (1 / 0.6540 - 1) == pytest.approx(4.46e5, abs=500)


def test_compressible_flow_examples_3_11_and_3_13_pipe_flow_with_friction():
    # Modern Compressible Flow, 3rd ed., Example 3.11: D = 0.15 m, L = 30 m, M1 = 0.3,
    # p1 = 1 atm, T1 = 273 K and Fanning f = 0.005 (Darcy 0.02) give M2 = 0.475,
    # p2 = 0.624 atm, T2 = 266 K and p02 = 0.728 atm. Example 3.13: the M1 = 3 duct of
    # Example 3.12 (D = 0.4 ft) chokes after 10.44 ft.
    D, f = 0.15, 0.02
    remaining = a.fanno_choking_length(0.3, D, f) - 30.0
    M2 = brentq(lambda M: a.fanno_choking_length(M, D, f) - remaining, 0.3, 0.99)
    inlet, outlet = a.fanno_flow(0.3), a.fanno_flow(M2)
    assert M2 == pytest.approx(0.475, abs=1e-3)
    assert outlet.pressure_ratio / inlet.pressure_ratio == pytest.approx(0.624, abs=1e-3)
    assert 273 * outlet.temperature_ratio / inlet.temperature_ratio == pytest.approx(266, abs=0.5)
    p01 = a.stagnation_pressure(1.0, 0.3)
    assert p01 * outlet.stagnation_pressure_ratio / inlet.stagnation_pressure_ratio == pytest.approx(0.728, abs=1e-3)
    assert a.fanno_choking_length(3.0, 0.4, f) == pytest.approx(10.44, abs=5e-3)
