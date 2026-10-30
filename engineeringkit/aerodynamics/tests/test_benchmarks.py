import numpy as np
import pytest
from scipy.integrate import quad

from engineeringkit import aerodynamics as a


def test_mach_two_normal_shock_mass_momentum_energy():
    M, gamma, R, T, p = 2.0, 1.4, 287.0, 300.0, 1e5
    M2 = a.normal_shock_mach2(M, gamma)
    ratio = a.normal_shock_pressure_ratio(M, gamma)
    assert np.isclose(M2, np.sqrt(1 / 3))
    assert np.isclose(ratio, 4.5)
    density_ratio = (gamma + 1) * M * M / ((gamma - 1) * M * M + 2)
    T2 = T * ratio / density_ratio
    p2 = p * ratio
    rho1, rho2 = p / (R * T), p2 / (R * T2)
    u1 = M * a.speed_of_sound(gamma, R, T)
    u2 = M2 * a.speed_of_sound(gamma, R, T2)
    assert np.isclose(rho1 * u1, rho2 * u2)
    assert np.isclose(p + rho1 * u1 * u1, p2 + rho2 * u2 * u2)
    cp = gamma * R / (gamma - 1)
    assert np.isclose(cp * T + u1 * u1 / 2, cp * T2 + u2 * u2 / 2)


def test_sonic_and_infinite_wing_limits():
    assert a.area_mach_ratio(1) == 1
    assert a.normal_shock_mach2(1) == 1
    assert a.normal_shock_pressure_ratio(1) == 1
    assert np.isclose(a.finite_wing_lift_slope(1e10), 2 * np.pi)
    assert a.area_mach_ratio(0.5) > 1 and a.area_mach_ratio(2) > 1


def test_best_glide_is_stationary():
    cl = a.best_glide_cl(0.02, 8, 0.8)
    optimum = a.glide_ratio(cl, a.drag_polar(cl, 0.02, 8, 0.8))
    assert all(optimum > a.glide_ratio(x, a.drag_polar(x, 0.02, 8, 0.8)) for x in (cl * 0.9, cl * 1.1))


@pytest.mark.parametrize("M", [-1, 1, 1.1, np.nan])
def test_prandtl_glauert_invalid(M):
    with pytest.raises(ValueError):
        a.prandtl_glauert_cp(-0.5, M)


@pytest.mark.parametrize("M", [0, 0.5])
def test_subsonic_is_not_a_normal_shock(M):
    with pytest.raises(ValueError):
        a.normal_shock_mach2(M)


def test_oblique_shock_conserves_mass_momentum_and_energy():
    gamma, R, p1, T1 = 1.4, 287.05, 1e5, 250.0
    shock = a.oblique_shock(3.0, np.radians(15.0), gamma)
    rho1 = p1 / (R * T1)
    V1 = 3.0 * a.speed_of_sound(gamma, R, T1)
    p2, rho2, T2 = p1 * shock.pressure_ratio, rho1 * shock.density_ratio, T1 * shock.temperature_ratio
    V2 = shock.mach2 * a.speed_of_sound(gamma, R, T2)
    b, t = shock.beta, shock.beta - shock.theta
    u1, w1, u2, w2 = V1 * np.sin(b), V1 * np.cos(b), V2 * np.sin(t), V2 * np.cos(t)
    cp = gamma * R / (gamma - 1)
    assert np.isclose(rho1 * u1, rho2 * u2)
    assert np.isclose(p1 + rho1 * u1**2, p2 + rho2 * u2**2)
    assert np.isclose(w1, w2)
    assert np.isclose(cp * T1 + V1**2 / 2, cp * T2 + V2**2 / 2)


def test_prandtl_meyer_function_is_the_integral_of_its_ode():
    gamma = 1.3
    for M in (1.4, 2.5, 6.0):
        integral = quad(lambda m: np.sqrt(m**2 - 1) / (1 + (gamma - 1) / 2 * m**2) / m, 1, M)[0]
        assert np.isclose(a.prandtl_meyer_angle(M, gamma), integral)
    limit = np.pi / 2 * (np.sqrt(2.4 / 0.4) - 1)
    assert np.isclose(a.prandtl_meyer_angle(1e8), limit, rtol=1e-6)


def test_fanno_friction_parameter_integrates_the_friction_ode():
    gamma = 1.4

    def rate(m):
        return 2 * (1 - m**2) / (gamma * m**3 * (1 + (gamma - 1) / 2 * m**2))

    for M in (0.3, 0.8, 1.5, 3.0):
        assert np.isclose(a.fanno_flow(M, gamma).friction_parameter, 4 * quad(rate, M, 1)[0])


def test_duct_flows_conserve_mass_momentum_and_energy():
    gamma = 1.4
    M = np.linspace(0.2, 3.0, 15)
    fanno, rayleigh = a.fanno_flow(M, gamma), a.rayleigh_flow(M, gamma)
    velocity = M * np.sqrt(fanno.temperature_ratio)
    assert np.allclose(fanno.density_ratio * velocity, 1)
    assert np.allclose(fanno.temperature_ratio * (1 + (gamma - 1) / 2 * M**2), (gamma + 1) / 2)
    assert np.allclose(rayleigh.density_ratio * M * np.sqrt(rayleigh.temperature_ratio), 1)
    assert np.allclose(rayleigh.pressure_ratio * (1 + gamma * M**2), 1 + gamma)
    assert np.all(rayleigh.stagnation_temperature_ratio <= 1 + 1e-12)


def test_panel_lift_agrees_with_kutta_joukowski_and_pressure_integration():
    for alpha in np.radians([-4.0, 2.0, 8.0]):
        result = a.panel_method(a.naca4_airfoil("4412", 121), alpha)
        assert np.isclose(result.cl, result.cl_pressure, atol=0.01)


def test_neutral_point_is_where_the_built_up_moment_slope_vanishes():
    h_ac, a_wb, a_t, V_H, eps, eta = 0.25, 4.9, 4.2, 0.55, 0.35, 0.9
    h_n = a.neutral_point(h_ac, a_wb, a_t, V_H, eps, eta)

    def moment(alpha, h):
        wing = a_wb * alpha * (h - h_ac)
        tail = -eta * V_H * a_t * (alpha * (1 - eps) - 0.02)
        return 0.05 * -1 + wing + tail

    for h in (0.2, h_n, 0.6):
        slope = (moment(0.01, h) - moment(-0.01, h)) / 0.02
        assert np.isclose(slope, a.pitching_moment_slope(h, h_ac, a_wb, a_t, V_H, eps, eta), atol=1e-12)
    assert np.isclose((moment(0.01, h_n) - moment(-0.01, h_n)) / 0.02, 0, atol=1e-12)


def test_takeoff_roll_integrates_newtons_law_to_lift_off():
    W, rho, S, Clmax, T, R, g = 10.7e3, 1.225, 16.2, 1.6, 2400.0, 400.0, 9.80665
    v_lo = 1.2 * a.stall_speed(W, rho, S, Clmax)
    acceleration = (T - R) * g / W
    assert np.isclose(a.takeoff_ground_roll(W, rho, S, Clmax, T, R), v_lo**2 / (2 * acceleration))


def test_propeller_breguet_range_matches_time_stepped_cruise():
    g, eta, psfc, LD, W0, Wf, V = 9.80665, 0.8, 7.6e-8, 13.0, 13.1e3, 11.7e3, 60.0
    W, x, dt = W0, 0.0, 1.0
    while W > Wf:
        W -= g * psfc * (W / LD) * V / eta * dt
        x += V * dt
    assert np.isclose(x, a.breguet_range_propeller(eta, psfc, LD, W0, Wf), rtol=1e-3)
