import numpy as np
import pytest

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
