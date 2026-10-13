import numpy as np
import pytest

from engineeringkit import aerodynamics as a


def test_loads_coefficients_and_speed_recover_weight():
    q = a.dynamic_pressure(1.2, 50)
    assert q == 1500
    assert a.lift(q, 2, 0.5) == 1500
    assert a.drag(q, 2, 0.02) == 60
    assert a.lift_coefficient(1500, 1.2, 50, 2) == 0.5
    assert a.drag_coefficient(60, 1.2, 50, 2) == 0.02
    assert a.aspect_ratio(10, 20) == 5
    speed = a.stall_speed(12000, 1.2, 20, 1.5)
    assert np.isclose(a.lift(a.dynamic_pressure(1.2, speed), 20, 1.5), 12000)
    assert a.mach_number(170, 340) == 0.5
    assert np.isclose(a.finite_wing_cl(0.1, 8), a.finite_wing_lift_slope(8) * 0.1)


def test_compressible_state_and_flat_plate_models():
    assert np.isclose(a.stagnation_temperature(300, 2), 540)
    assert np.isclose(a.stagnation_pressure(100, 2) / 100, 1.8**3.5)
    assert np.isclose(a.prandtl_glauert_cp(-1, 0.6), -1.25)
    assert np.isclose(a.reynolds_number(1, 10, 1, 1e-5), 1e6)
    assert np.isclose(a.skin_friction_laminar(1e6), 0.001328)
    assert np.isclose(a.skin_friction_turbulent(1e5), 0.0074)
    assert np.isclose(a.boundary_layer_laminar(1, 1e6), 0.005)
    assert 0 < a.oswald_efficiency_estimate(8) < 1


@pytest.mark.parametrize(
    "fn,args",
    [
        (a.best_glide_cl, (0.02, 8, 0)),
        (a.normal_shock_mach2, (2, 1)),
        (a.normal_shock_pressure_ratio, (2, 1)),
        (a.normal_shock_pressure_ratio, (0.5, 1.4)),
        (a.area_mach_ratio, (2, 1)),
        (a.stagnation_pressure, (1, 2, 1)),
    ],
)
def test_invalid_compressible_and_wing_parameters(fn, args):
    with pytest.raises(ValueError):
        fn(*args)
