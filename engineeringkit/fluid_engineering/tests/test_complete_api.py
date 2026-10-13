import numpy as np
import pytest

from engineeringkit import fluid_engineering as f


def test_velocity_pressure_and_loss_conversions():
    assert f.colebrook_friction_factor(1000, 0) == 0.064
    assert np.isclose(f.swamee_jain_friction_factor(1e5, 0.0001), f.colebrook_friction_factor(1e5, 0.0001), rtol=0.01)
    assert f.dynamic_pressure(1000, 2) == 2000
    assert f.pitot_velocity(2000, 1000) == 2
    assert f.cavitation_number(3000, 1000, 1000, 2) == 1
    assert f.mach_number(170, 340) == 0.5
    assert f.hydraulic_diameter(np.pi, 0.5 * 4 * np.pi) == 2
    assert np.isclose(f.pressure_head(1000 * f.G, 1000), 1)
    assert f.euler_number(4000, 1000, 2) == 1
    assert np.isclose(f.pump_specific_speed(10, 4, 1), 20 / f.G**0.75)
    assert np.isclose(f.minor_head_loss(2, 3), 9 / f.G)
    with pytest.raises(ValueError):
        f.pump_input_power(1000, 0.1, 10, 0)
