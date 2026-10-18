import numpy as np

from engineeringkit import mechanical_engineering as m


def test_gear_power_with_efficiency():
    omega, torque = 100.0, 10.0
    out_omega = m.gear_output_speed(omega, 20, 60)
    out_torque = m.gear_output_torque(torque, 20, 60, 0.9)
    assert np.isclose(m.power_torque(out_torque, out_omega), 0.9 * m.power_torque(torque, omega))


def test_parallel_axis_and_slider_crank_strokes():
    assert m.parallel_axis(m.disk_inertia(2, 1), 2, 1) == 3
    positions = m.slider_crank_displacement(1, 3, np.array([0, np.pi]))
    np.testing.assert_allclose(positions, [4, 2])
    assert np.isclose(m.damped_frequency(4, 1, 4), 0)


def test_forced_static_displacement_and_dissipative_decay():
    r = m.simulate_sdof(1, 4, 4, t=np.linspace(0, 10, 1001), force=lambda t: 8.0)
    assert np.isclose(r.displacement[-1], 2, atol=1e-6)
    r = m.simulate_sdof(1, 1, 4, x0=1, t=np.linspace(0, 10, 1001))
    energy = 0.5 * r.velocity**2 + 2 * r.displacement**2
    assert np.all(np.diff(energy) <= 1e-8)
