import numpy as np
import pytest

from engineeringkit import mechanical_engineering as m


def test_rigid_body_and_bearing_reference_values():
    assert m.kinetic_energy(2, 3) == 9
    assert m.rotational_energy(2, 3) == 9
    assert np.isclose(m.work_constant_force(10, 2, np.pi / 3), 10)
    assert m.angular_momentum(2, 3) == 6
    assert m.centripetal_force(2, 3, 1) == 18
    assert np.isclose(m.belt_speed(1, 60), np.pi)
    assert m.bearing_l10_life(2, 1) == 8e6
    assert m.flywheel_energy(2, 3, 1) == 8
    assert np.isclose(m.natural_frequency_hz(4, 1), 1 / np.pi)
    assert m.transmissibility(0, 0.1) == 1
    assert np.isinf(m.transmissibility(1, 0))
    assert m.rotating_unbalance_force(2, 0.1, 10) == 20
    assert m.fourbar_grashof([1, 2, 2, 3])
    assert not m.fourbar_grashof([1, 2, 2, 4])
    assert m.solid_sphere_inertia(5, 1) == 2


def test_square_thread_screw_and_geometry_errors():
    # No friction: work per revolution equals load * lead.
    torque = m.screw_torque_raise(100, 0.02, 0.01, 0)
    assert np.isclose(2 * np.pi * torque, 100 * 0.01)
    with pytest.raises(ValueError):
        m.screw_torque_raise(100, 0.01, 1, 100)
    with pytest.raises(ValueError):
        m.fourbar_grashof([1, 2, 3])


def test_sdof_summary_plot_and_animation(tmp_path):
    import matplotlib.pyplot as plt

    r = m.simulate_sdof(1, 0, 4, x0=0.1, t=np.linspace(0, 1, 5))
    assert r.summary()["peak_displacement_m"] == 0.1
    fig, ax = r.plot()
    assert len(ax.lines) == 1
    plt.close(fig)
    fig, animation = r.animate(interval=100)
    animation.save(tmp_path / "oscillator.gif", writer="pillow", fps=5)
    assert (tmp_path / "oscillator.gif").stat().st_size > 0
    plt.close(fig)
