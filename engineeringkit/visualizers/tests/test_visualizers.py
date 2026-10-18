import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.integrators import integrate_ode
from engineeringkit.visualizers import animate_displacement, plot_history


def test_plot_uses_supplied_axes(tmp_path):
    fig, ax = plt.subplots()
    out, out_ax = plot_history([0, 1], [[0, 1], [1, 0]], ["x [m]", "v [m/s]"], ax)
    assert out is fig and out_ax is ax and len(ax.lines) == 2
    out.savefig(tmp_path / "history.png")
    plt.close(fig)


def test_animation_and_result_plot(tmp_path):
    r = integrate_ode(lambda t, y: -y, [1], np.linspace(0, 1, 5))
    fig, ax = r.plot(["x [m]"])
    assert len(ax.lines) == 1
    plt.close(fig)
    fig, animation = animate_displacement(r.t, r.history[:, 0])
    animation.save(tmp_path / "motion.gif", writer="pillow", fps=5)
    assert (tmp_path / "motion.gif").stat().st_size > 0
    plt.close(fig)
