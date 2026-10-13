"""Reusable customizable plots and animations, independent of numerical solvers."""

import numpy as np


def plot_history(t, history, labels=None, ax=None):
    """Plot scalar/vector histories; return (figure, axes). Time is in seconds.

    State labels should include units, for example ``x [m]``.
    """
    import matplotlib.pyplot as plt

    t = np.asarray(t, dtype=float)
    values = np.asarray(history)
    if values.ndim == 1:
        values = values[:, None]
    if t.ndim != 1 or values.ndim != 2 or len(t) != len(values) or not len(t):
        raise ValueError("history must have shape (len(t), states)")
    if labels is not None and len(labels) != values.shape[1]:
        raise ValueError("one label is required per state")
    if ax is None:
        _, ax = plt.subplots()
    for i in range(values.shape[1]):
        ax.plot(t, values[:, i], label=None if labels is None else labels[i])
    ax.set(xlabel="time [s]", ylabel="state")
    if labels is not None:
        ax.legend()
    return ax.figure, ax


def animate_displacement(t, displacement, *, interval=30):
    """Animate 1D displacement [m]; return (figure, FuncAnimation).

    Input samples play at interval milliseconds/frame, independently of physical
    time. Retain the animation until displaying or saving it.
    """
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    t = np.asarray(t, dtype=float)
    x = np.asarray(displacement, dtype=float)
    if t.ndim != 1 or x.shape != t.shape or not len(t):
        raise ValueError("t and displacement must be matching nonempty vectors")
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(x)) or interval <= 0:
        raise ValueError("data must be finite and interval positive")
    fig, ax = plt.subplots()
    margin = max(float(np.ptp(x)) * 0.1, 0.01)
    ax.set(xlim=(x.min() - margin, x.max() + margin), ylim=(-1, 1), xlabel="displacement [m]")
    (marker,) = ax.plot([], [], "o", markersize=16)
    label = ax.text(0.02, 0.9, "", transform=ax.transAxes)

    def update(i):
        marker.set_data([x[i]], [0])
        label.set_text(f"t = {t[i]:.3f} s")
        return marker, label

    animation = FuncAnimation(fig, update, frames=len(t), interval=interval, blit=True)
    return fig, animation


__all__ = ["plot_history", "animate_displacement"]
