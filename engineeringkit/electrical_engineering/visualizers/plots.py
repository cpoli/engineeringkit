"""Plotting helpers for engineeringkit.electrical_engineering: phasor diagrams and the Smith chart."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = ["plot_phasors", "plot_smith_chart"]


def plot_phasors(phasors, labels=None, ax=None):
    """Draw phasors as arrows from the origin of the complex plane.

    Parameters
    ----------
    phasors : array_like of complex
        Phasors to draw [V, A or other], at least one.
    labels : sequence of str, optional
        One legend label per phasor.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    z = np.asarray(phasors, dtype=complex).ravel()
    if not z.size or not np.all(np.isfinite(z)):
        raise ValueError("phasors must be a nonempty array of finite values")
    if labels is not None and len(labels) != z.size:
        raise ValueError("labels must have one entry per phasor")
    if ax is None:
        _, ax = plt.subplots(figsize=(4, 4))
    for k, p in enumerate(z):
        ax.annotate("", xy=(p.real, p.imag), xytext=(0, 0), arrowprops={"arrowstyle": "->", "color": f"C{k % 10}", "lw": 1.5})
        ax.plot([], [], color=f"C{k % 10}", label=None if labels is None else labels[k])
    r = 1.15 * float(np.max(np.abs(z))) or 1.0
    ax.set(xlim=(-r, r), ylim=(-r, r), xlabel="real", ylabel="imaginary")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.set_aspect("equal")
    if labels is not None:
        ax.legend(loc="upper left", fontsize="small")
    return ax


def plot_smith_chart(gamma=None, ax=None, label=None):
    """Draw a Smith chart grid and, optionally, reflection coefficients on it.

    Circles of constant normalized resistance :math:`r` and arcs of constant
    normalized reactance :math:`x` are drawn in the :math:`\\Gamma` plane.

    Parameters
    ----------
    gamma : array_like of complex, optional
        Reflection coefficients [-] to plot as a line (a single value as a
        point), all with :math:`|\\Gamma| \\le 1`.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label for ``gamma``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))
    t = np.linspace(0, 2 * np.pi, 361)
    grid = {"color": "0.75", "lw": 0.6}
    ax.plot(np.cos(t), np.sin(t), color="0.3", lw=1)
    for r in (0.2, 0.5, 1, 2, 5):
        ax.plot(r / (1 + r) + np.cos(t) / (1 + r), np.sin(t) / (1 + r), **grid)
    for x in (0.2, 0.5, 1, 2, 5):
        zx = 1j * x + np.concatenate(([0], np.logspace(-3, 3, 400)))
        for sign in (1, -1):
            g = (zx.real + 1j * sign * zx.imag - 1) / (zx.real + 1j * sign * zx.imag + 1)
            ax.plot(g.real, g.imag, **grid)
    ax.plot([-1, 1], [0, 0], **grid)
    if gamma is not None:
        g = np.asarray(gamma, dtype=complex).ravel()
        if not g.size or not np.all(np.isfinite(g)) or np.any(np.abs(g) > 1 + 1e-12):
            raise ValueError("gamma must be finite reflection coefficients with |gamma| <= 1")
        ax.plot(g.real, g.imag, "o" if g.size == 1 else "-", label=label)
        if label:
            ax.legend(loc="upper left", fontsize="small")
    ax.set(xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))
    ax.set_aspect("equal")
    ax.axis("off")
    return ax
