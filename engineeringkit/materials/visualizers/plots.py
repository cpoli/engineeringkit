"""Plotting helpers for engineeringkit.materials: stress-strain curves, yield loci, S-N curves and crack growth."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit._validation import positive as _positive

__all__ = ["plot_stress_strain", "plot_yield_surfaces", "plot_sn_curve", "plot_crack_growth"]


def _matching(x, y, names):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 1 or not x.size:
        raise ValueError(f"{names} must be matching nonempty vectors")
    return x, y


def plot_stress_strain(strain, stress, ax=None, label=None):
    """Plot stress in MPa against strain.

    Parameters
    ----------
    strain : array_like
        Strains [-].
    stress : array_like
        Matching stresses [Pa].
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.

    Returns
    -------
    matplotlib.axes.Axes
    """
    x, y = _matching(strain, stress, "strain and stress")
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(x, y / 1e6, label=label)
    ax.set(xlabel="strain", ylabel="stress [MPa]")
    ax.grid(True, alpha=0.3)
    if label:
        ax.legend()
    return ax


def plot_yield_surfaces(yield_strength, ax=None):
    """Draw the Tresca hexagon and the von Mises ellipse in the plane-stress principal plane.

    Parameters
    ----------
    yield_strength : float
        Uniaxial yield strength [Pa], positive; the axes are scaled by it.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    _positive(yield_strength=yield_strength)
    if ax is None:
        _, ax = plt.subplots()
    hexagon = np.array([[1, 0], [1, 1], [0, 1], [-1, 0], [-1, -1], [0, -1], [1, 0]])
    t = np.linspace(0, 2 * np.pi, 361)
    # sigma1^2 - sigma1 sigma2 + sigma2^2 = 1 in rotated coordinates.
    u, v = np.sqrt(2) * np.cos(t), np.sqrt(2 / 3) * np.sin(t)
    ax.plot(hexagon[:, 0], hexagon[:, 1], label="Tresca")
    ax.plot((u - v) / np.sqrt(2), (u + v) / np.sqrt(2), label="von Mises")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.axvline(0, color="0.6", lw=0.8)
    ax.set(xlabel="sigma_1 / sigma_y", ylabel="sigma_2 / sigma_y", aspect="equal")
    ax.legend()
    return ax


def plot_sn_curve(cycles, stress_amplitude, ax=None, label=None, style="-"):
    """Log-log plot of stress amplitude in MPa against cycles to failure.

    Parameters
    ----------
    cycles : array_like
        Cycles to failure, positive.
    stress_amplitude : array_like
        Matching stress amplitudes [Pa], positive.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.
    style : str
        Matplotlib format string, e.g. ``"o"`` for test points.

    Returns
    -------
    matplotlib.axes.Axes
    """
    x, y = _matching(cycles, stress_amplitude, "cycles and stress_amplitude")
    if np.any(x <= 0) or np.any(y <= 0):
        raise ValueError("cycles and stress_amplitude must be positive")
    if ax is None:
        _, ax = plt.subplots()
    ax.loglog(x, y / 1e6, style, label=label)
    ax.set(xlabel="cycles to failure N", ylabel="stress amplitude [MPa]")
    ax.grid(True, which="both", alpha=0.3)
    if label:
        ax.legend()
    return ax


def plot_crack_growth(cycles, crack_length, ax=None, label=None):
    """Plot crack length in mm against cycles.

    Parameters
    ----------
    cycles : array_like
        Cycles.
    crack_length : array_like
        Matching crack lengths [m].
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.

    Returns
    -------
    matplotlib.axes.Axes
    """
    x, y = _matching(cycles, crack_length, "cycles and crack_length")
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(x, 1e3 * y, label=label)
    ax.set(xlabel="cycles", ylabel="crack length [mm]")
    ax.grid(True, alpha=0.3)
    if label:
        ax.legend()
    return ax
