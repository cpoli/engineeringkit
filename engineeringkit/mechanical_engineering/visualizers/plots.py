"""Plotting helpers for engineeringkit.mechanical_engineering: frequency responses, mode shapes and linkages."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = ["plot_frequency_response", "plot_mode_shapes", "plot_linkage"]


def plot_frequency_response(r, amplitude, ax=None, label=None):
    """Semilog plot of an amplitude ratio against a frequency ratio.

    Parameters
    ----------
    r : array_like
        Frequency ratios [-].
    amplitude : array_like
        Matching amplitude ratios [-], positive; infinite or nonpositive
        values are left out of the log axis.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.

    Returns
    -------
    matplotlib.axes.Axes
    """
    x, y = np.asarray(r, dtype=float), np.asarray(amplitude, dtype=float)
    if x.shape != y.shape or x.ndim != 1 or not x.size:
        raise ValueError("r and amplitude must be matching nonempty vectors")
    if ax is None:
        _, ax = plt.subplots()
    ax.semilogy(x, np.where(np.isfinite(y) & (y > 0), y, np.nan), label=label)
    ax.set(xlabel="frequency ratio", ylabel="amplitude ratio")
    ax.grid(True, which="both", alpha=0.3)
    if label:
        ax.legend()
    return ax


def plot_mode_shapes(mode_shapes, frequencies, ax=None):
    """Plot each mode's disk amplitudes against disk number.

    Parameters
    ----------
    mode_shapes : array_like
        Shape ``(modes, disks)``.
    frequencies : array_like
        One natural frequency per mode [rad/s], used in the legend.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    shapes, w = np.atleast_2d(np.asarray(mode_shapes, dtype=float)), np.asarray(frequencies, dtype=float).ravel()
    if shapes.shape[0] != w.size or not w.size:
        raise ValueError("one frequency is required per mode shape")
    if ax is None:
        _, ax = plt.subplots()
    disks = np.arange(1, shapes.shape[1] + 1)
    for shape, frequency in zip(shapes, w):
        ax.plot(disks, shape, "o-", label=f"{frequency / (2 * np.pi):.3g} Hz")
    ax.axhline(0, color="gray", lw=0.5)
    ax.set(xlabel="disk", ylabel="relative twist", xticks=disks)
    ax.legend()
    return ax


def plot_linkage(points, ax=None, **style):
    """Draw a chain of rigid links through successive joint positions.

    Parameters
    ----------
    points : array_like
        Joint coordinates [m], shape ``(n, 2)`` with ``n >= 2``.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    **style
        Passed to :meth:`matplotlib.axes.Axes.plot`.

    Returns
    -------
    matplotlib.axes.Axes
    """
    p = np.asarray(points, dtype=float)
    if p.ndim != 2 or p.shape[1] != 2 or len(p) < 2:
        raise ValueError("points must have shape (n, 2) with n >= 2")
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(p[:, 0], p[:, 1], **{"marker": "o", "lw": 2, **style})
    ax.set_aspect("equal")
    ax.set(xlabel="x [m]", ylabel="y [m]")
    return ax
