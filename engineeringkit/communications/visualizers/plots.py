"""Plotting helpers for engineeringkit.communications: constellations and error-rate curves."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = ["plot_constellation", "plot_ber_curve"]


def plot_constellation(symbols, ax=None, label=None):
    """Scatter complex symbols in the in-phase/quadrature plane.

    Parameters
    ----------
    symbols : array_like
        Complex (or real) symbols.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.

    Returns
    -------
    matplotlib.axes.Axes
    """
    s = np.asarray(symbols, dtype=complex).ravel()
    if not s.size:
        raise ValueError("symbols must be nonempty")
    if ax is None:
        _, ax = plt.subplots(figsize=(4, 4))
    ax.plot(s.real, s.imag, ".", label=label)
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.set_aspect("equal")
    ax.set(xlabel="in-phase", ylabel="quadrature")
    if label:
        ax.legend()
    return ax


def plot_ber_curve(ebn0_db, ber, ax=None, label=None):
    """Semilog plot of bit-error rate against :math:`E_b/N_0` in decibels.

    Parameters
    ----------
    ebn0_db : array_like
        :math:`E_b/N_0` values [dB].
    ber : array_like
        Matching bit-error rates [-]; zeros are omitted from the log axis.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.

    Returns
    -------
    matplotlib.axes.Axes
    """
    x, y = np.asarray(ebn0_db, dtype=float), np.asarray(ber, dtype=float)
    if x.shape != y.shape or x.ndim != 1 or not x.size:
        raise ValueError("ebn0_db and ber must be matching nonempty vectors")
    if ax is None:
        _, ax = plt.subplots()
    ax.semilogy(x, np.where(y > 0, y, np.nan), label=label)
    ax.set(xlabel="$E_b/N_0$ [dB]", ylabel="bit-error rate")
    ax.grid(True, which="both", alpha=0.3)
    if label:
        ax.legend()
    return ax
