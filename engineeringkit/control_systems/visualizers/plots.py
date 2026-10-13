"""Plotting helpers for engineeringkit.control_systems: time responses, Bode, Nyquist and root-locus plots."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.control_systems.systems.frequency import bode, nyquist
from engineeringkit.control_systems.systems.root_locus import root_locus
from engineeringkit.control_systems.utils.polynomials import as_polynomial

__all__ = ["plot_response", "plot_bode", "plot_nyquist", "plot_root_locus"]


def plot_response(t, y, ax=None, label=None):
    """Plot a time response.

    Parameters
    ----------
    t, y : array_like
        Matching time samples [s] and output values.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.

    Returns
    -------
    matplotlib.axes.Axes
    """
    t, y = np.asarray(t, float), np.asarray(y, float)
    if t.shape != y.shape or t.ndim != 1 or not t.size:
        raise ValueError("t and y must be matching nonempty vectors")
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(t, y, label=label)
    ax.set(xlabel="time [s]", ylabel="response")
    if label:
        ax.legend()
    return ax


def plot_bode(sys, w=None, axes=None, label=None):
    """Bode magnitude and phase plots on a shared logarithmic frequency axis.

    Parameters
    ----------
    sys : scipy.signal.lti
        Continuous-time model.
    w : array_like, optional
        Angular frequencies [rad/s].
    axes : sequence of two matplotlib.axes.Axes, optional
        Magnitude and phase axes; a new figure is created if omitted.
    label : str, optional
        Legend label.

    Returns
    -------
    tuple of matplotlib.axes.Axes
        ``(magnitude_axes, phase_axes)``.
    """
    w, mag, phase = bode(sys, w)
    if axes is None:
        _, axes = plt.subplots(2, 1, sharex=True)
    ax_mag, ax_phase = axes
    ax_mag.semilogx(w, mag, label=label)
    ax_phase.semilogx(w, phase, label=label)
    ax_mag.set(ylabel="magnitude [dB]")
    ax_phase.set(xlabel="angular frequency [rad/s]", ylabel="phase [deg]")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.3)
    if label:
        ax_mag.legend()
    return ax_mag, ax_phase


def plot_nyquist(sys, w=None, ax=None, label=None):
    """Nyquist plot of a loop for positive and negative frequencies, marking -1.

    Parameters
    ----------
    sys : scipy.signal.lti
        SISO open-loop model.
    w : array_like, optional
        Positive angular frequencies [rad/s]; defaults to a log grid from
        :math:`10^{-3}` to :math:`10^3`.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.

    Returns
    -------
    matplotlib.axes.Axes
    """
    w = np.logspace(-3, 3, 2000) if w is None else np.asarray(w, float)
    h = nyquist(sys, w)
    if ax is None:
        _, ax = plt.subplots()
    (line,) = ax.plot(h.real, h.imag, label=label)
    ax.plot(h.real, -h.imag, "--", color=line.get_color(), lw=0.8)
    ax.plot([-1], [0], "r+", ms=12, mew=2)
    ax.set(xlabel="real", ylabel="imaginary")
    ax.grid(True, alpha=0.3)
    if label:
        ax.legend()
    return ax


def plot_root_locus(num, den, gains, ax=None):
    """Plot closed-loop root loci with open-loop poles (x) and zeros (o).

    Parameters
    ----------
    num, den : array_like
        Strictly proper open-loop numerator and denominator.
    gains : array_like
        Loop gains to evaluate.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
    """
    roots = root_locus(num, den, gains)
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(roots.real, roots.imag, ".", ms=2, color="C0")
    poles, zeros = np.roots(as_polynomial(den, "den")), np.roots(as_polynomial(num, "num"))
    ax.plot(poles.real, poles.imag, "kx", ms=9, label="open-loop poles")
    if zeros.size:
        ax.plot(zeros.real, zeros.imag, "ko", mfc="none", ms=9, label="open-loop zeros")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.set(xlabel="real [1/s]", ylabel="imaginary [1/s]")
    ax.legend()
    return ax
