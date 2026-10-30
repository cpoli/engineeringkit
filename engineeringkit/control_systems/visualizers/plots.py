"""Plotting helpers for engineeringkit.control_systems: time responses, Bode, Nyquist, Nichols, root-locus and phase-plane plots."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.control_systems.systems.frequency import bode, nichols, nichols_closed_loop, nyquist
from engineeringkit.control_systems.systems.root_locus import root_locus
from engineeringkit.control_systems.utils.polynomials import as_polynomial

__all__ = ["plot_response", "plot_bode", "plot_nyquist", "plot_nichols", "plot_root_locus", "plot_phase_portrait"]


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


NICHOLS_CONTOURS_DB = (-20.0, -12.0, -6.0, -3.0, -1.0, 0.0, 0.5, 1.0, 3.0, 6.0)


def plot_nichols(sys, w=None, ax=None, label=None, contours=NICHOLS_CONTOURS_DB):
    """Nichols chart: open-loop gain against phase over constant closed-loop gain (M) contours.

    Parameters
    ----------
    sys : scipy.signal.lti
        SISO open-loop model.
    w : array_like, optional
        Positive angular frequencies [rad/s], as in :func:`nichols`.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.
    contours : sequence of float or None, optional
        Closed-loop gains [dB] to draw as contours; ``None`` draws none.

    Returns
    -------
    matplotlib.axes.Axes
    """
    n = nichols(sys, w)
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(n["phase_deg"], n["gain_db"], label=label)
    if contours is not None:
        low = min(-360.0, 360 * np.floor(np.min(n["phase_deg"]) / 360))
        high = max(0.0, 360 * np.ceil(np.max(n["phase_deg"]) / 360))
        phase, gain = np.meshgrid(np.linspace(low, high, 721), np.linspace(-40, 40, 400))
        closed, _ = nichols_closed_loop(gain, phase)
        cs = ax.contour(phase, gain, closed, levels=sorted(contours), colors="gray", linewidths=0.6)
        ax.clabel(cs, fmt="%g dB", fontsize=7)
    ax.plot([-180], [0], "r+", ms=12, mew=2)
    ax.set(xlabel="open-loop phase [deg]", ylabel="open-loop gain [dB]")
    ax.grid(True, alpha=0.3)
    if label:
        ax.legend()
    return ax


def plot_phase_portrait(x, v, acceleration=None, ax=None, label=None):
    """Phase-plane trajectory, optionally over the direction field of :math:`\\ddot x = f(x, \\dot x)`.

    Parameters
    ----------
    x, v : array_like
        Matching position and velocity samples.
    acceleration : callable, optional
        ``acceleration(x, v)``, evaluated on arrays, to draw normalised
        direction arrows over the plotted region.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.

    Returns
    -------
    matplotlib.axes.Axes
    """
    x, v = np.asarray(x, float), np.asarray(v, float)
    if x.shape != v.shape or x.ndim != 1 or not x.size:
        raise ValueError("x and v must be matching nonempty vectors")
    if ax is None:
        _, ax = plt.subplots()
    if acceleration is not None:
        pad_x, pad_v = 0.1 * np.ptp(x) + 1e-9, 0.1 * np.ptp(v) + 1e-9
        X, V = np.meshgrid(np.linspace(x.min() - pad_x, x.max() + pad_x, 21), np.linspace(v.min() - pad_v, v.max() + pad_v, 21))
        dx, dv = V / (np.ptp(X) or 1.0), np.broadcast_to(acceleration(X, V), X.shape) / (np.ptp(V) or 1.0)
        norm = np.hypot(dx, dv)
        norm[norm == 0] = 1.0
        ax.quiver(X, V, dx / norm, dv / norm, color="0.7", angles="xy", pivot="mid", width=0.002)
    ax.plot(x, v, label=label)
    ax.plot(x[:1], v[:1], "ko", ms=4)
    ax.set(xlabel="x", ylabel="dx/dt")
    if label:
        ax.legend()
    return ax
