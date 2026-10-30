"""Plotting helpers for engineeringkit.reliability: survival and hazard curves and Weibull probability plots."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = ["plot_survival", "plot_hazard", "plot_weibull_probability"]


def _curve(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 1 or not x.size:
        raise ValueError("x and y must be matching nonempty vectors")
    return x, y


def plot_survival(times, survival, ax=None, label=None, step=True):
    """Plot a survival (reliability) curve against time.

    Parameters
    ----------
    times : array_like
        Times [s].
    survival : array_like
        Matching survival probabilities [-].
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.
    step : bool, optional
        Draw a right-continuous step function (for empirical estimates)
        instead of a line.

    Returns
    -------
    matplotlib.axes.Axes
    """
    x, y = _curve(times, survival)
    if ax is None:
        _, ax = plt.subplots()
    if step:
        ax.step(x, y, where="post", label=label)
    else:
        ax.plot(x, y, label=label)
    ax.set(xlabel="time", ylabel="survival probability", ylim=(0, 1.05))
    if label:
        ax.legend()
    return ax


def plot_hazard(times, hazard, ax=None, label=None, log=False):
    """Plot a hazard rate against time.

    Parameters
    ----------
    times : array_like
        Times [s].
    hazard : array_like
        Matching hazard rates [1/s].
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label.
    log : bool, optional
        Use a logarithmic hazard axis.

    Returns
    -------
    matplotlib.axes.Axes
    """
    x, y = _curve(times, hazard)
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(x, y, label=label)
    if log:
        ax.set_yscale("log")
    ax.set(xlabel="time", ylabel="hazard rate")
    if label:
        ax.legend()
    return ax


def plot_weibull_probability(times, failure_probability, eta=None, beta=None, ax=None, label=None):
    r"""Plot life data on Weibull probability axes, with an optional fitted line.

    The axes are :math:`\ln t` and :math:`\ln(-\ln(1 - F))`, the scales of
    Weibull probability paper, on which a Weibull distribution is a straight
    line of slope :math:`\beta` crossing zero at :math:`t = \eta`. The
    right axis is labelled in failure probability.

    Parameters
    ----------
    times : array_like
        Failure times [s], positive.
    failure_probability : array_like
        Matching plotting positions (median ranks) [-], in (0, 1).
    eta, beta : float, optional
        Scale [s] and shape [-] of a fitted line to draw; both or neither.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    label : str, optional
        Legend label of the fitted line.

    Returns
    -------
    matplotlib.axes.Axes
    """
    t, f = _curve(times, failure_probability)
    if np.any(t <= 0) or np.any((f <= 0) | (f >= 1)):
        raise ValueError("times must be positive and failure probabilities in (0, 1)")
    if (eta is None) != (beta is None):
        raise ValueError("give both eta and beta, or neither")
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(np.log(t), np.log(-np.log1p(-f)), "o", label="data")
    if eta is not None:
        x = np.linspace(np.log(t.min()) - 0.5, np.log(t.max()) + 0.5, 2)
        ax.plot(x, beta * (x - np.log(eta)), "-", label=label or "fit")
        ax.legend()
    ax.set(xlabel="ln(time)", ylabel="ln(-ln(1 - F))", title="Weibull probability plot")
    ticks = np.array([0.01, 0.05, 0.1, 0.3, 0.5, 0.632, 0.9, 0.99])
    right = ax.secondary_yaxis("right", functions=(lambda y: y, lambda y: y))
    right.set_yticks(np.log(-np.log1p(-ticks)), [f"{q:g}" for q in ticks])
    right.set_ylabel("failure probability F")
    return ax
