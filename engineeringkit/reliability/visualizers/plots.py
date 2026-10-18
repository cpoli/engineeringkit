"""Plotting helpers for engineeringkit.reliability: survival and hazard curves."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = ["plot_survival", "plot_hazard"]


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
