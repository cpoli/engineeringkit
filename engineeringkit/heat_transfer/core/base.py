"""Result containers for engineeringkit.heat_transfer.

Most heat-transfer models are closed forms that return plain floats or
arrays. The steady composite wall and the finite-difference conduction
solver return more than one number, so they get frozen dataclasses with the
usual ``summary()`` and ``plot()``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["WallResult", "ConductionResult"]


@dataclass(frozen=True)
class WallResult:
    """Output of :func:`engineeringkit.heat_transfer.composite_wall`."""

    heat_rate: float
    """float: Heat rate through the wall [W], positive from hot to cold."""

    interface_temperatures: np.ndarray
    """ndarray: Temperature at each node of the resistance chain, hot side first [K]."""

    def summary(self) -> dict:
        """Heat rate [W] and interface temperatures [K]."""
        return {"heat_rate_W": self.heat_rate, "interfaces_K": self.interface_temperatures.tolist()}

    def plot(self, ax=None):
        """Draw the temperature at each interface; return ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        ax.plot(np.arange(len(self.interface_temperatures)), self.interface_temperatures, marker="o")
        ax.set(xlabel="interface", ylabel="temperature [K]")
        return ax.figure, ax


@dataclass(frozen=True)
class ConductionResult:
    """Output of :func:`engineeringkit.heat_transfer.conduction_1d`."""

    x: np.ndarray
    """ndarray: Node positions [m], shape ``(n,)``."""

    t: np.ndarray
    """ndarray: Times [s], shape ``(steps + 1,)``."""

    T: np.ndarray
    """ndarray: Temperatures [K], shape ``(steps + 1, n)``; row ``i`` is the profile at ``t[i]``."""

    method: str
    """str: ``"explicit"``, ``"implicit"`` or ``"crank_nicolson"``."""

    mesh_fourier: float
    r"""float: Mesh Fourier number :math:`r = \alpha\,\Delta t/\Delta x^2` [-]; the explicit scheme needs :math:`r \le 1/2`."""

    def summary(self) -> dict:
        """Grid size, scheme, mesh Fourier number and final temperature range."""
        return {
            "nodes": int(self.x.size),
            "steps": int(self.t.size - 1),
            "method": self.method,
            "mesh_fourier": self.mesh_fourier,
            "final_min_K": float(self.T[-1].min()),
            "final_max_K": float(self.T[-1].max()),
        }

    def plot(self, ax=None, profiles=6):
        """Draw evenly spaced temperature profiles; return ``(fig, ax)``.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure is created if omitted.
        profiles : int, optional
            Number of time instants to draw, first and last included.

        Returns
        -------
        tuple
            ``(fig, ax)``.
        """
        from engineeringkit.heat_transfer.visualizers.plots import plot_temperature_profiles

        rows = np.unique(np.linspace(0, self.t.size - 1, profiles).round().astype(int))
        ax = plot_temperature_profiles(self.x, self.T[rows], self.t[rows], ax=ax)
        return ax.figure, ax
