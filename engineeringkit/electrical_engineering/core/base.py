"""Result containers for engineeringkit.electrical_engineering.

Most circuit models are closed-form relations that return plain floats,
complex numbers or arrays. The three solvers that produce structured output
(a transient response, a nodal analysis and a load flow) return the frozen
dataclasses below, each with ``summary()`` and ``plot()``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["CircuitResult", "NodalResult", "LoadFlowResult"]


@dataclass(frozen=True)
class CircuitResult:
    """Time response of a circuit, from :meth:`~engineeringkit.electrical_engineering.SeriesRLC.step_response`."""

    t: np.ndarray
    """ndarray: Time samples [s]."""

    response: np.ndarray
    """ndarray: Loop current at each sample [A]."""

    def summary(self) -> dict:
        """Peak and final current [A]."""
        return {"peak": float(np.max(self.response)), "final": float(self.response[-1])}

    def plot(self, ax=None):
        """Plot current against time; return ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        ax.plot(self.t, self.response)
        ax.set(xlabel="time [s]", ylabel="current [A]")
        return ax.figure, ax


@dataclass(frozen=True)
class NodalResult:
    """Output of :func:`engineeringkit.electrical_engineering.nodal_analysis`."""

    voltages: np.ndarray
    """ndarray: Node voltages [V], index 0 is the ground node (always 0)."""

    resistor_currents: np.ndarray
    """ndarray: Current in each resistor from its first node to its second [A]."""

    source_currents: np.ndarray
    """ndarray: Current delivered by each voltage source out of its + terminal [A]."""

    kcl_residual: float
    """float: Largest net current leaving any non-ground node [A]; zero up to rounding."""

    def summary(self) -> dict:
        """Number of nodes, largest node voltage magnitude [V] and the KCL residual [A]."""
        return {
            "nodes": int(self.voltages.size),
            "max_voltage": float(np.max(np.abs(self.voltages))),
            "kcl_residual": float(self.kcl_residual),
        }

    def plot(self, ax=None):
        """Bar chart of node voltages; return ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        ax.bar(np.arange(self.voltages.size), self.voltages)
        ax.set(xlabel="node", ylabel="voltage [V]")
        return ax.figure, ax


@dataclass(frozen=True)
class LoadFlowResult:
    """Output of :func:`engineeringkit.electrical_engineering.gauss_seidel_load_flow`.

    All quantities are per unit on the system base.
    """

    voltages: np.ndarray
    """ndarray of complex: Bus voltages, bus 0 is the slack bus."""

    power: np.ndarray
    """ndarray of complex: Net complex power injected at each bus, :math:`S_i = V_i I_i^*`."""

    iterations: int
    """int: Gauss-Seidel sweeps performed."""

    converged: bool
    """bool: Whether the largest voltage update fell below the tolerance."""

    history: np.ndarray
    """ndarray: Largest voltage update in each sweep."""

    @property
    def losses(self) -> complex:
        """Total network losses, the sum of all bus injections (per unit)."""
        return complex(np.sum(self.power))

    def summary(self) -> dict:
        """Convergence, smallest voltage magnitude, slack power and losses (per unit)."""
        return {
            "converged": bool(self.converged),
            "iterations": int(self.iterations),
            "min_voltage": float(np.min(np.abs(self.voltages))),
            "slack_power": complex(self.power[0]),
            "losses": self.losses,
        }

    def plot(self, ax=None):
        """Bar chart of bus voltage magnitudes; return ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        ax.bar(np.arange(self.voltages.size), np.abs(self.voltages))
        ax.set(xlabel="bus", ylabel="voltage magnitude [pu]")
        return ax.figure, ax
