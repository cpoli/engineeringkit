"""Result containers for engineeringkit.control_systems.

Most control models return SciPy LTI objects, arrays or dictionaries. A
simulated closed-loop response, a compensator design and a phase-plane
trajectory are returned as the frozen dataclasses below so they can be
summarised and plotted consistently.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np

__all__ = ["ResponseResult", "CompensatorResult", "PhasePortraitResult"]


@dataclass(frozen=True)
class ResponseResult:
    """Time response of a simulated control loop."""

    t: np.ndarray
    """ndarray: Time samples [s]."""

    y: np.ndarray
    """ndarray: Output at each sample, in output units per unit reference."""

    u: np.ndarray | None = None
    """ndarray or None: Actuator output at each sample, when the simulation records it."""

    def summary(self) -> dict:
        """Final value and peak of the response, and the input range when recorded."""
        out = {"final": float(self.y[-1]), "peak": float(np.max(self.y))}
        if self.u is not None:
            out.update(input_min=float(np.min(self.u)), input_max=float(np.max(self.u)))
        return out

    def plot(self, ax=None):
        """Plot the response against time; return ``(fig, ax)``."""
        from engineeringkit.control_systems.visualizers.plots import plot_response

        ax = plot_response(self.t, self.y, ax=ax)
        return ax.figure, ax


@dataclass(frozen=True)
class CompensatorResult:
    """A lead, lag or lead-lag compensator designed for a loop, with the margins it achieves."""

    kind: str
    """str: ``"lead"``, ``"lag"`` or ``"lead-lag"``."""

    compensator: Any
    """scipy.signal.TransferFunction: Compensator :math:`C(s)`, unit gain at :math:`s = 0`."""

    loop: Any
    """scipy.signal.TransferFunction: Uncompensated loop :math:`L(s)`."""

    compensated_loop: Any
    """scipy.signal.TransferFunction: Compensated loop :math:`C(s)L(s)`."""

    phase_margin_deg: float
    """float: Phase margin of the compensated loop [deg]."""

    gain_cross_frequency: float
    """float: Gain-crossover frequency of the compensated loop [rad/s]."""

    gain_margin: float
    """float: Gain margin of the compensated loop [-]."""

    def summary(self) -> dict:
        """Achieved margins, crossover and the compensator's high-frequency gain."""
        num, den = np.atleast_1d(np.squeeze(self.compensator.num)), np.atleast_1d(np.squeeze(self.compensator.den))
        return {
            "phase_margin_deg": float(self.phase_margin_deg),
            "gain_cross_frequency": float(self.gain_cross_frequency),
            "gain_margin": float(self.gain_margin),
            "high_frequency_gain": float(num[0] / den[0]),
        }

    def plot(self, w=None, axes=None):
        """Bode plots of the loop before and after compensation; return ``(fig, axs)``."""
        from engineeringkit.control_systems.visualizers.plots import plot_bode

        axes = plot_bode(self.loop, w, axes=axes, label="uncompensated")
        plot_bode(self.compensated_loop, w, axes=axes, label=f"{self.kind} compensated")
        return axes[0].figure, axes


@dataclass(frozen=True)
class PhasePortraitResult:
    r"""Trajectory of a second-order system :math:`\ddot x = f(x, \dot x)` in its phase plane."""

    t: np.ndarray
    """ndarray: Time samples [s]."""

    x: np.ndarray
    """ndarray: Position at each sample, in the model's units."""

    v: np.ndarray
    """ndarray: Velocity at each sample, in position units per second."""

    acceleration: Callable
    """callable: ``acceleration(x, v)``, the right-hand side :math:`f`."""

    def summary(self) -> dict:
        """Final state and the largest position excursion."""
        return {"final_x": float(self.x[-1]), "final_v": float(self.v[-1]), "max_abs_x": float(np.max(np.abs(self.x)))}

    def plot(self, ax=None, field=True):
        """Plot the trajectory over the direction field; return ``(fig, ax)``."""
        from engineeringkit.control_systems.visualizers.plots import plot_phase_portrait

        ax = plot_phase_portrait(self.x, self.v, self.acceleration if field else None, ax=ax)
        return ax.figure, ax
