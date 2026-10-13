"""Result containers for engineeringkit.control_systems.

Most control models return SciPy LTI objects, arrays or dictionaries. A
simulated closed-loop response is returned as the frozen dataclass below so
it can be summarised and plotted consistently.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["ResponseResult"]


@dataclass(frozen=True)
class ResponseResult:
    """Time response of a simulated control loop."""

    t: np.ndarray
    """ndarray: Time samples [s]."""

    y: np.ndarray
    """ndarray: Output at each sample, in output units per unit reference."""

    def summary(self) -> dict:
        """Final value and peak of the response."""
        return {"final": float(self.y[-1]), "peak": float(np.max(self.y))}

    def plot(self, ax=None):
        """Plot the response against time; return ``(fig, ax)``."""
        from engineeringkit.control_systems.visualizers.plots import plot_response

        ax = plot_response(self.t, self.y, ax=ax)
        return ax.figure, ax
