"""Result containers for engineeringkit.fluid_engineering.

Most fluid models are closed forms that return plain floats or arrays.
Buckingham's dimensional analysis and Hardy Cross's network iteration return
more than one number, so they get frozen dataclasses with the usual
``summary()`` and ``plot()``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["PiGroups", "NetworkResult"]


@dataclass(frozen=True)
class PiGroups:
    """Output of :func:`engineeringkit.fluid_engineering.buckingham_pi`."""

    names: tuple
    """tuple of str: Variable names, repeating variables first."""

    bases: tuple
    """tuple of str: Base dimensions, one per row of :attr:`dimension_matrix`."""

    dimension_matrix: np.ndarray
    """ndarray of int: Exponent of each base dimension (rows) in each variable (columns)."""

    exponents: np.ndarray
    """ndarray of int: One Pi group per row, the exponent of each variable in it; shape ``(k, n)``."""

    repeating: tuple
    """tuple of str: The repeating variables, which appear in every group."""

    @property
    def rank(self) -> int:
        """int: Rank of the dimension matrix, so ``k = n - rank`` groups."""
        return len(self.names) - len(self.exponents)

    def groups(self) -> list:
        """Each Pi group as a product such as ``"dp^1 rho^-1 V^-2"``."""
        return [" ".join(f"{name}^{int(e)}" for name, e in zip(self.names, row) if e) for row in self.exponents]

    def summary(self) -> dict:
        """Number of variables, rank, number of groups and the groups themselves."""
        return {
            "variables": len(self.names),
            "rank": self.rank,
            "groups": len(self.exponents),
            "pi": self.groups(),
        }

    def plot(self, ax=None):
        """Draw the exponent table as a colored grid; return ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        limit = max(1, int(np.abs(self.exponents).max(initial=1)))
        ax.imshow(self.exponents, cmap="RdBu_r", vmin=-limit, vmax=limit)
        for (i, j), e in np.ndenumerate(self.exponents):
            if e:
                ax.text(j, i, f"{int(e)}", ha="center", va="center")
        ax.set_xticks(range(len(self.names)), self.names)
        ax.set_yticks(range(len(self.exponents)), [f"Π{i + 1}" for i in range(len(self.exponents))])
        ax.set_title("exponents of the Pi groups")
        return ax.figure, ax


@dataclass(frozen=True)
class NetworkResult:
    """Output of :func:`engineeringkit.fluid_engineering.hardy_cross`."""

    flows: np.ndarray
    """ndarray: Volume flow in each pipe [m^3/s], signed along the pipe's reference direction."""

    head_losses: np.ndarray
    """ndarray: Head loss in each pipe [m], signed like :attr:`flows`."""

    loop_residuals: np.ndarray
    """ndarray: Signed sum of head losses round each loop [m]; zero at convergence."""

    corrections: np.ndarray
    """ndarray: Largest loop correction :math:`|\\Delta Q|` of each sweep [m^3/s]."""

    def summary(self) -> dict:
        """Sweeps taken, final correction and largest loop residual."""
        return {
            "iterations": int(self.corrections.size),
            "last_correction_m3_s": float(self.corrections[-1]),
            "max_loop_residual_m": float(np.abs(self.loop_residuals).max()),
            "flows_m3_s": self.flows.tolist(),
        }

    def plot(self, ax=None):
        """Draw the loop correction against the sweep number; return ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        ax.semilogy(np.arange(1, self.corrections.size + 1), self.corrections, marker="o")
        ax.set(xlabel="sweep", ylabel="largest correction [m³/s]", title="Hardy Cross convergence")
        return ax.figure, ax
