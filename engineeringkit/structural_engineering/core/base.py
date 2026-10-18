"""Result containers for engineeringkit.structural_engineering.

Closed-form member and section formulas return plain floats or arrays. The
solvers (beams, trusses, the method of joints, moment distribution and plane
frames) return the frozen dataclasses below, each with ``summary()`` and,
where a picture helps, ``plot()``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["BeamResult", "TrussResult", "JointMethodResult", "MomentDistributionResult", "FrameResult"]


@dataclass(frozen=True)
class BeamResult:
    """Shear, moment and deflection along a simply supported beam.

    Positive loads act downward, reactions upward, sagging moment is
    positive and downward deflection is negative.
    """

    x: np.ndarray
    """ndarray: Stations along the span [m]."""

    shear: np.ndarray
    """ndarray: Shear force at each station [N]."""

    moment: np.ndarray
    """ndarray: Bending moment at each station [N m], sagging positive."""

    deflection: np.ndarray
    """ndarray: Transverse deflection at each station [m], downward negative."""

    reactions: tuple[float, float]
    """tuple of float: Left and right support reactions [N], upward positive."""

    def summary(self) -> dict:
        """Largest absolute shear [N], moment [N m] and deflection [m], and the reactions [N]."""
        return {
            "max_abs_shear": float(np.max(abs(self.shear))),
            "max_abs_moment": float(np.max(abs(self.moment))),
            "max_abs_deflection": float(np.max(abs(self.deflection))),
            "reactions": self.reactions,
        }

    def plot(self):
        """Draw the shear, moment and deflection diagrams; return ``(fig, axes)``."""
        from engineeringkit.structural_engineering.visualizers.plots import plot_beam_diagrams

        axes = plot_beam_diagrams(self)
        return axes[0].figure, axes


@dataclass(frozen=True)
class TrussResult:
    """Displacements and member forces of a pin-jointed plane truss.

    Degrees of freedom are ordered ``x0, y0, x1, y1, ...``.
    """

    displacements: np.ndarray
    """ndarray: Nodal displacements [m], two per node."""

    reactions: np.ndarray
    """ndarray: Support forces [N], two per node (zero at free degrees of freedom)."""

    member_forces: np.ndarray
    """ndarray: Axial force in each member [N], tension positive."""

    member_stresses: np.ndarray
    """ndarray: Axial stress in each member [Pa], tension positive."""

    condition_number: float = float("nan")
    """float: 2-norm condition number of the reduced stiffness matrix ``K_ff``."""

    def summary(self) -> dict:
        """Largest absolute displacement [m] and member force [N], and ``cond(K_ff)``."""
        return {
            "max_abs_displacement": float(np.max(np.abs(self.displacements))),
            "max_abs_member_force": float(np.max(np.abs(self.member_forces))),
            "condition_number": float(self.condition_number),
        }

    def plot(self, nodes, members, scale=1.0):
        """Draw the undeformed (dashed) and deformed truss; return ``(fig, ax)``.

        Parameters
        ----------
        nodes : array_like, shape (n, 2)
            Node coordinates [m] used to solve the truss.
        members : sequence of (int, int)
            Member connectivity used to solve the truss.
        scale : float, optional
            Displacement magnification [-].
        """
        from engineeringkit.structural_engineering.visualizers.plots import plot_truss

        ax = plot_truss(nodes, members, displacements=self.displacements, scale=scale)
        return ax.figure, ax


@dataclass(frozen=True)
class JointMethodResult:
    """Output of :func:`engineeringkit.structural_engineering.method_of_joints`."""

    member_forces: np.ndarray
    """ndarray: Axial force in each member [N], tension positive."""

    reactions: np.ndarray
    """ndarray: Support forces [N], two per node (zero at free degrees of freedom)."""

    order: tuple[int, ...]
    """tuple of int: Joints in the order their equilibrium equations were used."""

    def summary(self) -> dict:
        """Largest absolute member force [N], numbers of members in tension and compression, and the joint order."""
        return {
            "max_abs_member_force": float(np.max(np.abs(self.member_forces))),
            "tension_members": int(np.count_nonzero(self.member_forces > 0)),
            "compression_members": int(np.count_nonzero(self.member_forces < 0)),
            "order": self.order,
        }


@dataclass(frozen=True)
class MomentDistributionResult:
    """Output of :func:`engineeringkit.structural_engineering.moment_distribution`.

    End moments follow Hardy Cross's convention: a moment acting on the end
    of a member is positive clockwise. For a horizontal member drawn from
    joint ``i`` on the left to joint ``j`` on the right, the sagging-positive
    bending moment is ``end_moments[k, 0]`` at ``i`` and ``-end_moments[k, 1]``
    at ``j``.
    """

    end_moments: np.ndarray
    """ndarray, shape (m, 2): Moments at the ``i`` and ``j`` ends of each member [N m], clockwise positive."""

    iterations: int
    """int: Number of balance-and-carry-over cycles performed."""

    unbalance: np.ndarray
    """ndarray: Largest unbalanced joint moment [N m] before each cycle, and after the last one."""

    def summary(self) -> dict:
        """Number of cycles, largest end moment [N m] and remaining unbalance [N m]."""
        return {
            "iterations": self.iterations,
            "max_abs_end_moment": float(np.max(np.abs(self.end_moments))),
            "final_unbalance": float(self.unbalance[-1]),
        }

    def plot(self, ax=None):
        """Plot the unbalanced moment against cycle number on a log scale; return ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        ax.semilogy(np.arange(self.unbalance.size), np.where(self.unbalance > 0, self.unbalance, np.nan), "o-")
        ax.set(xlabel="cycle", ylabel="largest unbalanced moment [N m]")
        ax.grid(True, which="both", alpha=0.3)
        return ax.figure, ax


@dataclass(frozen=True)
class FrameResult:
    """Displacements and member end forces of a rigid-jointed plane frame.

    Degrees of freedom are ordered ``u0, v0, theta0, u1, v1, theta1, ...``,
    with rotations and moments counterclockwise positive.
    """

    nodes: np.ndarray
    """ndarray, shape (n, 2): Node coordinates [m]."""

    members: tuple
    """tuple of (int, int): Member connectivity."""

    displacements: np.ndarray
    """ndarray: Nodal displacements [m] and rotations [rad], three per node."""

    reactions: np.ndarray
    """ndarray: Support forces [N] and moments [N m], three per node (zero at free degrees of freedom)."""

    member_end_forces: np.ndarray
    """ndarray, shape (m, 6): ``[N_i, V_i, M_i, N_j, V_j, M_j]`` exerted on each member by its joints, in member axes [N, N m]."""

    condition_number: float
    """float: 2-norm condition number of the reduced stiffness matrix ``K_ff``."""

    def axial_forces(self) -> np.ndarray:
        """Axial force in each member [N], tension positive."""
        return self.member_end_forces[:, 3].copy()

    def summary(self) -> dict:
        """Largest translation [m], rotation [rad], end moment [N m] and axial force [N], and ``cond(K_ff)``."""
        d = self.displacements.reshape(-1, 3)
        return {
            "max_abs_translation": float(np.max(np.hypot(d[:, 0], d[:, 1]))),
            "max_abs_rotation": float(np.max(np.abs(d[:, 2]))),
            "max_abs_end_moment": float(np.max(np.abs(self.member_end_forces[:, [2, 5]]))),
            "max_abs_axial_force": float(np.max(np.abs(self.member_end_forces[:, 3]))),
            "condition_number": float(self.condition_number),
        }

    def plot(self, scale=1.0, ax=None):
        """Draw the undeformed (dashed) and deformed frame; return ``(fig, ax)``.

        Members are drawn as cubic Hermite curves through their end
        displacements and rotations, so the bending shape is visible.

        Parameters
        ----------
        scale : float, optional
            Factor applied to the displacements so they are visible [-].
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure is created if omitted.

        Returns
        -------
        tuple
            ``(fig, ax)``.
        """
        from engineeringkit.structural_engineering.visualizers.plots import plot_frame

        ax = plot_frame(self, scale=scale, ax=ax)
        return ax.figure, ax
