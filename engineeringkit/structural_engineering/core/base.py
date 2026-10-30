"""Result containers for engineeringkit.structural_engineering.

Closed-form member formulas return plain floats or arrays. Section
properties and the solvers (beams, trusses, the method of joints, moment
distribution, slope-deflection, plane frames and three-hinged arches) return
the frozen dataclasses below, each with ``summary()`` and, where a picture
helps, ``plot()``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "BeamResult",
    "TrussResult",
    "JointMethodResult",
    "MomentDistributionResult",
    "FrameResult",
    "SectionProperties",
    "SlopeDeflectionResult",
    "ArchResult",
]


@dataclass(frozen=True)
class BeamResult:
    """Shear, moment and deflection along a beam.

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

    reactions: tuple
    """tuple of float: Support reactions [N], upward positive, left to right
    (left and right for :class:`SimplySupportedBeam`)."""

    support_moments: tuple = ()
    """tuple of float: Reaction moments of the fixed supports [N m], clockwise
    positive, left to right; empty when no support is fixed."""

    slope: np.ndarray | None = None
    """ndarray or None: Slope :math:`dy/dx` at each station [-], counterclockwise
    positive; ``None`` when the solver does not return it."""

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


@dataclass(frozen=True)
class SectionProperties:
    r"""Area, centroid and second moments of a cross-section built from rectangles.

    The section lies in the :math:`y`-:math:`z` plane with :math:`z`
    horizontal and :math:`y` vertical (upward), as in Hibbeler and Gere:
    :math:`I_z` is the second moment about the horizontal centroidal axis
    (the strong axis of an I-beam), :math:`I_y` about the vertical one, and
    :math:`I_{yz} = \int yz\,dA`. Each rectangle carries a weight: 1 for
    solid material, -1 for a hole and the modular ratio
    :math:`n = E/E_\mathrm{ref}` for a transformed section, in which case
    every property is that of the transformed section in the reference
    material.
    """

    area: float
    """float: (Transformed) area [m^2]."""

    y_centroid: float
    """float: Vertical position of the centroid [m]."""

    z_centroid: float
    """float: Horizontal position of the centroid [m]."""

    Iz: float
    """float: Second moment about the horizontal centroidal axis [m^4]."""

    Iy: float
    """float: Second moment about the vertical centroidal axis [m^4]."""

    Iyz: float
    """float: Product of inertia about the centroidal axes [m^4]."""

    rectangles: tuple
    """tuple: The parts as ``(b, h, z, y, weight)`` [m, m, m, m, -], with ``(z, y)`` each part's centre."""

    def principal_moments(self) -> tuple[float, float, float]:
        r"""Principal second moments and the direction of the major axis.

        Returns
        -------
        tuple of float
            ``(I_max [m^4], I_min [m^4], angle [rad])``; the major principal
            axis is the :math:`z` axis turned counterclockwise by ``angle``,
            with :math:`\tan 2\theta = -2I_{yz}/(I_z - I_y)`.
        """
        mean, half = (self.Iz + self.Iy) / 2, (self.Iz - self.Iy) / 2
        radius = float(np.hypot(half, self.Iyz))
        return mean + radius, mean - radius, float(0.5 * np.arctan2(-2 * self.Iyz, self.Iz - self.Iy))

    def first_moment(self, y: float) -> float:
        r"""First moment :math:`Q` of the part above height ``y`` about the centroidal :math:`z` axis.

        Parameters
        ----------
        y : float
            Height of the cut [m], in the coordinates of the rectangles.

        Returns
        -------
        float
            :math:`Q = \int_{A'} (y - \bar y)\,dA` [m^3], weighted like the
            area; the part below the cut has the opposite value.
        """
        Q = 0.0
        for b, h, _, yc, weight in self.rectangles:
            top, bottom = yc + h / 2, max(yc - h / 2, y)
            if top > bottom:
                Q += weight * b * (top - bottom) * ((top + bottom) / 2 - self.y_centroid)
        return float(Q)

    def width_at(self, y: float) -> float:
        """Width of material cut by the horizontal line at height ``y``.

        Holes subtract their width and modular ratios are ignored, so this
        is the true thickness :math:`t` in Jourawski's formula.

        Parameters
        ----------
        y : float
            Height of the cut [m].

        Returns
        -------
        float
            Width [m].
        """
        return float(sum(np.sign(weight) * b for b, h, _, yc, weight in self.rectangles if yc - h / 2 <= y <= yc + h / 2))

    def summary(self) -> dict:
        """Area [m^2], centroid [m] and centroidal second moments [m^4]."""
        return {
            "area": self.area,
            "y_centroid": self.y_centroid,
            "z_centroid": self.z_centroid,
            "Iz": self.Iz,
            "Iy": self.Iy,
            "Iyz": self.Iyz,
        }

    def plot(self, ax=None):
        """Draw the parts (holes white, transformed parts orange) and the centroid; return ``(fig, ax)``.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure is created if omitted.

        Returns
        -------
        tuple
            ``(fig, ax)``.
        """
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle

        if ax is None:
            _, ax = plt.subplots()
        for b, h, z, y, weight in sorted(self.rectangles, key=lambda part: part[4] < 0):
            color = "white" if weight < 0 else "C0" if weight == 1 else "C1"
            ax.add_patch(Rectangle((z - b / 2, y - h / 2), b, h, facecolor=color, edgecolor="k"))
        ax.plot(self.z_centroid, self.y_centroid, "+", color="C3", ms=12)
        ax.autoscale_view()
        ax.set_aspect("equal")
        ax.set(xlabel="z [m]", ylabel="y [m]")
        return ax.figure, ax


@dataclass(frozen=True)
class SlopeDeflectionResult:
    """Output of :func:`engineeringkit.structural_engineering.slope_deflection`.

    Moments and rotations are clockwise positive, as in Maney's method and
    Hardy Cross's moment distribution, so ``end_moments`` compare directly
    with :class:`MomentDistributionResult`.
    """

    end_moments: np.ndarray
    """ndarray, shape (m, 2): Moments at the ``i`` and ``j`` ends of each member [N m], clockwise positive."""

    rotations: np.ndarray
    """ndarray: Rotation of each joint [rad], clockwise positive (zero at fixed joints)."""

    def summary(self) -> dict:
        """Largest end moment [N m] and joint rotation [rad]."""
        return {
            "max_abs_end_moment": float(np.max(np.abs(self.end_moments))),
            "max_abs_rotation": float(np.max(np.abs(self.rotations))),
        }


@dataclass(frozen=True)
class ArchResult:
    """Shape, bending moment and axial force of a three-hinged arch.

    Loads act downward when positive and sagging moment (intrados in
    tension) is positive.
    """

    x: np.ndarray
    """ndarray: Stations along the span [m]."""

    y: np.ndarray
    """ndarray: Height of the arch axis above the springings [m]."""

    moment: np.ndarray
    """ndarray: Bending moment at each station [N m], sagging positive."""

    axial_force: np.ndarray
    """ndarray: Axial force along the arch axis [N], tension positive."""

    reactions: tuple
    """tuple of float: Left and right vertical reactions [N], upward positive."""

    thrust: float
    """float: Horizontal thrust on each springing [N], pushing the supports outward."""

    def summary(self) -> dict:
        """Thrust [N], vertical reactions [N], largest moment [N m] and largest compression [N]."""
        return {
            "thrust": self.thrust,
            "reactions": self.reactions,
            "max_abs_moment": float(np.max(np.abs(self.moment))),
            "max_compression": float(np.max(-self.axial_force)),
        }

    def plot(self):
        """Draw the arch with its three hinges above its bending-moment diagram; return ``(fig, axes)``."""
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 1, sharex=True)
        axes[0].plot(self.x, self.y, color="k")
        span = self.x[-1] - self.x[0]
        axes[0].plot([self.x[0], self.x[0] + span / 2, self.x[-1]], [0, np.max(self.y), 0], "o", mfc="white", color="k")
        axes[0].set_aspect("equal")
        axes[0].set_ylabel("y [m]")
        axes[1].plot(self.x, self.moment)
        axes[1].axhline(0, color="gray", lw=0.5)
        axes[1].set(xlabel="x [m]", ylabel="moment [N m]")
        return fig, axes
