r"""Rigid-jointed frames: Maney's slope-deflection, Hardy Cross's moment distribution and the plane-frame stiffness method.

Slope-deflection writes each member-end moment in terms of the joint
rotations, :math:`M_{ij} = (2EI/L)(2\theta_i + \theta_j - 3\psi) + M^F_{ij}`, and
solves the joint equilibrium equations for the rotations at once.

Moment distribution locks every joint, applies the fixed-end moments, then
releases the joints in turn: each unbalanced joint moment is shared among
the members in proportion to their stiffness :math:`4EI/L` and half of each
share is carried over to the member's far end. The direct stiffness method
solves the same problem in one step, with three degrees of freedom per node
``(u, v, theta)``, and also handles sway and axial deformation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.structural_engineering.core.base import FrameResult, MomentDistributionResult, SlopeDeflectionResult
from engineeringkit.structural_engineering.utils.stiffness import (
    fixed_dof_vector,
    frame_rotation_matrix,
    frame_stiffness_matrix,
    mesh_geometry,
    solve_reduced,
)

__all__ = ["slope_deflection", "moment_distribution", "PlanarFrame"]


def _braced_members(members, fixed_joints):
    """Validate ``(i, j, length, EI, w)`` members; return end joints, properties, joints and the fixed set."""
    if not len(members):
        raise ValueError("at least one member is required")
    ends = np.zeros((len(members), 2), int)
    props = np.zeros((len(members), 3))
    for k, member in enumerate(members):
        if len(member) != 5:
            raise ValueError("each member must be (i, j, length, EI, w)")
        i, j, length, EI, w = member
        _integer("member joint", i, 0)
        _integer("member joint", j, 0)
        if i == j:
            raise ValueError("a member must join two different joints")
        _positive(length=length, EI=EI)
        _finite(w=w)
        ends[k] = i, j
        props[k] = length, EI, w
    joints = np.unique(ends)
    fixed = set()
    for joint in fixed_joints:
        _integer("fixed joint", joint, 0)
        if joint not in joints:
            raise ValueError(f"fixed joint {joint} is not connected to any member")
        fixed.add(int(joint))
    return ends, props, joints, fixed


def slope_deflection(members, fixed_joints=(), chord_rotations=None):
    r"""Member-end moments of a continuous beam or braced frame by Maney's slope-deflection method.

    With clockwise moments and rotations positive, the moment at end
    :math:`i` of a member :math:`ij` is

    .. math::

        M_{ij} = \frac{2EI}{L}\left(2\theta_i + \theta_j - 3\psi\right) + M^F_{ij},

    where :math:`\psi` is the clockwise rotation of the chord (a support
    settlement :math:`\Delta` of the far end gives :math:`\psi = \Delta/L`)
    and :math:`M^F = \mp wL^2/12` the fixed-end moments of a uniform load.
    The moments at each free joint must sum to zero; these equations are
    solved for the unknown rotations, and the rotations give the moments.
    It is the same problem :func:`moment_distribution` solves by iteration.

    Parameters
    ----------
    members : sequence of (int, int, float, float, float)
        ``(i, j, length, EI, w)`` for each member, as in
        :func:`moment_distribution`.
    fixed_joints : sequence of int, optional
        Joints restrained against rotation.
    chord_rotations : array_like, optional
        Known chord rotation :math:`\psi` of each member [rad], clockwise
        positive; zero by default (no settlement, no sway).

    Returns
    -------
    SlopeDeflectionResult
        End moments [N m] and joint rotations [rad], clockwise positive.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. Joint translations other than the given chord
    rotations are not solved for: sway frames need one more equilibrium
    equation per storey, as in :class:`PlanarFrame`.

    References
    ----------
    G. A. Maney, "Secondary stresses in rigid frames," *Studies in
    Engineering* (Minneapolis: University of Minnesota, 1915).
    """
    ends, props, joints, fixed = _braced_members(members, fixed_joints)
    L, EI, w = props.T
    psi = np.zeros(len(members)) if chord_rotations is None else _vector("chord_rotations", chord_rotations)
    if psi.size != len(members):
        raise ValueError("chord_rotations must have one entry per member")
    free = [int(joint) for joint in joints if joint not in fixed]
    index = {joint: k for k, joint in enumerate(free)}
    k = 2 * EI / L
    # end moments with every joint locked: fixed-end moments plus the chord rotation
    locked = np.column_stack([-w * L**2 / 12 - 3 * k * psi, w * L**2 / 12 - 3 * k * psi])
    K = np.zeros((len(free), len(free)))
    rhs = np.zeros(len(free))
    for m, (i, j) in enumerate(ends):
        for near, far, end in ((i, j, 0), (j, i, 1)):
            if near in index:
                K[index[near], index[near]] += 2 * k[m]
                if far in index:
                    K[index[near], index[far]] += k[m]
                rhs[index[near]] -= locked[m, end]
    theta = np.zeros(int(joints.max()) + 1)
    theta[free] = np.linalg.solve(K, rhs)
    ti, tj = theta[ends[:, 0]], theta[ends[:, 1]]
    end_moments = locked + np.column_stack([k * (2 * ti + tj), k * (2 * tj + ti)])
    return SlopeDeflectionResult(end_moments, theta)


def moment_distribution(members, fixed_joints=(), *, tol=1e-10, max_iterations=10_000):
    r"""Member-end moments of a continuous beam or braced frame by Hardy Cross's method.

    Every joint is assumed restrained against translation (no sway). A
    joint listed in ``fixed_joints`` is also restrained against rotation;
    every other joint rotates freely, so a joint with a single member is a
    pinned end.

    The fixed-end moments of a member carrying a uniform load :math:`w` are
    :math:`\mp wL^2/12`. Each cycle computes the unbalanced moment at every
    free joint, distributes its negative to the member ends in proportion to
    :math:`4EI/L`, and carries half of each distributed moment to the far
    end. Cycles stop when the largest unbalance falls below ``tol`` times
    the largest fixed-end moment.

    Parameters
    ----------
    members : sequence of (int, int, float, float, float)
        ``(i, j, length, EI, w)`` for each member: joint indices, length [m],
        flexural rigidity [N m^2] and uniform transverse load [N/m]. For a
        horizontal member drawn from ``i`` on the left to ``j`` on the
        right, positive ``w`` acts downward.
    fixed_joints : sequence of int, optional
        Joints restrained against rotation (built-in supports).
    tol : float, optional
        Relative convergence tolerance [-].
    max_iterations : int, optional
        Maximum number of cycles.

    Returns
    -------
    MomentDistributionResult
        End moments [N m], clockwise positive on the member end, with the
        number of cycles and the unbalance history.

    Raises
    ------
    RuntimeError
        If the unbalance has not converged after ``max_iterations`` cycles.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.

    References
    ----------
    H. Cross, "Analysis of continuous frames by distributing fixed-end
    moments," *Proc. ASCE* 56 (1930) 919-928.
    """
    ends, props, joints, fixed = _braced_members(members, fixed_joints)
    _positive(tol=tol)
    _integer("max_iterations", max_iterations, 1)
    L, EI, w = props.T
    stiffness = 4 * EI / L
    moments = np.column_stack([-w * L**2 / 12, w * L**2 / 12])
    free = [joint for joint in joints if joint not in fixed]
    end_stiffness = np.column_stack([stiffness, stiffness])
    factors = np.zeros_like(moments)
    for joint in free:
        at = ends == joint
        factors[at] = end_stiffness[at] / end_stiffness[at].sum()
    scale = np.max(np.abs(moments))
    history = []
    for iteration in range(max_iterations + 1):
        unbalanced = {joint: float(np.sum(moments[ends == joint])) for joint in free}
        worst = max((abs(value) for value in unbalanced.values()), default=0.0)
        history.append(worst)
        if worst <= tol * scale:
            return MomentDistributionResult(moments, iteration, np.array(history))
        if iteration == max_iterations:
            break
        distributed = np.zeros_like(moments)
        for joint, value in unbalanced.items():
            at = ends == joint
            distributed[at] = -factors[at] * value
        moments = moments + distributed + 0.5 * distributed[:, ::-1]
    raise RuntimeError(f"moment distribution did not converge in {max_iterations} cycles")


@dataclass(frozen=True)
class PlanarFrame:
    """Rigid-jointed plane frame solved by the direct stiffness method.

    Members are Euler-Bernoulli beam-columns with axial and bending
    stiffness, joined rigidly at the nodes.

    Parameters
    ----------
    nodes : array_like, shape (n, 2)
        Node coordinates [m].
    members : sequence of (int, int)
        Zero-based node indices of each member; the member's local x axis
        runs from the first node to the second.
    young_modulus : float
        Young's modulus [Pa], positive.
    area : float or array_like
        Cross-sectional area [m^2], one value or one per member, positive.
    second_moment : float or array_like
        Second moment of area [m^4], one value or one per member, positive.
    """

    nodes: np.ndarray
    members: tuple
    young_modulus: float
    area: float | np.ndarray
    second_moment: float | np.ndarray

    def solve(self, loads, fixed_dofs, member_loads=None, cond_limit=1e12):
        r"""Solve for displacements, reactions and member end forces.

        Parameters
        ----------
        loads : array_like
            Nodal loads ``(Fx [N], Fy [N], Mz [N m])`` per node, ``3n`` values
            or shape ``(n, 3)``; moments counterclockwise positive.
        fixed_dofs : array_like of int
            Restrained degrees of freedom, ``3 * node + (0, 1, 2)`` for
            ``(u, v, theta)``.
        member_loads : float or array_like, optional
            Uniform transverse load on each member [N/m], acting along the
            member's local -y axis (downward on a member drawn left to
            right). Converted to equivalent nodal loads
            :math:`(wL/2, wL^2/12)`.
        cond_limit : float, optional
            A :class:`RuntimeWarning` is issued if :math:`\operatorname{cond}(K_{ff})`
            exceeds this value.

        Returns
        -------
        FrameResult

        Raises
        ------
        ValueError
            If :math:`K_{ff}` is singular (a mechanism, or a frame not
            restrained against rigid-body motion).

        References
        ----------
        R. K. Livesley, "The application of an electronic digital computer to
        some problems of structural analysis," *Struct. Eng.* 34 (1956) 1-12.
        """
        xy, geometry = mesh_geometry(self.nodes, self.members)
        m = len(self.members)
        _positive(E=self.young_modulus, area=self.area, second_moment=self.second_moment)
        A = np.broadcast_to(np.asarray(self.area, float), (m,))
        I = np.broadcast_to(np.asarray(self.second_moment, float), (m,))
        ndof = 3 * len(xy)
        f = np.asarray(loads, float).reshape(-1)
        if f.size != ndof:
            raise ValueError("loads must contain three components per node")
        _finite(loads=f)
        w = np.zeros(m) if member_loads is None else np.asarray(member_loads, float)
        _finite(member_loads=w)
        w = np.broadcast_to(w, (m,))
        fixed = fixed_dof_vector(fixed_dofs, ndof)
        K = np.zeros((ndof, ndof))
        total = f.copy()
        elements = []
        for k, ((i, j), (L, c, s)) in enumerate(zip(self.members, geometry)):
            dofs = [3 * i, 3 * i + 1, 3 * i + 2, 3 * j, 3 * j + 1, 3 * j + 2]
            T = frame_rotation_matrix(c, s)
            k_local = frame_stiffness_matrix(self.young_modulus, A[k], I[k], L)
            equivalent = w[k] * np.array([0.0, -L / 2, -(L**2) / 12, 0.0, -L / 2, L**2 / 12])
            K[np.ix_(dofs, dofs)] += T.T @ k_local @ T
            total[dofs] += T.T @ equivalent
            elements.append((dofs, T, k_local, equivalent))
        u, cond = solve_reduced(K, total, fixed, cond_limit, name="frame")
        reactions = K @ u - total
        reactions[np.setdiff1d(np.arange(ndof), fixed)] = 0.0
        end_forces = np.array([k_local @ T @ u[dofs] - equivalent for dofs, T, k_local, equivalent in elements])
        return FrameResult(xy, tuple(self.members), u, reactions, end_forces, cond)
