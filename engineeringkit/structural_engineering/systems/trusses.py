r"""Pin-jointed plane trusses: Maxwell's count, Whipple's method of joints and direct stiffness.

A plane truss with :math:`j` joints, :math:`m` members and :math:`r` support
reactions has :math:`2j` equilibrium equations in :math:`m + r` unknown
forces. Maxwell's count :math:`m + r - 2j` is negative for a mechanism, zero
for a statically determinate truss and positive for an indeterminate one.
The count is necessary but not sufficient: members can be arranged so that
some panels are over-braced and others are free to move.

Degrees of freedom are ordered ``x0, y0, x1, y1, ...``; member forces are
positive in tension.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive
from engineeringkit.structural_engineering.core.base import JointMethodResult, TrussResult
from engineeringkit.structural_engineering.utils.stiffness import bar_stiffness_matrix, fixed_dof_vector, mesh_geometry, solve_reduced

__all__ = ["truss_determinacy", "method_of_joints", "PlanarTruss"]


def truss_determinacy(n_joints, n_members, n_reactions):
    r"""Maxwell's count for a plane truss, :math:`m + r - 2j`.

    Parameters
    ----------
    n_joints : int
        Number of joints :math:`j`, at least 2.
    n_members : int
        Number of members :math:`m`, at least 1.
    n_reactions : int
        Number of restrained degrees of freedom :math:`r`, nonnegative.

    Returns
    -------
    int
        Negative: a mechanism, whatever the arrangement. Zero: statically
        determinate if the truss is also geometrically stable. Positive:
        that many redundant members or reactions.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. A zero or positive count does not prove stability;
    :class:`PlanarTruss` checks that separately through the stiffness matrix.

    References
    ----------
    J. C. Maxwell, "On the calculation of the equilibrium and stiffness of
    frames," *Phil. Mag.* 27 (1864) 294-299.
    """
    _integer("n_joints", n_joints, 2)
    _integer("n_members", n_members, 1)
    _integer("n_reactions", n_reactions, 0)
    return int(n_members + n_reactions - 2 * n_joints)


def _equilibrium_matrix(xy, geometry, members, fixed):
    """Columns: member tensions, then reactions; rows: x and y equilibrium of each joint."""
    ndof, m = 2 * len(xy), len(members)
    A = np.zeros((ndof, m + fixed.size))
    for k, ((i, j), (_, c, s)) in enumerate(zip(members, geometry)):
        A[[2 * i, 2 * i + 1], k] = c, s
        A[[2 * j, 2 * j + 1], k] = -c, -s
    A[fixed, m + np.arange(fixed.size)] = 1.0
    return A


def method_of_joints(nodes, members, loads, fixed_dofs):
    r"""Member forces of a statically determinate plane truss, joint by joint.

    This is the hand method. When the truss has three reactions, they are
    first found from the equilibrium of the whole truss. Then, repeatedly, a
    joint with at most two unknown forces is chosen and its two equilibrium
    equations, :math:`\sum F_x = 0` and :math:`\sum F_y = 0`, are solved for
    them, until every joint has been used. No stiffness enters: the forces
    of a determinate truss follow from statics alone.

    Parameters
    ----------
    nodes : array_like, shape (n, 2)
        Joint coordinates [m].
    members : sequence of (int, int)
        Zero-based joint indices of each member.
    loads : array_like
        Applied joint loads [N], ``2n`` values or shape ``(n, 2)``.
    fixed_dofs : array_like of int
        Restrained degrees of freedom; each carries one reaction.

    Returns
    -------
    JointMethodResult
        Member forces [N] (tension positive), reactions [N] at every degree
        of freedom (zero where free) and the order in which joints were used.

    Raises
    ------
    ValueError
        If Maxwell's count is not zero, or if the procedure stalls because
        every remaining joint has three or more unknowns, or two unknown
        members that are collinear (use :class:`PlanarTruss`, or the method
        of sections, for such trusses).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.

    References
    ----------
    S. Whipple, *A Work on Bridge Building* (Utica, N.Y., 1847).
    """
    xy, geometry = mesh_geometry(nodes, members)
    ndof, m = 2 * len(xy), len(members)
    f = np.asarray(loads, float).reshape(-1)
    if f.size != ndof:
        raise ValueError("loads must contain two components per node")
    _finite(loads=f)
    fixed = fixed_dof_vector(fixed_dofs, ndof)
    count = truss_determinacy(len(xy), m, fixed.size)
    if count:
        raise ValueError(f"the method of joints needs a statically determinate truss; Maxwell's count m + r - 2j = {count}")
    A = _equilibrium_matrix(xy, geometry, members, fixed)
    rhs = -f
    x = np.zeros(m + fixed.size)
    known = np.zeros(m + fixed.size, bool)
    if fixed.size == 3:
        # Whole-truss equilibrium: sum Fx, sum Fy and moments about the origin.
        rows = np.zeros((3, ndof))
        rows[0, 0::2] = 1.0
        rows[1, 1::2] = 1.0
        rows[2, 0::2], rows[2, 1::2] = -xy[:, 1], xy[:, 0]
        G = rows @ A[:, m:]
        if np.linalg.matrix_rank(G) == 3:
            x[m:] = np.linalg.solve(G, rows @ rhs)
            known[m:] = True
    order: list[int] = []
    remaining = list(range(len(xy)))
    while remaining:
        for joint in remaining:
            eq = [2 * joint, 2 * joint + 1]
            unknown = np.flatnonzero(~known & np.any(A[eq] != 0, axis=0))
            if unknown.size > 2:
                continue
            if unknown.size:
                B = A[np.ix_(eq, unknown)]
                if np.linalg.matrix_rank(B) < unknown.size:
                    continue
                x[unknown] = np.linalg.lstsq(B, rhs[eq] - A[eq] @ x, rcond=None)[0]
                known[unknown] = True
            order.append(joint)
            remaining.remove(joint)
            break
        else:
            raise ValueError(f"the method of joints stalls at joints {remaining}: each has more than two independent unknowns")
    reactions = np.zeros(ndof)
    reactions[fixed] = x[m:]
    return JointMethodResult(x[:m].copy(), reactions, tuple(order))


@dataclass(frozen=True)
class PlanarTruss:
    """Pin-jointed plane truss solved by the direct stiffness method.

    Parameters
    ----------
    nodes : array_like, shape (n, 2)
        Joint coordinates [m].
    members : sequence of (int, int)
        Zero-based joint indices of each member.
    young_modulus : float
        Young's modulus of every member [Pa], positive.
    area : float
        Cross-sectional area of every member [m^2], positive.
    """

    nodes: np.ndarray
    members: tuple
    young_modulus: float
    area: float

    def solve(self, loads, fixed_dofs, cond_limit=1e12):
        r"""Solve for displacements, reactions and member forces.

        Each member contributes :math:`(EA/L)[\mathbf{nn}^T, -\mathbf{nn}^T;
        -\mathbf{nn}^T, \mathbf{nn}^T]` to the global stiffness matrix
        :math:`K`; the free displacements solve :math:`K_{ff}u_f = F_f`.

        Parameters
        ----------
        loads : array_like
            Joint loads [N], ``2n`` values or shape ``(n, 2)``.
        fixed_dofs : array_like of int
            Restrained degrees of freedom (zero displacement).
        cond_limit : float, optional
            A :class:`RuntimeWarning` is issued if :math:`\operatorname{cond}(K_{ff})`
            exceeds this value, a sign that the truss is close to a mechanism.

        Returns
        -------
        TrussResult
            Includes ``condition_number``, :math:`\operatorname{cond}(K_{ff})`.

        Raises
        ------
        ValueError
            If Maxwell's count is negative, or :math:`K_{ff}` is singular
            (a mechanism, or a truss not restrained against rigid-body motion).

        References
        ----------
        M. J. Turner, R. W. Clough, H. C. Martin and L. J. Topp, "Stiffness
        and deflection analysis of complex structures," *J. Aeronaut. Sci.*
        23 (1956) 805-823.
        """
        xy, geometry = mesh_geometry(self.nodes, self.members)
        _positive(E=self.young_modulus, area=self.area)
        ndof = 2 * len(xy)
        f = np.asarray(loads, float).reshape(-1)
        if f.size != ndof:
            raise ValueError("loads must contain two components per node")
        _finite(loads=f)
        fixed = fixed_dof_vector(fixed_dofs, ndof)
        count = truss_determinacy(len(xy), len(self.members), fixed.size)
        if count < 0:
            raise ValueError(f"truss is unstable: Maxwell's count m + r - 2j = {count} < 0, so it is a mechanism")
        K = np.zeros((ndof, ndof))
        for (i, j), (L, c, s) in zip(self.members, geometry):
            dofs = [2 * i, 2 * i + 1, 2 * j, 2 * j + 1]
            K[np.ix_(dofs, dofs)] += bar_stiffness_matrix(self.young_modulus, self.area, L, c, s)
        u, cond = solve_reduced(K, f, fixed, cond_limit, name="truss")
        reactions = K @ u - f
        reactions[np.setdiff1d(np.arange(ndof), fixed)] = 0.0
        forces = np.array(
            [
                self.young_modulus * self.area / L * np.dot([-c, -s, c, s], u[[2 * i, 2 * i + 1, 2 * j, 2 * j + 1]])
                for (i, j), (L, c, s) in zip(self.members, geometry)
            ]
        )
        return TrussResult(u, reactions, forces, forces / self.area, cond)
