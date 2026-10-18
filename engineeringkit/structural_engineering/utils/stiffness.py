r"""Element matrices and the reduced solve shared by the truss and frame solvers.

The direct stiffness method assembles one small matrix per member into a
global stiffness matrix :math:`K`, deletes the rows and columns of the
supported degrees of freedom and solves :math:`K_{ff}u_f = F_f`. The solve
below also reports :math:`\operatorname{cond}(K_{ff})`, since a near-mechanism
gives a matrix that is invertible in floating point but meaningless.
"""

from __future__ import annotations

import warnings

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive

__all__ = [
    "mesh_geometry",
    "fixed_dof_vector",
    "bar_stiffness_matrix",
    "frame_stiffness_matrix",
    "frame_rotation_matrix",
    "solve_reduced",
]


def mesh_geometry(nodes, members):
    """Validate a plane mesh and return its coordinates and member geometry.

    Parameters
    ----------
    nodes : array_like, shape (n, 2)
        Node coordinates [m].
    members : sequence of (int, int)
        Zero-based node indices of each member's ends.

    Returns
    -------
    xy : ndarray, shape (n, 2)
        Node coordinates [m] as floats.
    geometry : ndarray, shape (m, 3)
        Length [m] and direction cosines ``c, s`` [-] of each member, from
        its first node to its second.
    """
    xy = np.asarray(nodes, float)
    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError("nodes must have shape (n, 2)")
    _finite(nodes=xy)
    if len(xy) < 2 or not len(members):
        raise ValueError("a structure requires at least two nodes and one member")
    geometry = []
    for i, j in members:
        _integer("member node", i, 0)
        _integer("member node", j, 0)
        if i >= len(xy) or j >= len(xy):
            raise ValueError("member node index out of range")
        dx, dy = xy[j] - xy[i]
        length = float(np.hypot(dx, dy))
        if length <= 0:
            raise ValueError("members must have positive length")
        geometry.append((length, dx / length, dy / length))
    return xy, np.array(geometry)


def fixed_dof_vector(fixed_dofs, ndof):
    """Validate supported degrees of freedom and return them sorted and unique.

    Parameters
    ----------
    fixed_dofs : array_like of int
        Indices of the degrees of freedom held at zero displacement.
    ndof : int
        Total number of degrees of freedom.

    Returns
    -------
    ndarray of int
    """
    raw = np.asarray(fixed_dofs)
    if raw.size == 0:
        raw = raw.astype(int)
    if raw.ndim != 1 or not np.issubdtype(raw.dtype, np.integer):
        raise ValueError("fixed_dofs must be an integer vector")
    if np.any(raw < 0) or np.any(raw >= ndof):
        raise ValueError("fixed degree of freedom out of range")
    return np.unique(raw.astype(int))


def bar_stiffness_matrix(E, A, L, c, s):
    r"""Global 4x4 stiffness matrix of a pin-ended bar.

    Parameters
    ----------
    E : float
        Young's modulus [Pa].
    A : float
        Cross-sectional area [m^2].
    L : float
        Length [m].
    c, s : float
        Direction cosines of the bar axis [-].

    Returns
    -------
    ndarray, shape (4, 4)
        Stiffness [N/m] relating ``(x_i, y_i, x_j, y_j)`` forces to
        displacements: :math:`(EA/L)\,[\mathbf{n}\mathbf{n}^T, -\mathbf{n}\mathbf{n}^T;
        -\mathbf{n}\mathbf{n}^T, \mathbf{n}\mathbf{n}^T]` with
        :math:`\mathbf{n} = (c, s)`.
    """
    _positive(E=E, A=A, L=L)
    n = np.array([[c * c, c * s], [c * s, s * s]])
    return E * A / L * np.block([[n, -n], [-n, n]])


def frame_stiffness_matrix(E, A, I, L):
    """Local 6x6 stiffness matrix of a plane Euler-Bernoulli frame member.

    Parameters
    ----------
    E : float
        Young's modulus [Pa].
    A : float
        Cross-sectional area [m^2].
    I : float
        Second moment of area about the bending axis [m^4].
    L : float
        Length [m].

    Returns
    -------
    ndarray, shape (6, 6)
        Stiffness relating ``(N_i, V_i, M_i, N_j, V_j, M_j)`` to
        ``(u_i, v_i, theta_i, u_j, v_j, theta_j)`` in member axes, with the
        local x axis from node i to node j.
    """
    _positive(E=E, A=A, I=I, L=L)
    a = E * A / L
    b, c, d, e = 12 * E * I / L**3, 6 * E * I / L**2, 4 * E * I / L, 2 * E * I / L
    return np.array(
        [
            [a, 0, 0, -a, 0, 0],
            [0, b, c, 0, -b, c],
            [0, c, d, 0, -c, e],
            [-a, 0, 0, a, 0, 0],
            [0, -b, -c, 0, b, -c],
            [0, c, e, 0, -c, d],
        ]
    )


def frame_rotation_matrix(c, s):
    """6x6 rotation from global to member axes for a plane frame member.

    Parameters
    ----------
    c, s : float
        Direction cosines of the member axis [-].

    Returns
    -------
    ndarray, shape (6, 6)
        ``T`` such that local displacements are ``T @ global``.
    """
    r = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    z = np.zeros((3, 3))
    return np.block([[r, z], [z, r]])


def solve_reduced(K, loads, fixed, cond_limit, name="structure"):
    r"""Solve :math:`K_{ff}u_f = F_f` with a conditioning check.

    Parameters
    ----------
    K : ndarray, shape (ndof, ndof)
        Assembled global stiffness matrix.
    loads : ndarray, shape (ndof,)
        Applied nodal loads.
    fixed : ndarray of int
        Supported degrees of freedom (zero displacement).
    cond_limit : float
        A :class:`RuntimeWarning` is issued when :math:`\operatorname{cond}(K_{ff})`
        exceeds this value.
    name : str, optional
        Word used in messages ("truss", "frame").

    Returns
    -------
    u : ndarray, shape (ndof,)
        Displacements, zero at the fixed degrees of freedom.
    cond : float
        2-norm condition number of :math:`K_{ff}` (1 if every degree of
        freedom is fixed).

    Raises
    ------
    ValueError
        If :math:`K_{ff}` is singular to working precision: the structure is
        a mechanism or is not restrained against rigid-body motion.
    """
    _positive(cond_limit=cond_limit)
    free = np.setdiff1d(np.arange(len(loads)), fixed)
    u = np.zeros(len(loads))
    if not free.size:
        return u, 1.0
    Kff = K[np.ix_(free, free)]
    sv = np.linalg.svd(Kff, compute_uv=False)
    if sv[-1] <= sv[0] * len(free) * np.finfo(float).eps:
        raise ValueError(f"{name} is unstable or insufficiently constrained: K_ff is singular")
    cond = float(sv[0] / sv[-1])
    if cond > cond_limit:
        warnings.warn(
            f"{name} stiffness matrix is ill-conditioned (cond(K_ff) = {cond:.2e}); the structure is close to a mechanism",
            RuntimeWarning,
            stacklevel=3,
        )
    u[free] = np.linalg.solve(Kff, loads[free])
    return u, cond
