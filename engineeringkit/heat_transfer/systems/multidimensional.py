r"""Two-dimensional steady conduction: finite differences and conduction shape factors.

Laplace's equation :math:`\nabla^2 T = 0` has closed-form solutions only for
simple shapes. The five-point finite-difference stencil turns it into a
sparse linear system for any rectangle, and for common buried and
eccentric geometries the conduction shape factor :math:`S` packs the
solution into :math:`\dot Q = Sk\,\Delta T`.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.heat_transfer.core.base import Conduction2DResult

__all__ = [
    "conduction_2d",
    "shape_factor_buried_sphere",
    "shape_factor_buried_cylinder",
    "shape_factor_eccentric_cylinders",
    "shape_factor_parallel_cylinders",
]


def conduction_2d(width, height, nx, ny, T_left, T_right, T_bottom, T_top, k=1.0):
    r"""Steady two-dimensional conduction in a rectangle with fixed edge temperatures, by finite differences.

    The interior nodes satisfy the five-point form of Laplace's equation,

    .. math::

        \frac{T_{i+1,j} - 2T_{i,j} + T_{i-1,j}}{\Delta x^2}
        + \frac{T_{i,j+1} - 2T_{i,j} + T_{i,j-1}}{\Delta y^2} = 0,

    solved as one sparse linear system. The heat entering through each edge
    is summed from the conduction between each edge node and its interior
    neighbour, so the four edge heat rates add up to zero to round-off.

    Parameters
    ----------
    width, height : float
        Sides of the rectangle along :math:`x` and :math:`y` [m], positive.
    nx, ny : int
        Number of nodes along :math:`x` and :math:`y`, edges included; at least three.
    T_left, T_right, T_bottom, T_top : float
        Temperatures of the edges :math:`x = 0`, :math:`x = ` ``width``,
        :math:`y = 0` and :math:`y = ` ``height`` [K], nonnegative.
    k : float, optional
        Thermal conductivity [W/(m K)], positive.

    Returns
    -------
    Conduction2DResult
        Node coordinates, the temperature field (rows along :math:`y`) and
        the heat rate per unit depth entering through each edge.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Constant conductivity, no generation. The
    corner nodes belong to two edges and take the mean of their temperatures;
    they enter no equation. The error falls as :math:`\Delta x^2`.

    Examples
    --------
    >>> result = conduction_2d(1.0, 1.0, 41, 41, 300.0, 300.0, 300.0, 400.0)
    >>> round(float(result.T[20, 20]), 2)
    325.0
    """
    _positive(width=width, height=height, k=k)
    _integer("nx", nx, 3)
    _integer("ny", ny, 3)
    _finite(T_left=T_left, T_right=T_right, T_bottom=T_bottom, T_top=T_top)
    _nonnegative(T_left=T_left, T_right=T_right, T_bottom=T_bottom, T_top=T_top)
    dx, dy = width / (nx - 1), height / (ny - 1)
    T = np.zeros((ny, nx))
    T[:, 0], T[:, -1], T[0, :], T[-1, :] = T_left, T_right, T_bottom, T_top
    for row, col in ((0, 0), (0, -1), (-1, 0), (-1, -1)):
        T[row, col] = 0.5 * ((T_left if col == 0 else T_right) + (T_bottom if row == 0 else T_top))

    mx, my = nx - 2, ny - 2
    ax, ay = 1 / dx**2, 1 / dy**2
    Dx = sparse.diags([np.ones(mx - 1), -2 * np.ones(mx), np.ones(mx - 1)], [-1, 0, 1])
    Dy = sparse.diags([np.ones(my - 1), -2 * np.ones(my), np.ones(my - 1)], [-1, 0, 1])
    A = ax * sparse.kron(sparse.identity(my), Dx) + ay * sparse.kron(Dy, sparse.identity(mx))
    rhs = np.zeros((my, mx))
    rhs[:, 0] -= ax * T[1:-1, 0]
    rhs[:, -1] -= ax * T[1:-1, -1]
    rhs[0, :] -= ay * T[0, 1:-1]
    rhs[-1, :] -= ay * T[-1, 1:-1]
    T[1:-1, 1:-1] = spsolve(A.tocsc(), rhs.ravel()).reshape(my, mx)

    heat_rates = {
        "left": float(k * dy / dx * np.sum(T[1:-1, 0] - T[1:-1, 1])),
        "right": float(k * dy / dx * np.sum(T[1:-1, -1] - T[1:-1, -2])),
        "bottom": float(k * dx / dy * np.sum(T[0, 1:-1] - T[1, 1:-1])),
        "top": float(k * dx / dy * np.sum(T[-1, 1:-1] - T[-2, 1:-1])),
    }
    return Conduction2DResult(np.linspace(0.0, width, nx), np.linspace(0.0, height, ny), T, heat_rates)


def shape_factor_buried_sphere(D, z):
    r"""Conduction shape factor of an isothermal sphere buried in a semi-infinite medium, :math:`S = 2\pi D/(1 - D/4z)`.

    Parameters
    ----------
    D : float or array_like
        Sphere diameter [m], positive.
    z : float or array_like
        Depth of the centre below the isothermal surface [m], greater than ``D/2``.

    Returns
    -------
    float or ndarray
        Shape factor [m]: :math:`\dot Q = Sk(T_1 - T_2)`. Tends to the
        isolated sphere's :math:`2\pi D` as the depth grows.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. From the method of images; exact as
    :math:`z \to \infty`, approximate near the surface.

    Examples
    --------
    >>> round(float(shape_factor_buried_sphere(1.0, 1.0)), 4)
    8.3776
    """
    _positive(D=D, z=z)
    D, z = np.asarray(D, dtype=float), np.asarray(z, dtype=float)
    if np.any(z <= D / 2):
        raise ValueError("z must exceed D/2")
    return 2 * np.pi * D / (1 - D / (4 * z))


def shape_factor_buried_cylinder(D, z, L=1.0):
    r"""Conduction shape factor of a long isothermal cylinder buried in a semi-infinite medium, :math:`S = 2\pi L/\cosh^{-1}(2z/D)`.

    Parameters
    ----------
    D : float or array_like
        Cylinder diameter [m], positive.
    z : float or array_like
        Depth of the axis below the isothermal surface [m], greater than ``D/2``.
    L : float or array_like, optional
        Length [m], positive and much larger than ``D``; the default 1 gives
        the shape factor per metre.

    Returns
    -------
    float or ndarray
        Shape factor [m]; for :math:`z \gg D` it tends to :math:`2\pi L/\ln(4z/D)`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Exact in two dimensions (method of images).

    Examples
    --------
    >>> round(float(shape_factor_buried_cylinder(0.5, 1.0)), 4)
    3.045
    """
    _positive(D=D, z=z, L=L)
    D, z = np.asarray(D, dtype=float), np.asarray(z, dtype=float)
    if np.any(z <= D / 2):
        raise ValueError("z must exceed D/2")
    return 2 * np.pi * np.asarray(L, dtype=float) / np.arccosh(2 * z / D)


def shape_factor_eccentric_cylinders(d, D, z, L=1.0):
    r"""Conduction shape factor between a cylinder and an enclosing eccentric cylinder.

    .. math::

        S = \frac{2\pi L}{\cosh^{-1}\left[(D^2 + d^2 - 4z^2)/(2Dd)\right]}

    Parameters
    ----------
    d : float or array_like
        Inner diameter [m], positive.
    D : float or array_like
        Outer diameter [m], greater than ``d + 2z``.
    z : float or array_like
        Offset between the two axes [m], nonnegative.
    L : float or array_like, optional
        Length [m], positive; the default 1 gives the shape factor per metre.

    Returns
    -------
    float or ndarray
        Shape factor [m]; at :math:`z = 0` it is the concentric
        :math:`2\pi L/\ln(D/d)` of :func:`cylinder_resistance`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Exact (bipolar coordinates).

    Examples
    --------
    >>> round(float(shape_factor_eccentric_cylinders(0.1, 0.3, 0.05)), 4)
    6.5285
    """
    _positive(d=d, D=D, L=L)
    _nonnegative(z=z)
    d, D, z = np.asarray(d, dtype=float), np.asarray(D, dtype=float), np.asarray(z, dtype=float)
    if np.any(D <= d + 2 * z):
        raise ValueError("D must exceed d + 2z (the inner cylinder must lie inside)")
    return 2 * np.pi * np.asarray(L, dtype=float) / np.arccosh((D**2 + d**2 - 4 * z**2) / (2 * D * d))


def shape_factor_parallel_cylinders(D1, D2, w, L=1.0):
    r"""Conduction shape factor between two long parallel cylinders in an infinite medium.

    .. math::

        S = \frac{2\pi L}{\cosh^{-1}\left[(4w^2 - D_1^2 - D_2^2)/(2D_1D_2)\right]}

    Parameters
    ----------
    D1, D2 : float or array_like
        Diameters [m], positive.
    w : float or array_like
        Distance between the axes [m], greater than ``(D1 + D2)/2``.
    L : float or array_like, optional
        Length [m], positive; the default 1 gives the shape factor per metre.

    Returns
    -------
    float or ndarray
        Shape factor [m]. For equal cylinders the midplane is isothermal, so
        it is half that of one cylinder buried at depth :math:`w/2`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Exact in two dimensions.

    Examples
    --------
    >>> round(float(shape_factor_parallel_cylinders(0.5, 0.5, 2.0)), 4)
    1.5225
    """
    _positive(D1=D1, D2=D2, w=w, L=L)
    D1, D2, w = np.asarray(D1, dtype=float), np.asarray(D2, dtype=float), np.asarray(w, dtype=float)
    if np.any(2 * w <= D1 + D2):
        raise ValueError("w must exceed (D1 + D2)/2")
    return 2 * np.pi * np.asarray(L, dtype=float) / np.arccosh((4 * w**2 - D1**2 - D2**2) / (2 * D1 * D2))
