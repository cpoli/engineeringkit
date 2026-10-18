r"""Steady and transient conduction: Fourier's law, thermal resistances, fins and finite differences.

Fourier (1822) made heat flux proportional to the temperature gradient,
which turns steady one-dimensional conduction into a chain of resistances.
Harper and Brown (1922) solved the straight fin for NACA's air-cooled
engines, and Crank and Nicolson (1947) gave the unconditionally stable,
second-order scheme still used for the transient heat equation.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.heat_transfer.core.base import ConductionResult, WallResult

__all__ = [
    "fourier_heat_rate",
    "plane_wall_resistance",
    "cylinder_resistance",
    "sphere_resistance",
    "overall_u",
    "critical_insulation_radius_cylinder",
    "critical_insulation_radius_sphere",
    "composite_wall",
    "fin_parameter",
    "fin_efficiency_adiabatic",
    "fin_heat_rate_adiabatic",
    "conduction_1d",
]


def fourier_heat_rate(k, A, dT, L):
    r"""Steady conduction through a plane slab, :math:`\dot Q = kA\,\Delta T/L` (Fourier's law).

    Parameters
    ----------
    k : float or array_like
        Thermal conductivity [W/(m K)], positive.
    A : float or array_like
        Area normal to the heat flow [m^2], positive.
    dT : float or array_like
        Temperature difference across the slab, hot minus cold [K].
    L : float or array_like
        Slab thickness [m], positive.

    Returns
    -------
    float or ndarray
        Heat rate [W], positive from the hot face to the cold face.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    References
    ----------
    J. Fourier, *Théorie analytique de la chaleur* (Paris: Firmin Didot, 1822).

    Examples
    --------
    >>> float(fourier_heat_rate(1.0, 2.0, 30.0, 0.1))
    600.0
    """
    _finite(k=k, A=A, dT=dT, L=L)
    _positive(k=k, A=A, L=L)
    return k * A * dT / L


def plane_wall_resistance(L, k, A=1.0):
    r"""Conduction resistance of a plane wall, :math:`R = L/(kA)`.

    Parameters
    ----------
    L : float or array_like
        Thickness [m], nonnegative.
    k : float or array_like
        Thermal conductivity [W/(m K)], positive.
    A : float or array_like, optional
        Area [m^2], positive; the default 1 gives the resistance per unit area.

    Returns
    -------
    float or ndarray
        Thermal resistance [K/W].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> float(plane_wall_resistance(0.2, 0.8, 5.0))
    0.05
    """
    _finite(L=L, k=k, A=A)
    _nonnegative(L=L)
    _positive(k=k, A=A)
    return L / (k * A)


def cylinder_resistance(r1, r2, k, L=1.0):
    r"""Radial conduction resistance of a cylindrical shell, :math:`R = \ln(r_2/r_1)/(2\pi k L)`.

    Parameters
    ----------
    r1 : float or array_like
        Inner radius [m], positive.
    r2 : float or array_like
        Outer radius [m], greater than ``r1``.
    k : float
        Thermal conductivity [W/(m K)], positive.
    L : float, optional
        Length [m], positive; the default 1 gives the resistance per metre.

    Returns
    -------
    float
        Thermal resistance [K/W].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> round(float(cylinder_resistance(0.01, 0.02, 1.0)), 5)
    0.11032
    """
    _finite(r1=r1, r2=r2, k=k, L=L)
    _positive(r1=r1, r2=r2, k=k, L=L)
    r1, r2 = np.asarray(r1, dtype=float), np.asarray(r2, dtype=float)
    if np.any(r2 <= r1):
        raise ValueError("outer radius must exceed inner radius")
    return np.log(r2 / r1) / (2 * np.pi * k * L)


def sphere_resistance(r1, r2, k):
    r"""Radial conduction resistance of a spherical shell, :math:`R = (1/r_1 - 1/r_2)/(4\pi k)`.

    Parameters
    ----------
    r1 : float or array_like
        Inner radius [m], positive.
    r2 : float or array_like
        Outer radius [m], greater than ``r1``.
    k : float
        Thermal conductivity [W/(m K)], positive.

    Returns
    -------
    float
        Thermal resistance [K/W].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> round(float(sphere_resistance(1.0, 2.0, 1.0)), 6)
    0.039789
    """
    _finite(r1=r1, r2=r2, k=k)
    _positive(r1=r1, r2=r2, k=k)
    r1, r2 = np.asarray(r1, dtype=float), np.asarray(r2, dtype=float)
    if np.any(r2 <= r1):
        raise ValueError("outer radius must exceed inner radius")
    return (1 / r1 - 1 / r2) / (4 * np.pi * k)


def overall_u(*resistances):
    r"""Conductance :math:`1/\sum R_i` of thermal resistances in series (legacy name).

    Parameters
    ----------
    *resistances : float
        Full thermal resistances, each including its area [K/W], positive.

    Returns
    -------
    float
        Conductance :math:`UA` [W/K]. Divide by the reference area for the
        area-normalised coefficient :math:`U` [W/(m^2 K)].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> overall_u(0.1, 0.4)
    2.0
    """
    _positive(resistances=resistances)
    return 1 / sum(resistances)


def critical_insulation_radius_cylinder(k, h):
    r"""Critical insulation radius of a cylinder, :math:`r_{cr} = k/h`.

    Below this outer radius, adding insulation *increases* the heat loss,
    because the convective area grows faster than the conduction resistance.

    Parameters
    ----------
    k : float or array_like
        Insulation conductivity [W/(m K)], positive.
    h : float or array_like
        Outer convection coefficient [W/(m^2 K)], positive.

    Returns
    -------
    float or ndarray
        Radius of maximum heat loss [m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> critical_insulation_radius_cylinder(0.05, 10.0)
    0.005
    """
    _finite(k=k, h=h)
    _positive(k=k, h=h)
    return k / h


def critical_insulation_radius_sphere(k, h):
    r"""Critical insulation radius of a sphere, :math:`r_{cr} = 2k/h`.

    Parameters
    ----------
    k : float or array_like
        Insulation conductivity [W/(m K)], positive.
    h : float or array_like
        Outer convection coefficient [W/(m^2 K)], positive.

    Returns
    -------
    float or ndarray
        Radius of maximum heat loss [m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> critical_insulation_radius_sphere(0.05, 10.0)
    0.01
    """
    _finite(k=k, h=h)
    _positive(k=k, h=h)
    return 2 * k / h


def composite_wall(T_hot, T_cold, layers, A=1.0, h_hot=None, h_cold=None):
    r"""Steady heat flow through plane layers in series, with optional surface films.

    The resistances :math:`1/(h_{hot}A)`, :math:`L_i/(k_iA)` and
    :math:`1/(h_{cold}A)` add; the heat rate is
    :math:`\dot Q = (T_{hot} - T_{cold})/\sum R` and each node temperature
    drops by :math:`\dot Q R` from the one before.

    Parameters
    ----------
    T_hot : float
        Hot-side temperature [K]: the fluid if ``h_hot`` is given, else the surface.
    T_cold : float
        Cold-side temperature [K]: the fluid if ``h_cold`` is given, else the surface.
    layers : sequence of (float, float)
        ``(thickness [m], conductivity [W/(m K)])`` per layer, hot side first; at least one.
    A : float, optional
        Wall area [m^2], positive.
    h_hot, h_cold : float, optional
        Film coefficients [W/(m^2 K)], positive; ``None`` omits that film.

    Returns
    -------
    WallResult
        Heat rate [W] and the temperature at each node [K], hot side first.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> composite_wall(400, 300, [(0.1, 1.0), (0.2, 2.0)], A=2.0).interface_temperatures.tolist()
    [400.0, 350.0, 300.0]
    """
    _finite(T_hot=T_hot, T_cold=T_cold)
    _nonnegative(T_hot=T_hot, T_cold=T_cold)
    _positive(A=A)
    layers = list(layers)
    if not layers:
        raise ValueError("at least one wall layer is required")
    for L, k in layers:
        _positive(L=L, k=k)
    if h_hot is not None:
        _positive(h_hot=h_hot)
    if h_cold is not None:
        _positive(h_cold=h_cold)
    rs = []
    if h_hot is not None:
        rs.append(1 / (h_hot * A))
    rs += [L / (k * A) for L, k in layers]
    if h_cold is not None:
        rs.append(1 / (h_cold * A))
    q = (T_hot - T_cold) / sum(rs)
    temps = [T_hot]
    for r in rs:
        temps.append(temps[-1] - q * r)
    return WallResult(float(q), np.array(temps, dtype=float))


def fin_parameter(h, P, k, Ac):
    r"""Fin parameter :math:`m = \sqrt{hP/(kA_c)}`.

    Parameters
    ----------
    h : float or array_like
        Convection coefficient [W/(m^2 K)], nonnegative.
    P : float or array_like
        Fin perimeter [m], positive.
    k : float or array_like
        Fin conductivity [W/(m K)], positive.
    Ac : float or array_like
        Fin cross-sectional area [m^2], positive.

    Returns
    -------
    float or ndarray
        Fin parameter [1/m]; :math:`1/m` is the decay length of the excess temperature.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> float(fin_parameter(25.0, 0.04, 200.0, 5e-5))
    10.0
    """
    _finite(h=h, P=P, k=k, Ac=Ac)
    _nonnegative(h=h)
    _positive(P=P, k=k, Ac=Ac)
    return np.sqrt(h * P / (k * Ac))


def fin_efficiency_adiabatic(m, L):
    r"""Efficiency of a straight fin with an adiabatic tip, :math:`\eta = \tanh(mL)/(mL)`.

    Parameters
    ----------
    m : float or array_like
        Fin parameter [1/m], nonnegative (see :func:`fin_parameter`).
    L : float or array_like
        Fin length [m], nonnegative. For a convecting tip, Harper and
        Brown's corrected length :math:`L + A_c/P` is a good approximation.

    Returns
    -------
    float or ndarray
        Fin efficiency [-]: the heat rate over that of an isothermal fin at
        the base temperature; one in the limit :math:`mL = 0`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    References
    ----------
    D. R. Harper and W. B. Brown, "Mathematical equations for heat
    conduction in the fins of air-cooled engines," NACA Report 158 (1922).

    Examples
    --------
    >>> round(float(fin_efficiency_adiabatic(10.0, 0.1)), 4)
    0.7616
    """
    _nonnegative(m=m, L=L)
    z = np.asarray(m, dtype=float) * np.asarray(L, dtype=float)
    safe = np.where(z == 0, 1.0, z)
    return np.where(z == 0, 1.0, np.tanh(safe) / safe)[()]


def fin_heat_rate_adiabatic(h, P, k, Ac, L, Tb, Tinf):
    r"""Heat rate of a straight fin with an adiabatic tip, :math:`\sqrt{hPkA_c}\,(T_b - T_\infty)\tanh(mL)`.

    Parameters
    ----------
    h : float
        Convection coefficient [W/(m^2 K)], nonnegative.
    P : float
        Fin perimeter [m], positive.
    k : float
        Fin conductivity [W/(m K)], positive.
    Ac : float
        Cross-sectional area [m^2], positive.
    L : float or array_like
        Fin length [m], nonnegative.
    Tb : float
        Base temperature [K].
    Tinf : float
        Surrounding fluid temperature [K].

    Returns
    -------
    float or ndarray
        Heat rate leaving the base [W].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> round(float(fin_heat_rate_adiabatic(25.0, 0.04, 200.0, 5e-5, 0.1, 400.0, 300.0)), 3)
    7.616
    """
    _finite(h=h, P=P, k=k, Ac=Ac, L=L, Tb=Tb, Tinf=Tinf)
    _nonnegative(L=L, Tb=Tb, Tinf=Tinf)
    m = fin_parameter(h, P, k, Ac)
    return np.sqrt(h * P * k * Ac) * (Tb - Tinf) * np.tanh(m * L)


_THETA = {"explicit": 0.0, "crank_nicolson": 0.5, "implicit": 1.0}


def conduction_1d(T_initial, length, alpha, dt, steps, T_left=None, T_right=None, method="crank_nicolson"):
    r"""Transient one-dimensional conduction by finite differences (explicit, implicit or Crank-Nicolson).

    Solves :math:`\partial T/\partial t = \alpha\,\partial^2 T/\partial x^2`
    on evenly spaced nodes with the :math:`\theta`-scheme

    .. math::

        (I - \theta r A)\,T^{n+1} = (I + (1-\theta) r A)\,T^{n},
        \qquad r = \frac{\alpha\,\Delta t}{\Delta x^2},

    where :math:`A` is the second-difference matrix and :math:`\theta` is 0
    (explicit, stable only for :math:`r \le 1/2`), 1 (fully implicit) or 1/2
    (Crank-Nicolson, unconditionally stable and second order in time).

    Parameters
    ----------
    T_initial : array_like
        Initial temperature at each node [K], at least three nodes spanning
        ``0`` to ``length``.
    length : float
        Domain length [m], positive.
    alpha : float
        Thermal diffusivity [m^2/s], positive.
    dt : float
        Time step [s], positive.
    steps : int
        Number of time steps, at least one.
    T_left, T_right : float, optional
        Fixed surface temperature [K] at ``x = 0`` and ``x = length``;
        ``None`` (the default) makes that face adiabatic.
    method : {"crank_nicolson", "implicit", "explicit"}, optional
        Time-stepping scheme.

    Returns
    -------
    ConductionResult
        Node positions, times, the temperature history, the scheme and
        the mesh Fourier number :math:`r`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. The explicit scheme is not refused above
    :math:`r = 1/2`; it is left to grow so that the instability can be seen.

    References
    ----------
    J. Crank and P. Nicolson, "A practical method for numerical evaluation
    of solutions of partial differential equations of the heat-conduction
    type," *Proc. Cambridge Phil. Soc.* 43 (1947) 50-67.

    Examples
    --------
    >>> result = conduction_1d([300.0] * 11, 0.1, 1e-5, 1.0, 100, T_left=400.0, T_right=400.0)
    >>> round(float(result.T[-1].min()), 1)
    352.6
    """
    T = _vector("T_initial", T_initial, min_length=3)
    _nonnegative(T_initial=T)
    _positive(length=length, alpha=alpha, dt=dt)
    _integer("steps", steps, 1)
    if method not in _THETA:
        raise ValueError(f"method must be one of {sorted(_THETA)}")
    theta = _THETA[method]
    n = T.size
    dx = length / (n - 1)
    r = alpha * dt / dx**2

    # Second-difference matrix A in banded storage: rows are upper, main, lower diagonals.
    upper, main, lower = np.ones(n - 1), -2.0 * np.ones(n), np.ones(n - 1)
    for end, T_fixed, inward in ((0, T_left, upper), (n - 1, T_right, lower)):
        index = 0 if end == 0 else n - 2
        if T_fixed is None:
            inward[index] = 2.0  # mirror node across an adiabatic face
        else:
            _nonnegative(**{"T_left" if end == 0 else "T_right": T_fixed})
            main[end], inward[index] = 0.0, 0.0
            T[end] = T_fixed

    def apply_A(v):
        out = main * v
        out[:-1] += upper * v[1:]
        out[1:] += lower * v[:-1]
        return out

    ab = np.zeros((3, n))
    ab[0, 1:], ab[1], ab[2, :-1] = -theta * r * upper, 1.0 - theta * r * main, -theta * r * lower
    history = np.empty((steps + 1, n))
    history[0] = T
    for i in range(steps):
        rhs = T + (1.0 - theta) * r * apply_A(T)
        T = rhs if theta == 0 else solve_banded((1, 1), ab, rhs)
        history[i + 1] = T
    return ConductionResult(np.linspace(0.0, length, n), dt * np.arange(steps + 1), history, method, float(r))
