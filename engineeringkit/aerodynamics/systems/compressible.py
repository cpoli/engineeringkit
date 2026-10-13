r"""Compressible flow of a calorically perfect gas.

Sound speed, isentropic stagnation states, the Mach cone, Prandtl-Glauert
compressibility correction, the Rankine-Hugoniot normal shock, and the
area-Mach relation of the converging-diverging nozzle.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.aerodynamics.utils.gas import check_gamma, isentropic_factor

__all__ = [
    "speed_of_sound",
    "mach_number",
    "mach_angle",
    "stagnation_temperature",
    "stagnation_pressure",
    "prandtl_glauert_cp",
    "normal_shock_mach2",
    "normal_shock_pressure_ratio",
    "area_mach_ratio",
]


def speed_of_sound(gamma, R, T):
    r"""Speed of sound in a perfect gas, :math:`a = \sqrt{\gamma R T}` (Laplace).

    Parameters
    ----------
    gamma : float or array_like
        Ratio of specific heats [-], greater than one.
    R : float or array_like
        Specific gas constant [J/(kg K)], positive (287.05 for air).
    T : float or array_like
        Absolute temperature [K], positive.

    Returns
    -------
    float or ndarray
        Speed of sound [m/s].

    Notes
    -----
    Newton's isothermal estimate is :math:`\sqrt{RT}`; Laplace's adiabatic
    correction multiplies it by :math:`\sqrt{\gamma}`. Inputs are raw SI
    numbers; see :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(speed_of_sound(1.4, 287.05, 288.15)), 2)
    340.29
    """
    check_gamma(gamma)
    _positive(R=R, T=T)
    return np.sqrt(gamma * R * np.asarray(T))


def mach_number(V, a):
    r"""Mach number :math:`M = V/a`.

    Parameters
    ----------
    V : float or array_like
        Flow speed [m/s], nonnegative.
    a : float or array_like
        Speed of sound [m/s], positive.

    Returns
    -------
    float or ndarray
        Mach number [-].

    Notes
    -----
    Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> float(mach_number(170.0, 340.0))
    0.5
    """
    _nonnegative(V=V)
    _positive(a=a)
    return np.asarray(V) / a


def mach_angle(M):
    r"""Half-angle of the Mach cone, :math:`\mu = \arcsin(1/M)`.

    Parameters
    ----------
    M : float or array_like
        Mach number [-], at least one.

    Returns
    -------
    float or ndarray
        Mach angle [rad]; :math:`\pi/2` at :math:`M = 1`.

    Notes
    -----
    Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    References
    ----------
    E. Mach and P. Salcher, "Photographische Fixirung der durch Projectile in
    der Luft eingeleiteten Vorgänge," *Sitzungsber. Akad. Wiss. Wien* 95 (1887).

    Examples
    --------
    >>> round(float(np.degrees(mach_angle(2.0))), 6)
    30.0
    """
    _interval("M", M, 1, np.inf, closed="left")
    return np.arcsin(1 / np.asarray(M))


def stagnation_temperature(T, M, gamma=1.4):
    r"""Isentropic stagnation temperature :math:`T_0 = T(1 + \tfrac{\gamma-1}{2}M^2)`.

    Parameters
    ----------
    T : float or array_like
        Static temperature [K], positive.
    M : float or array_like
        Mach number [-], nonnegative.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Stagnation (total) temperature [K].

    Notes
    -----
    Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(stagnation_temperature(300.0, 2.0)), 6)
    540.0
    """
    _positive(T=T)
    return np.asarray(T) * isentropic_factor(M, gamma)


def stagnation_pressure(p, M, gamma=1.4):
    r"""Isentropic stagnation pressure :math:`p_0 = p(1 + \tfrac{\gamma-1}{2}M^2)^{\gamma/(\gamma-1)}`.

    Parameters
    ----------
    p : float or array_like
        Static pressure [Pa], positive.
    M : float or array_like
        Mach number [-], nonnegative.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Stagnation (total) pressure [Pa].

    Notes
    -----
    Not conserved across a shock wave. Inputs are raw SI numbers; see
    :doc:`/guides/aerodynamics`.

    References
    ----------
    A. J. C. B. de Saint-Venant and L. Wantzel, "Mémoire et expériences sur
    l'écoulement de l'air," *J. École Polytech.* 16 (1839) 85-122.

    Examples
    --------
    >>> round(float(stagnation_pressure(100.0, 2.0)), 3)
    782.445
    """
    _positive(p=p)
    return np.asarray(p) * isentropic_factor(M, gamma) ** (gamma / (gamma - 1))


def prandtl_glauert_cp(cp0, M):
    r"""Prandtl-Glauert compressibility correction :math:`C_p = C_{p,0}/\sqrt{1 - M_\infty^2}`.

    Parameters
    ----------
    cp0 : float or array_like
        Incompressible pressure coefficient [-].
    M : float or array_like
        Free-stream Mach number [-] in ``[0, 1)``.

    Returns
    -------
    float or ndarray
        Compressible pressure coefficient [-].

    Notes
    -----
    Linearized subsonic theory; it diverges as :math:`M \to 1` and is
    unreliable above roughly :math:`M = 0.7`. Inputs are raw SI numbers; see
    :doc:`/guides/aerodynamics`.

    References
    ----------
    H. Glauert, "The effect of compressibility on the lift of an aerofoil,"
    *Proc. R. Soc. Lond. A* 118 (1928) 113-119.

    Examples
    --------
    >>> round(float(prandtl_glauert_cp(-1.0, 0.6)), 6)
    -1.25
    """
    _finite(cp0=cp0)
    _interval("M", M, 0, 1, closed="left")
    return np.asarray(cp0) / np.sqrt(1 - np.asarray(M) ** 2)


def normal_shock_mach2(M1, gamma=1.4):
    r"""Downstream Mach number of a normal shock.

    .. math::

        M_2^2 = \frac{1 + \tfrac{\gamma-1}{2}M_1^2}{\gamma M_1^2 - \tfrac{\gamma-1}{2}}

    Parameters
    ----------
    M1 : float or array_like
        Upstream Mach number [-], at least one.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Downstream Mach number [-], at most one.

    Notes
    -----
    Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    References
    ----------
    W. J. M. Rankine, *Phil. Trans. R. Soc.* 160 (1870) 277-288;
    H. Hugoniot, *J. École Polytech.* 57 (1887) 3-97; 58 (1889) 1-125.

    Examples
    --------
    >>> round(float(normal_shock_mach2(2.0)), 6)
    0.57735
    """
    _interval("M1", M1, 1, np.inf, closed="left")
    check_gamma(gamma)
    m1sq = np.asarray(M1) ** 2
    return np.sqrt((1 + (gamma - 1) * m1sq / 2) / (gamma * m1sq - (gamma - 1) / 2))


def normal_shock_pressure_ratio(M1, gamma=1.4):
    r"""Static pressure ratio across a normal shock, :math:`p_2/p_1 = 1 + \tfrac{2\gamma}{\gamma+1}(M_1^2 - 1)`.

    Parameters
    ----------
    M1 : float or array_like
        Upstream Mach number [-], at least one.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Pressure ratio [-], at least one.

    Notes
    -----
    Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(normal_shock_pressure_ratio(2.0)), 6)
    4.5
    """
    _interval("M1", M1, 1, np.inf, closed="left")
    check_gamma(gamma)
    return 1 + 2 * gamma / (gamma + 1) * (np.asarray(M1) ** 2 - 1)


def area_mach_ratio(M, gamma=1.4):
    r"""Isentropic area ratio :math:`A/A^*` of a quasi-one-dimensional nozzle.

    .. math::

        \frac{A}{A^*} = \frac{1}{M}\left[\frac{2}{\gamma+1}\left(1 + \frac{\gamma-1}{2}M^2\right)\right]^{\frac{\gamma+1}{2(\gamma-1)}}

    Parameters
    ----------
    M : float or array_like
        Mach number [-], positive.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Area relative to the sonic throat [-]; one at :math:`M = 1`. Each
        ratio above one corresponds to one subsonic and one supersonic Mach
        number; this function does not invert.

    Notes
    -----
    Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(area_mach_ratio(2.0)), 4)
    1.6875
    """
    _positive(M=M)
    M = np.asarray(M)
    return (1 / M) * (2 / (gamma + 1) * isentropic_factor(M, gamma)) ** ((gamma + 1) / (2 * (gamma - 1)))
