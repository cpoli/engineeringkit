r"""Convection: Newton's law of cooling and the classic Nusselt-number correlations.

Newton (1701) made the heat lost by a hot body proportional to its excess
temperature over the surroundings; the coefficient :math:`h` absorbs all
the fluid mechanics. Dittus and Boelter (1930) correlated it for turbulent
flow in tubes, and Churchill and Chu (1975) for natural convection on a
vertical plate over every Rayleigh number.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.heat_transfer.utils.constants import STANDARD_GRAVITY

__all__ = [
    "convection_heat_rate",
    "convection_resistance",
    "dittus_boelter_nusselt",
    "rayleigh_number",
    "churchill_chu_nusselt",
]


def convection_heat_rate(h, A, Ts, Tinf):
    r"""Newton's law of cooling, :math:`\dot Q = hA(T_s - T_\infty)`.

    Parameters
    ----------
    h : float or array_like
        Convection coefficient [W/(m^2 K)], nonnegative.
    A : float or array_like
        Surface area [m^2], positive.
    Ts : float or array_like
        Surface temperature [K].
    Tinf : float or array_like
        Fluid temperature far from the surface [K].

    Returns
    -------
    float or ndarray
        Heat rate from the surface to the fluid [W]; negative if the fluid is hotter.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    References
    ----------
    [I. Newton], "Scala graduum caloris," *Phil. Trans. R. Soc.* 22 (1701) 824-829.

    Examples
    --------
    >>> convection_heat_rate(10.0, 2.0, 400.0, 300.0)
    2000.0
    """
    _finite(h=h, A=A, Ts=Ts, Tinf=Tinf)
    _nonnegative(h=h, Ts=Ts, Tinf=Tinf)
    _positive(A=A)
    return h * A * (Ts - Tinf)


def convection_resistance(h, A=1.0):
    r"""Convection (film) resistance, :math:`R = 1/(hA)`.

    Parameters
    ----------
    h : float or array_like
        Convection coefficient [W/(m^2 K)], positive.
    A : float or array_like, optional
        Surface area [m^2], positive; the default 1 gives the resistance per unit area.

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
    >>> convection_resistance(10.0, 2.0)
    0.05
    """
    _finite(h=h, A=A)
    _positive(h=h, A=A)
    return 1 / (h * A)


def dittus_boelter_nusselt(Re, Pr, heating=True):
    r"""Dittus-Boelter correlation for fully developed turbulent tube flow, :math:`Nu = 0.023\,Re^{0.8} Pr^{n}`.

    Parameters
    ----------
    Re : float or array_like
        Reynolds number on the tube diameter [-], positive; the correlation
        is meant for :math:`Re \gtrsim 10^4` (not checked).
    Pr : float or array_like
        Prandtl number [-], positive; meant for :math:`0.6 \lesssim Pr \lesssim 160`.
    heating : bool, optional
        ``True`` (the fluid is heated, :math:`n = 0.4`) or ``False``
        (cooled, :math:`n = 0.3`).

    Returns
    -------
    float or ndarray
        Nusselt number :math:`hD/k` [-]; typical error is about 25%.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. The coefficient 0.023 is McAdams's
    rounding; Dittus and Boelter's paper gave 0.0243 (heating) and 0.0265
    (cooling).

    References
    ----------
    F. W. Dittus and L. M. K. Boelter, "Heat transfer in automobile radiators
    of the tubular type," *Univ. Calif. Publ. Eng.* 2 (1930) 443-461;
    R. H. S. Winterton, "Where did the Dittus and Boelter equation come
    from?", *Int. J. Heat Mass Transfer* 41 (1998) 809-810.

    Examples
    --------
    >>> round(float(dittus_boelter_nusselt(1e5, 0.7)), 1)
    199.4
    """
    _positive(Re=Re, Pr=Pr)
    n = 0.4 if heating else 0.3
    return 0.023 * np.asarray(Re, dtype=float) ** 0.8 * np.asarray(Pr, dtype=float) ** n


def rayleigh_number(beta, dT, L, nu, alpha, g=STANDARD_GRAVITY):
    r"""Rayleigh number :math:`Ra = g\beta\,\Delta T L^3/(\nu\alpha)` for natural convection.

    Parameters
    ----------
    beta : float or array_like
        Volumetric expansion coefficient [1/K], positive (:math:`1/T` for an ideal gas).
    dT : float or array_like
        Magnitude of the surface-to-fluid temperature difference [K], nonnegative.
    L : float or array_like
        Characteristic length (plate height) [m], positive.
    nu : float or array_like
        Kinematic viscosity [m^2/s], positive.
    alpha : float or array_like
        Thermal diffusivity [m^2/s], positive.
    g : float, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Rayleigh number [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> f"{float(rayleigh_number(1 / 300, 20.0, 0.5, 1.6e-5, 2.2e-5)):.3e}"
    '2.322e+08'
    """
    _positive(beta=beta, L=L, nu=nu, alpha=alpha, g=g)
    _nonnegative(dT=dT)
    return g * beta * np.asarray(dT, dtype=float) * np.asarray(L, dtype=float) ** 3 / (nu * alpha)


def churchill_chu_nusselt(Ra, Pr, laminar=False):
    r"""Churchill-Chu average Nusselt number for natural convection on a vertical plate.

    .. math::

        Nu = \left\{0.825 + \frac{0.387\,Ra^{1/6}}{[1 + (0.492/Pr)^{9/16}]^{8/27}}\right\}^2

    holds for every Rayleigh number. With ``laminar=True`` the more accurate
    laminar form :math:`Nu = 0.68 + 0.670\,Ra^{1/4}/[1 + (0.492/Pr)^{9/16}]^{4/9}`
    (for :math:`Ra \lesssim 10^9`) is returned instead.

    Parameters
    ----------
    Ra : float or array_like
        Rayleigh number on the plate height [-], nonnegative (see :func:`rayleigh_number`).
    Pr : float or array_like
        Prandtl number [-], positive.
    laminar : bool, optional
        Use the laminar-only correlation.

    Returns
    -------
    float or ndarray
        Average Nusselt number :math:`\bar h L/k` [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Fluid properties are taken at the film
    temperature :math:`(T_s + T_\infty)/2`.

    References
    ----------
    S. W. Churchill and H. H. S. Chu, "Correlating equations for laminar and
    turbulent free convection from a vertical plate," *Int. J. Heat Mass
    Transfer* 18 (1975) 1323-1329.

    Examples
    --------
    >>> round(float(churchill_chu_nusselt(1e9, 0.71)), 1)
    122.9
    """
    _nonnegative(Ra=Ra)
    _positive(Pr=Pr)
    Ra = np.asarray(Ra, dtype=float)
    psi = 1 + (0.492 / np.asarray(Pr, dtype=float)) ** (9 / 16)
    if laminar:
        return 0.68 + 0.670 * Ra**0.25 / psi ** (4 / 9)
    return (0.825 + 0.387 * Ra ** (1 / 6) / psi ** (8 / 27)) ** 2
