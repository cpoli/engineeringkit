r"""The International Standard Atmosphere (troposphere and lower stratosphere).

Hydrostatic balance :math:`dp/dh = -\rho g` with the ideal-gas law and a
piecewise-linear temperature profile gives closed-form pressure and density:
a power law in the troposphere (constant lapse rate) and an exponential in the
isothermal layer above the tropopause.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import interval as _interval
from engineeringkit.aerodynamics.core.base import AtmosphereState

__all__ = ["isa_atmosphere", "R_AIR", "GAMMA_AIR"]

R_AIR = 287.05287
"""Specific gas constant of dry air [J/(kg K)] used by the ISA."""

GAMMA_AIR = 1.4
"""Ratio of specific heats of air [-]."""

_G0 = 9.80665
_T0, _P0, _LAPSE = 288.15, 101325.0, 0.0065
_H_TROPOPAUSE, _H_TOP = 11000.0, 20000.0


def isa_atmosphere(altitude) -> AtmosphereState:
    r"""Standard-atmosphere temperature, pressure, density and sound speed.

    .. math::

        T = T_0 - L h,\quad p = p_0 (T/T_0)^{g_0/(R L)} \quad (h \le 11\,\mathrm{km}),
        \qquad p = p_{11}\, e^{-g_0 (h - 11\,\mathrm{km})/(R T_{11})} \quad (11 < h \le 20\,\mathrm{km})

    Parameters
    ----------
    altitude : float or array_like
        Geopotential altitude [m] in ``[0, 20000]``.

    Returns
    -------
    AtmosphereState

    Notes
    -----
    Sea-level values are 288.15 K, 101325 Pa and 1.225 kg/m^3 with a lapse
    rate of 6.5 K/km. Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    References
    ----------
    ISO 2533:1975, *Standard Atmosphere*; ICAO Doc 7488 (1952).

    Examples
    --------
    >>> state = isa_atmosphere(11000.0)
    >>> round(float(state.temperature), 2), round(float(state.pressure))
    (216.65, 22632)
    """
    _interval("altitude", altitude, 0, _H_TOP)
    h = np.asarray(altitude, dtype=float)
    exponent = _G0 / (R_AIR * _LAPSE)
    t11 = _T0 - _LAPSE * _H_TROPOPAUSE
    p11 = _P0 * (t11 / _T0) ** exponent
    troposphere = h <= _H_TROPOPAUSE
    T = np.where(troposphere, _T0 - _LAPSE * h, t11)
    p = np.where(troposphere, _P0 * (T / _T0) ** exponent, p11 * np.exp(-_G0 * (h - _H_TROPOPAUSE) / (R_AIR * t11)))
    rho = p / (R_AIR * T)
    a = np.sqrt(GAMMA_AIR * R_AIR * T)
    return AtmosphereState(h, T, p, rho, a)
