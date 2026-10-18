r"""Psychrometrics: moist air as an ideal mixture of dry air and water vapour.

Carrier's psychrometric chart plots the humidity ratio against dry-bulb
temperature, with curves of constant relative humidity. A parcel of air
cooled at constant pressure keeps its humidity ratio until it reaches its
dew point, where the vapour starts to condense.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.thermodynamics.systems.steam import saturation_pressure, saturation_temperature

__all__ = ["MOLAR_MASS_RATIO", "humidity_ratio", "relative_humidity", "dew_point"]

MOLAR_MASS_RATIO = 0.621945
"""Ratio of the molar masses of water and dry air, 18.015268/28.966 [-]."""


def _vapour_pressure(w, p):
    return w * p / (MOLAR_MASS_RATIO + w)


def humidity_ratio(T, phi, p=101325.0):
    r"""Humidity ratio of moist air, :math:`w = 0.622\,\phi p_s/(p - \phi p_s)`.

    Parameters
    ----------
    T : float or array_like
        Dry-bulb temperature [K], 273.15 K to 647 K.
    phi : float or array_like
        Relative humidity [-] in [0, 1].
    p : float, optional
        Total pressure [Pa], positive; sea level by default.

    Returns
    -------
    float or ndarray
        Mass of water vapour per mass of dry air [kg/kg].

    Notes
    -----
    The saturation pressure :math:`p_s(T)` over liquid water comes from
    :func:`saturation_pressure` (IF97). Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/thermodynamics`.

    References
    ----------
    W. H. Carrier, "Rational psychrometric formulae," *Trans. ASME* 33
    (1911) 1005-1053.

    Examples
    --------
    >>> round(float(humidity_ratio(298.15, 0.5)), 5)
    0.00988
    """
    _interval("phi", phi, 0, 1)
    _positive(p=p)
    p_v = np.asarray(phi) * saturation_pressure(T)
    _positive(**{"p - p_v": p - p_v})
    return MOLAR_MASS_RATIO * p_v / (p - p_v)


def relative_humidity(T, w, p=101325.0):
    r"""Relative humidity :math:`\phi = p_v/p_s(T)` from the humidity ratio.

    Parameters
    ----------
    T : float or array_like
        Dry-bulb temperature [K], 273.15 K to 647 K.
    w : float or array_like
        Humidity ratio [kg/kg], nonnegative.
    p : float, optional
        Total pressure [Pa], positive.

    Returns
    -------
    float or ndarray
        Relative humidity [-]; above one the air is supersaturated and
        condensate (fog) forms.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.
    """
    _nonnegative(w=w)
    _positive(p=p)
    return _vapour_pressure(np.asarray(w), p) / saturation_pressure(T)


def dew_point(w, p=101325.0):
    r"""Dew-point temperature: where air of humidity ratio ``w`` saturates on cooling.

    Parameters
    ----------
    w : float or array_like
        Humidity ratio [kg/kg], positive.
    p : float, optional
        Total pressure [Pa], positive.

    Returns
    -------
    float or ndarray
        Dew point [K], the saturation temperature of the vapour's partial
        pressure. Dew points below 0 °C (frost points) are not covered.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> round(float(dew_point(humidity_ratio(298.15, 0.5))) - 273.15, 2)
    13.86
    """
    _positive(w=w, p=p)
    return saturation_temperature(_vapour_pressure(np.asarray(w), p))
