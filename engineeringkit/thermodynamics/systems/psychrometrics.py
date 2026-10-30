r"""Psychrometrics: moist air as an ideal mixture of dry air and water vapour.

Carrier's psychrometric chart plots the humidity ratio against dry-bulb
temperature, with curves of constant relative humidity. A parcel of air
cooled at constant pressure keeps its humidity ratio until it reaches its
dew point, where the vapour starts to condense. Air-conditioning processes
(mixing, cooling and dehumidifying coils, evaporative cooling towers) are
steady-flow mass balances for dry air and water and an energy balance in
moist-air enthalpy per kilogram of dry air.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.thermodynamics.systems.steam import saturation_pressure, saturation_temperature, steam_properties

__all__ = [
    "MOLAR_MASS_RATIO",
    "humidity_ratio",
    "relative_humidity",
    "dew_point",
    "moist_air_enthalpy",
    "adiabatic_saturation_temperature",
    "adiabatic_mixing",
    "cooling_coil",
    "cooling_tower",
]

MOLAR_MASS_RATIO = 0.621945
"""Ratio of the molar masses of water and dry air, 18.015268/28.966 [-]."""


_CP_AIR, _CP_VAPOUR, _H_G0 = 1005.0, 1820.0, 2501e3  # J/(kg K), J/(kg K), J/kg at 0 °C


def _liquid_enthalpy(T):
    """Enthalpy of saturated liquid water at ``T`` [J/kg] (IF97), on the same reference as moist air."""
    # IF97's region 4 gives 611.2127 Pa at 273.15 K, a hair under the triple-point pressure x needs
    return steam_properties(np.maximum(saturation_pressure(T), 611.213), x=0.0).h


def _temperature(h, w):
    """Dry-bulb temperature [K] of moist air with enthalpy ``h`` [J/kg dry air] and humidity ratio ``w``."""
    return (h - _H_G0 * w) / (_CP_AIR + _CP_VAPOUR * w) + 273.15


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


def moist_air_enthalpy(T, w):
    r"""Enthalpy of moist air per kilogram of dry air, :math:`h = c_{p,a}t + w(h_{g,0} + c_{p,v}t)`.

    Parameters
    ----------
    T : float or array_like
        Dry-bulb temperature [K], positive.
    w : float or array_like
        Humidity ratio [kg/kg], nonnegative.

    Returns
    -------
    float or ndarray
        Enthalpy [J/kg dry air], zero for dry air at 0 °C, with
        :math:`c_{p,a} = 1005`, :math:`c_{p,v} = 1820` J/(kg K) and
        :math:`h_{g,0} = 2501` kJ/kg, :math:`t` in °C.

    Notes
    -----
    Liquid water at 0 °C has zero enthalpy on this scale, as in IF97, so
    condensate enthalpies come from :func:`steam_properties`. Inputs are
    raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> round(float(moist_air_enthalpy(298.15, 0.01)) / 1e3, 3)
    50.59
    """
    _positive(T=T)
    _nonnegative(w=w)
    t = np.asarray(T, dtype=float) - 273.15
    return (_CP_AIR * t + np.asarray(w, dtype=float) * (_H_G0 + _CP_VAPOUR * t))[()]


def adiabatic_saturation_temperature(T, w, p=101325.0):
    r"""Temperature reached by air saturated adiabatically with liquid water at that temperature.

    Parameters
    ----------
    T : float
        Dry-bulb temperature [K], above 273.15 K.
    w : float
        Humidity ratio [kg/kg], nonnegative and at most saturation.
    p : float, optional
        Total pressure [Pa], positive.

    Returns
    -------
    float
        Adiabatic saturation temperature [K], between the dew point and
        the dry-bulb temperature; for air-water it is practically the
        wet-bulb temperature.

    Notes
    -----
    Solves :math:`h(T, w) + (w_s - w)\,h_f(T_{as}) = h(T_{as}, w_s)` with
    :math:`w_s` the saturated humidity ratio at :math:`T_{as}`, by Brent's
    method. Raises ``ValueError`` when :math:`T_{as}` falls below 0 °C.
    Every parameter is scalar-only: the search runs for one air state.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> round(adiabatic_saturation_temperature(303.15, 0.01) - 273.15, 2)
    19.6
    """
    _finite(T=T, w=w)
    if np.ndim(T) or np.ndim(w):
        raise ValueError("T and w must be scalars")
    _interval("T", T, 273.15, 647.0, closed="right")
    _nonnegative(w=w)
    _positive(p=p)
    if relative_humidity(T, w, p) > 1:
        raise ValueError("the air is supersaturated")
    h = moist_air_enthalpy(T, w)

    def residual(T_as):
        w_s = humidity_ratio(T_as, 1.0, p)
        return h + (w_s - w) * _liquid_enthalpy(T_as) - moist_air_enthalpy(T_as, w_s)

    if residual(273.15) < 0:
        raise ValueError("the adiabatic saturation temperature is below 0 °C")
    return float(brentq(residual, 273.15, T, xtol=1e-9)) if residual(T) < 0 else float(T)


def adiabatic_mixing(m_a1, T1, w1, m_a2, T2, w2):
    r"""Adiabatic mixing of two moist-air streams.

    Parameters
    ----------
    m_a1, m_a2 : float or array_like
        Dry-air mass flow of each stream [kg/s], nonnegative, not both zero.
    T1, T2 : float or array_like
        Dry-bulb temperature of each stream [K], positive.
    w1, w2 : float or array_like
        Humidity ratio of each stream [kg/kg], nonnegative.

    Returns
    -------
    tuple of float or ndarray
        Outlet temperature [K] and humidity ratio [kg/kg]. A relative
        humidity above one at the outlet means fog forms.

    Notes
    -----
    Dry air, water and enthalpy are conserved, so the outlet lies on the
    straight line between the inlets on the psychrometric chart:
    :math:`w_3` and :math:`h_3` are the dry-air-weighted means. Inputs are
    raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> T3, w3 = adiabatic_mixing(1.0, 305.15, 0.015, 3.0, 287.15, 0.006)
    >>> round(T3 - 273.15, 2), round(w3, 5)
    (18.55, 0.00825)
    """
    _nonnegative(m_a1=m_a1, m_a2=m_a2, w1=w1, w2=w2)
    _positive(T1=T1, T2=T2, **{"m_a1 + m_a2": m_a1 + m_a2})
    m_a1, m_a2, w1, w2 = (np.asarray(v, dtype=float) for v in (m_a1, m_a2, w1, w2))
    total = m_a1 + m_a2
    w3 = (m_a1 * w1 + m_a2 * w2) / total
    h3 = (m_a1 * moist_air_enthalpy(T1, w1) + m_a2 * moist_air_enthalpy(T2, w2)) / total
    return _temperature(h3, w3)[()], w3[()]


def cooling_coil(T_in, phi_in, T_out, p=101325.0):
    r"""Cooling (and dehumidifying) coil: heat and condensate per kilogram of dry air.

    Parameters
    ----------
    T_in : float or array_like
        Inlet dry-bulb temperature [K], 273.15 K to 647 K.
    phi_in : float or array_like
        Inlet relative humidity [-] in [0, 1].
    T_out : float or array_like
        Outlet temperature [K], at least 273.15 K and at most ``T_in``.
    p : float, optional
        Total pressure [Pa], positive.

    Returns
    -------
    tuple of float or ndarray
        Outlet humidity ratio [kg/kg], heat removed [J/kg dry air] and
        condensate [kg/kg dry air].

    Notes
    -----
    Above the dew point the air is cooled at constant humidity ratio;
    below it, it leaves saturated at ``T_out`` and the excess water drains
    as liquid at ``T_out``:
    :math:`q = h_1 - h_2 - (w_1 - w_2)h_f(T_2)`. Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> w2, q, water = cooling_coil(303.15, 0.8, 287.15)
    >>> round(w2, 5), round(q / 1e3, 1), round(water, 5)
    (0.00997, 45.3, 0.0116)
    """
    _finite(T_in=T_in, phi_in=phi_in, T_out=T_out)
    _interval("T_out", T_out, 273.15, np.inf)
    _nonnegative(**{"T_in - T_out": np.asarray(T_in, dtype=float) - T_out})
    w1 = humidity_ratio(T_in, phi_in, p)
    w2 = np.minimum(w1, humidity_ratio(T_out, 1.0, p))
    water = w1 - w2
    q = moist_air_enthalpy(T_in, w1) - moist_air_enthalpy(T_out, w2) - water * _liquid_enthalpy(T_out)
    return w2[()], np.asarray(q)[()], water[()]


def cooling_tower(m_water, T_water_in, T_water_out, T_air_in, phi_air_in, T_air_out, phi_air_out=1.0, p=101325.0):
    r"""Wet cooling tower: dry-air flow needed and make-up water evaporated.

    Parameters
    ----------
    m_water : float or array_like
        Warm water flow entering the tower [kg/s], positive.
    T_water_in, T_water_out : float or array_like
        Water temperature entering (warm) and leaving (cooled) [K], with
        ``T_water_out < T_water_in``, both at least 273.15 K.
    T_air_in, phi_air_in : float or array_like
        Inlet air temperature [K] and relative humidity [-].
    T_air_out, phi_air_out : float or array_like, optional
        Outlet air temperature [K] and relative humidity [-]; the air
        usually leaves saturated.
    p : float, optional
        Total pressure [Pa], positive.

    Returns
    -------
    tuple of float or ndarray
        Dry-air mass flow [kg/s] and make-up water flow, the water
        evaporated into the air [kg/s].

    Notes
    -----
    From the water balance :math:`\dot m_{evap} = \dot m_a(w_2 - w_1)` and
    the energy balance
    :math:`\dot m_a(h_2 - h_1) = \dot m_w h_{f,in} - (\dot m_w - \dot m_{evap})h_{f,out}`,
    with liquid enthalpies from IF97. The air must gain enthalpy and
    moisture. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    Y. A. Çengel and M. A. Boles, *Thermodynamics: An Engineering
    Approach*, 9th ed., New York: McGraw-Hill (2019), Ch. 14.

    Examples
    --------
    >>> m_air, makeup = cooling_tower(100.0, 308.15, 295.15, 293.15, 0.6, 303.15)
    >>> round(m_air, 1), round(makeup, 2)
    (97.6, 1.8)
    """
    _positive(m_water=m_water)
    _finite(T_water_in=T_water_in, T_water_out=T_water_out)
    _interval("T_water_out", T_water_out, 273.15, np.inf)
    _positive(**{"T_water_in - T_water_out": np.asarray(T_water_in, dtype=float) - T_water_out})
    w1 = humidity_ratio(T_air_in, phi_air_in, p)
    w2 = humidity_ratio(T_air_out, phi_air_out, p)
    h1, h2 = moist_air_enthalpy(T_air_in, w1), moist_air_enthalpy(T_air_out, w2)
    hf_in, hf_out = _liquid_enthalpy(T_water_in), _liquid_enthalpy(T_water_out)
    gain = (h2 - h1) - (w2 - w1) * hf_out
    if np.any(w2 <= w1) or np.any(gain <= 0):
        raise ValueError("the air must leave with more moisture and enthalpy than it enters")
    m_air = m_water * (hf_in - hf_out) / gain
    return np.asarray(m_air)[()], np.asarray(m_air * (w2 - w1))[()]
