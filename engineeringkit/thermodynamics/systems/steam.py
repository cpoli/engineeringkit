r"""Properties of water and steam from the IAPWS-IF97 industrial formulation.

All five regions are implemented: 1 (compressed liquid), 2 (vapour), 3
(near the critical point), 4 (saturation line) and 5 (high-temperature
steam). That covers 273.15 K to 1073.15 K up to 100 MPa, and on to
2273.15 K up to 50 MPa: subcritical and supercritical Rankine cycles,
condensers, boilers, psychrometrics above 0 °C and steam in gas turbines.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.optimize import brentq

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive
from engineeringkit.thermodynamics.core.base import SteamState
from engineeringkit.thermodynamics.utils import if97

__all__ = ["saturation_pressure", "saturation_temperature", "steam_properties"]

_T_MIN, _T_25, _T_MAX, _T_13 = 273.15, 1073.15, 2273.15, 623.15
_P_MAX, _P_MAX_5 = 100e6, 50e6
_P_TRIPLE = 611.213
_P_SAT_13 = float(if97.saturation_pressure_mpa(_T_13)) * 1e6


def saturation_pressure(T):
    r"""Saturation (vapour) pressure of water, IF97 region 4.

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K], from 273.15 K to the critical 647.096 K.

    Returns
    -------
    float or ndarray
        Saturation pressure [Pa].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    W. Wagner et al., "The IAPWS Industrial Formulation 1997 for the
    thermodynamic properties of water and steam," *J. Eng. Gas Turbines
    Power* 122 (2000) 150-182.

    Examples
    --------
    >>> round(float(saturation_pressure(373.15)), 1)
    101418.0
    """
    _interval("T", T, _T_MIN, if97.T_CRIT)
    return if97.saturation_pressure_mpa(T) * 1e6


def saturation_temperature(p):
    r"""Saturation temperature of water, IF97 region 4.

    Parameters
    ----------
    p : float or array_like
        Pressure [Pa], from 611.213 Pa (0 °C) to the critical 22.064 MPa.

    Returns
    -------
    float or ndarray
        Saturation (boiling) temperature [K].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> round(float(saturation_temperature(101325.0)) - 273.15, 3)
    99.974
    """
    _interval("p", p, _P_TRIPLE, if97.P_CRIT)
    return if97.saturation_temperature_k(np.asarray(p) / 1e6)


def _state(region, p, T, props, x=math.nan):
    return SteamState(
        p=float(p),
        T=float(T),
        v=float(props["v"]),
        h=float(props["h"]),
        s=float(props["s"]),
        u=float(props["u"]),
        cp=float(props["cp"]) if region != 4 else math.nan,
        x=float(x),
        region=region,
    )


def _region3(phase):
    """Region 3 properties at ``(T, p_mpa)`` on the ``phase`` branch of the isotherm."""
    return lambda T, p_mpa: if97.region3(if97.region3_density(T, p_mpa, phase), T)


def _single_phase(p, T):
    p_mpa = p / 1e6
    if T > _T_25:
        return _state(5, p, T, if97.region5(T, p_mpa))
    vapour = T < if97.T_CRIT and p < float(if97.saturation_pressure_mpa(T)) * 1e6
    if T <= _T_13:
        if vapour:
            return _state(2, p, T, if97.region2(T, p_mpa))
        return _state(1, p, T, if97.region1(T, p_mpa))
    if T <= 863.15 and p > float(if97.b23_pressure_mpa(T)) * 1e6:
        return _state(3, p, T, _region3("vapour" if vapour else "liquid")(T, p_mpa))
    return _state(2, p, T, if97.region2(T, p_mpa))


def _saturated(p):
    """Saturation temperature [K] and the liquid and vapour properties at ``p`` [Pa]."""
    T = float(if97.saturation_temperature_k(p / 1e6))
    if p <= _P_SAT_13:
        return T, if97.region1(T, p / 1e6), if97.region2(T, p / 1e6)
    return T, _region3("liquid")(T, p / 1e6), _region3("vapour")(T, p / 1e6)


def _wet(p, x):
    T, liquid, vapour = _saturated(p)
    mix = {k: (1 - x) * liquid[k] + x * vapour[k] for k in ("v", "h", "s", "u")}
    mix["cp"] = math.nan
    return _state(4, p, T, mix, x)


def _from_property(p, name, value):
    """Invert ``name`` (``"h"`` or ``"s"``), which rises with temperature at fixed pressure."""
    # (region, properties(T, p_mpa), lowest T, highest T) along the isobar, in order of temperature
    if p < _P_TRIPLE:
        segments = [(2, if97.region2, _T_MIN, _T_25)]
    elif p <= _P_SAT_13:
        T_sat = float(if97.saturation_temperature_k(p / 1e6))
        segments = [(1, if97.region1, _T_MIN, T_sat), (2, if97.region2, T_sat, _T_25)]
    else:
        T_23 = float(if97.b23_temperature_k(p / 1e6))
        if p <= if97.P_CRIT:
            T_sat = float(if97.saturation_temperature_k(p / 1e6))
            near_critical = [(3, _region3("liquid"), _T_13, T_sat), (3, _region3("vapour"), T_sat, T_23)]
        else:
            near_critical = [(3, _region3("liquid"), _T_13, T_23)]
        segments = [(1, if97.region1, _T_MIN, _T_13), *near_critical, (2, if97.region2, T_23, _T_25)]
    if p <= _P_MAX_5:
        segments.append((5, if97.region5, _T_25, _T_MAX))
    if _P_TRIPLE <= p <= if97.P_CRIT:
        _, liquid, vapour = _saturated(p)
        f, g = float(liquid[name]), float(vapour[name])
        if f <= value <= g:
            return _wet(p, (value - f) / (g - f) if g > f else 0.0)
    for i, (region, props, low, high) in enumerate(segments):
        lo, hi = (float(props(T, p / 1e6)[name]) - value for T in (low, high))
        if hi >= 0:
            break
    if hi < 0 or (i == 0 and lo > 0):
        raise ValueError(f"{name} is outside the range covered by IF97 at this pressure")
    # a value between two segments (the regions meet to within the formulation's tolerance) takes the boundary
    T = low if lo >= 0 else brentq(lambda t: float(props(t, p / 1e6)[name]) - value, low, high, xtol=1e-10)
    return _state(region, p, T, props(T, p / 1e6))


def _scalar_state(p, name, value):
    if name == "T":
        if value > _T_25 and p > _P_MAX_5:
            raise ValueError("above 1073.15 K (IF97 region 5) p must be at most 50 MPa")
        return _single_phase(p, value)
    if name == "x":
        _interval("p", p, _P_TRIPLE, if97.P_CRIT)
        return _wet(p, value)
    return _from_property(p, name, value)


def steam_properties(p, T=None, *, x=None, s=None, h=None):
    r"""State of water or steam from pressure and one other property (IAPWS-IF97).

    Parameters
    ----------
    p : float or array_like
        Pressure [Pa], positive and at most 100 MPa (50 MPa above
        1073.15 K).
    T : float or array_like, optional
        Temperature [K], 273.15 K to 2273.15 K, for a single-phase state.
    x : float or array_like, optional
        Vapour quality [-] in [0, 1] for a saturated mixture; needs
        611.213 Pa <= ``p`` <= 22.064 MPa (the critical pressure).
    s : float or array_like, optional
        Specific entropy [J/(kg K)], e.g. after an isentropic turbine.
    h : float or array_like, optional
        Specific enthalpy [J/kg], e.g. after a throttling valve.

    Exactly one of ``T``, ``x``, ``s`` and ``h`` is given. Arrays broadcast
    against ``p``.

    Returns
    -------
    SteamState
        Pressure, temperature, specific volume, enthalpy, entropy, internal
        energy, isobaric heat capacity (NaN in the dome), quality (NaN outside
        it) and the IF97 region. Scalar inputs give float fields; array
        inputs give arrays of the broadcast shape, one state per element.

    Notes
    -----
    Given ``T``, the phase is liquid if ``p`` is at or above the saturation
    pressure and vapour below it; a state exactly on the saturation line
    needs ``x``. Region 3 (above 350 °C and above the B23 boundary
    pressure) is a Helmholtz function of density and temperature, so its
    density is found by Brent's method on the isotherm, taking the liquid or
    vapour root of its van der Waals loop below the critical temperature.
    Region 5 (above 1073.15 K) is a Gibbs function like region 2, valid to
    50 MPa. Given ``s`` or ``h``, the temperature is found by Brent's
    method, since both rise monotonically with temperature at fixed
    pressure. Array inputs are evaluated element by element. Inputs are raw
    SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    W. Wagner et al., "The IAPWS Industrial Formulation 1997 for the
    thermodynamic properties of water and steam," *J. Eng. Gas Turbines
    Power* 122 (2000) 150-182; IAPWS, *Revised Release on the IAPWS
    Industrial Formulation 1997*, R7-97 (2007, 2012), for region 5 to 50 MPa.

    Examples
    --------
    >>> st = steam_properties(10e6, 773.15)
    >>> round(st.h / 1e3, 1), round(st.s / 1e3, 4), st.region
    (3375.1, 6.5993, 2)
    >>> steam_properties(1e5, [300.0, 400.0, 1500.0]).region.tolist()
    [1, 2, 5]
    """
    _positive(p=p)
    if np.any(np.asarray(p) > _P_MAX):
        raise ValueError("p must be at most 100 MPa")
    given = {name: value for name, value in (("T", T), ("x", x), ("s", s), ("h", h)) if value is not None}
    if len(given) != 1:
        raise ValueError("give exactly one of T, x, s or h")
    ((name, value),) = given.items()
    _finite(**{name: value})
    if name == "T":
        _interval("T", value, _T_MIN, _T_MAX)
    elif name == "x":
        _interval("x", value, 0, 1)
    if not np.ndim(p) and not np.ndim(value):
        return _scalar_state(float(p), name, float(value))
    pp, vv = np.broadcast_arrays(np.asarray(p, dtype=float), np.asarray(value, dtype=float))
    states = [_scalar_state(float(a), name, float(b)) for a, b in zip(pp.ravel(), vv.ravel())]
    fields: dict[str, Any] = {f: np.array([getattr(st, f) for st in states]).reshape(pp.shape) for f in SteamState.__dataclass_fields__}
    return SteamState(**fields)
