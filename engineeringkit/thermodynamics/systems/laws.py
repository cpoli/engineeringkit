r"""The first and second laws: energy balances, entropy and reversible processes.

Joule's experiments made heat and work two forms of one conserved quantity,
:math:`Q = \Delta U + W`; Clausius's inequality :math:`\oint \delta Q/T \le 0`
introduced entropy and fixed which processes can run backwards.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.thermodynamics.systems.gases import R_UNIVERSAL

__all__ = [
    "first_law_heat",
    "specific_internal_energy_ideal_gas",
    "specific_enthalpy_ideal_gas",
    "entropy_change_isothermal_ideal_gas",
    "reversible_isothermal_work_ideal_gas",
    "isentropic_pressure_volume",
    "isentropic_temperature_ratio",
    "clausius_inequality",
]


def first_law_heat(delta_u, work_by_system):
    r"""Heat supplied to a closed system, :math:`Q = \Delta U + W`.

    Parameters
    ----------
    delta_u : float or array_like
        Change of internal energy [J].
    work_by_system : float or array_like
        Work done *by* the system [J]; negative when work is done on it.

    Returns
    -------
    float or ndarray
        Heat added to the system [J]; negative when heat leaves it.

    Notes
    -----
    Kinetic and potential energy changes need their own terms. Inputs are raw
    SI numbers with no dimensional metadata; see :doc:`/guides/thermodynamics`.

    References
    ----------
    J. P. Joule, "On the mechanical equivalent of heat," *Phil. Trans. R.
    Soc.* 140 (1850) 61-82.
    """
    return delta_u + work_by_system


def specific_internal_energy_ideal_gas(T, cv):
    r"""Specific internal energy :math:`u = c_v T` of a calorically perfect gas.

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K], nonnegative.
    cv : float or array_like
        Isochoric specific heat [J/(kg K)], positive.

    Returns
    -------
    float or ndarray
        Specific internal energy [J/kg], referred to zero at 0 K.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.
    """
    _finite(T=T, cv=cv)
    _nonnegative(T=T)
    _positive(cv=cv)
    return cv * T


def specific_enthalpy_ideal_gas(T, cp):
    r"""Specific enthalpy :math:`h = c_p T` of a calorically perfect gas.

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K], nonnegative.
    cp : float or array_like
        Isobaric specific heat [J/(kg K)], positive.

    Returns
    -------
    float or ndarray
        Specific enthalpy [J/kg], referred to zero at 0 K. With the same
        reference, :math:`h - u = (c_p - c_v)T = RT`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.
    """
    _finite(T=T, cp=cp)
    _nonnegative(T=T)
    _positive(cp=cp)
    return cp * T


def entropy_change_isothermal_ideal_gas(n, T, V1, V2):
    r"""Entropy change of an ideal gas at constant temperature, :math:`nR\ln(V_2/V_1)`.

    Parameters
    ----------
    n : float or array_like
        Amount of gas [mol], positive.
    T : float or array_like
        Absolute temperature [K], positive (it does not enter the result).
    V1, V2 : float or array_like
        Initial and final volumes [m^3], positive.

    Returns
    -------
    float or ndarray
        Entropy change [J/K].

    Notes
    -----
    Entropy is a state function, so the result holds for any process between
    the two states, reversible or not. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/thermodynamics`.
    """
    _finite(n=n, T=T, V1=V1, V2=V2)
    _positive(n=n, T=T, V1=V1, V2=V2)
    return n * R_UNIVERSAL * np.log(V2 / V1)


def reversible_isothermal_work_ideal_gas(n, T, V1, V2):
    r"""Work done by an ideal gas in a reversible isothermal process, :math:`nRT\ln(V_2/V_1)`.

    Parameters
    ----------
    n : float or array_like
        Amount of gas [mol], positive.
    T : float or array_like
        Absolute temperature [K], positive.
    V1, V2 : float or array_like
        Initial and final volumes [m^3], positive.

    Returns
    -------
    float or ndarray
        Work done by the gas [J], equal to the heat absorbed, :math:`T\Delta S`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.
    """
    _finite(n=n, T=T, V1=V1, V2=V2)
    _positive(n=n, T=T, V1=V1, V2=V2)
    return n * R_UNIVERSAL * T * np.log(V2 / V1)


def isentropic_pressure_volume(p1, v1, v2, gamma=1.4):
    r"""Pressure after a reversible adiabatic change of volume, :math:`p_1(v_1/v_2)^\gamma`.

    Parameters
    ----------
    p1 : float or array_like
        Initial pressure [Pa], positive.
    v1, v2 : float or array_like
        Initial and final volumes in consistent units (e.g. m^3/kg), positive.
    gamma : float or array_like, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Final pressure [Pa].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.
    """
    _finite(p1=p1, v1=v1, v2=v2, gamma=gamma)
    _positive(p1=p1, v1=v1, v2=v2)
    if np.any(np.asarray(gamma) <= 1):
        raise ValueError("gamma must exceed 1")
    p1, v1, v2, gamma = (np.asarray(x, dtype=float) for x in (p1, v1, v2, gamma))
    return (p1 * (v1 / v2) ** gamma)[()]


def isentropic_temperature_ratio(pressure_ratio, gamma=1.4):
    r"""Temperature ratio across an isentropic pressure change, :math:`T_2/T_1 = (p_2/p_1)^{(\gamma-1)/\gamma}`.

    Parameters
    ----------
    pressure_ratio : float or array_like
        :math:`p_2/p_1` [-], positive.
    gamma : float or array_like, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        :math:`T_2/T_1` [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.
    """
    _finite(pressure_ratio=pressure_ratio, gamma=gamma)
    _positive(pressure_ratio=pressure_ratio)
    if np.any(np.asarray(gamma) <= 1):
        raise ValueError("gamma must exceed 1")
    gamma = np.asarray(gamma, dtype=float)
    return (np.asarray(pressure_ratio, dtype=float) ** ((gamma - 1) / gamma))[()]


def clausius_inequality(cyclic_heat_over_temperature):
    r"""Whether a discretized cycle satisfies Clausius's inequality :math:`\sum \delta Q/T \le 0`.

    Parameters
    ----------
    cyclic_heat_over_temperature : array_like
        The terms :math:`\delta Q_i/T_i` [J/K] around one cycle, heat into the
        system positive and :math:`T_i` the temperature at the boundary where
        it crosses.

    Returns
    -------
    bool
        ``True`` if the sum is at most :math:`10^{-12}` J/K (zero within
        rounding for a reversible cycle, negative for an irreversible one).

    Notes
    -----
    The function checks a supplied sum and does not construct a path.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    R. Clausius, "Über eine veränderte Form des zweiten Hauptsatzes der
    mechanischen Wärmetheorie," *Ann. Phys.* 93 (1854) 481-506.
    """
    return float(np.sum(cyclic_heat_over_temperature)) <= 1e-12
