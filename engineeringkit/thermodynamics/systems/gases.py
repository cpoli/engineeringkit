r"""Equations of state: the ideal gas, van der Waals's real gas and the Joule-Thomson effect.

Clapeyron combined Boyle's and Gay-Lussac's laws into :math:`pV = nRT`;
van der Waals corrected it for the volume of the molecules and their mutual
attraction, which also explains the critical point and why a throttled gas
cools below its inversion temperature.
"""

from __future__ import annotations

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "R_UNIVERSAL",
    "ideal_gas_pressure",
    "ideal_gas_temperature",
    "van_der_waals_pressure",
    "van_der_waals_critical_point",
    "joule_thomson_coefficient_vdw",
]

R_UNIVERSAL = 8.31446261815324
"""Molar gas constant [J/(mol K)] (exact since the 2019 SI)."""


def ideal_gas_pressure(n, T, V):
    r"""Ideal-gas pressure :math:`p = nRT/V`.

    Parameters
    ----------
    n : float or array_like
        Amount of gas [mol], positive.
    T : float or array_like
        Absolute temperature [K], positive.
    V : float or array_like
        Volume [m^3], positive.

    Returns
    -------
    float or ndarray
        Pressure [Pa].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    É. Clapeyron, "Mémoire sur la puissance motrice de la chaleur," *J. École
    Polytechnique* 14 (1834) 153-190.

    Examples
    --------
    >>> round(float(ideal_gas_pressure(1.0, 273.15, 0.022414)), 0)
    101325.0
    """
    _finite(n=n, T=T, V=V)
    _positive(n=n, T=T, V=V)
    return n * R_UNIVERSAL * T / V


def ideal_gas_temperature(P, V, n):
    r"""Ideal-gas temperature :math:`T = pV/(nR)`.

    Parameters
    ----------
    P : float or array_like
        Pressure [Pa], positive.
    V : float or array_like
        Volume [m^3], positive.
    n : float or array_like
        Amount of gas [mol], positive.

    Returns
    -------
    float or ndarray
        Absolute temperature [K].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.
    """
    _finite(P=P, V=V, n=n)
    _positive(P=P, V=V, n=n)
    return P * V / (n * R_UNIVERSAL)


def van_der_waals_pressure(n, T, V, a, b):
    r"""Van der Waals pressure :math:`p = nRT/(V - nb) - an^2/V^2`.

    Parameters
    ----------
    n : float or array_like
        Amount of gas [mol], positive.
    T : float or array_like
        Absolute temperature [K], positive.
    V : float or array_like
        Volume [m^3], larger than the excluded volume :math:`nb`.
    a : float
        Attraction constant [Pa m^6/mol^2], nonnegative.
    b : float
        Excluded (co-)volume [m^3/mol], nonnegative.

    Returns
    -------
    float or ndarray
        Pressure [Pa]. With ``a = b = 0`` this is the ideal-gas pressure.

    Notes
    -----
    Below the critical temperature an isotherm has a loop in which
    :math:`\partial p/\partial V > 0`; those states are unstable, and the
    real fluid splits into liquid and vapour at the pressure given by
    Maxwell's equal-area rule. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/thermodynamics`.

    References
    ----------
    J. D. van der Waals, *Over de continuïteit van den gas- en
    vloeistoftoestand*, thesis, Leiden (1873).

    Examples
    --------
    >>> round(float(van_der_waals_pressure(1.0, 300.0, 1e-3, 0.3640, 4.267e-5)), -2)
    2241500.0
    """
    _positive(n=n, T=T, V=V)
    _nonnegative(a=a, b=b)
    excluded = V - n * b
    _positive(**{"V - n*b": excluded})
    return n * R_UNIVERSAL * T / excluded - a * n**2 / V**2


def van_der_waals_critical_point(a, b):
    r"""Critical temperature, pressure and molar volume of a van der Waals gas.

    Parameters
    ----------
    a : float
        Attraction constant [Pa m^6/mol^2], positive.
    b : float
        Excluded volume [m^3/mol], positive.

    Returns
    -------
    tuple of float
        :math:`T_c = 8a/(27Rb)` [K], :math:`p_c = a/(27b^2)` [Pa] and
        :math:`v_c = 3b` [m^3/mol], the point where the critical isotherm has
        a horizontal inflection.

    Notes
    -----
    The critical compressibility :math:`p_c v_c/(RT_c) = 3/8` is the same for
    every van der Waals gas; real gases are nearer 0.27-0.29. Inputs are raw
    SI numbers with no dimensional metadata; see :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> Tc, pc, vc = van_der_waals_critical_point(0.3640, 4.267e-5)
    >>> round(Tc, 1), round(pc / 1e6, 2)
    (304.0, 7.4)
    """
    _positive(a=a, b=b)
    return 8 * a / (27 * R_UNIVERSAL * b), a / (27 * b**2), 3 * b


def joule_thomson_coefficient_vdw(T, a, b, cp_molar):
    r"""Low-pressure Joule-Thomson coefficient of a van der Waals gas, :math:`(2a/RT - b)/c_{p,m}`.

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K], positive.
    a : float
        Attraction constant [Pa m^6/mol^2], nonnegative.
    b : float
        Excluded volume [m^3/mol], nonnegative.
    cp_molar : float
        Molar isobaric heat capacity [J/(mol K)], positive.

    Returns
    -------
    float or ndarray
        :math:`\mu_{JT} = (\partial T/\partial p)_h` [K/Pa]. Positive values
        mean the gas cools when throttled; the sign changes at the inversion
        temperature :math:`2a/(Rb)`.

    Notes
    -----
    This is the zero-pressure limit of
    :math:`\mu_{JT} = [T(\partial v/\partial T)_p - v]/c_p` for the van der
    Waals equation; at higher pressures the inversion temperature falls.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    J. P. Joule and W. Thomson, "On the thermal effects experienced by air in
    rushing through small apertures," *Phil. Mag.* 4 (1852) 481-492.

    Examples
    --------
    >>> round(float(joule_thomson_coefficient_vdw(300.0, 0.137, 3.87e-5, 29.1)) * 1e5, 3)
    0.244
    """
    _positive(T=T, cp_molar=cp_molar)
    _nonnegative(a=a, b=b)
    return (2 * a / (R_UNIVERSAL * T) - b) / cp_molar
