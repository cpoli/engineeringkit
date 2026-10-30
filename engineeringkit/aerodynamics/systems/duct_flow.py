r"""Constant-area duct flow with friction (Fanno) or heat addition (Rayleigh).

Both flows keep the mass flux fixed and drive the Mach number towards one:
friction does it adiabatically, along a *Fanno line*, and heat transfer does
it without friction, along a *Rayleigh line*. Each state is written relative
to the sonic state :math:`(\,)^*` on the same line, so the duct length or
heat that chokes the flow follows from the starred ratios.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import positive as _positive
from engineeringkit.aerodynamics.core.base import DuctFlowState
from engineeringkit.aerodynamics.systems.atmosphere import R_AIR
from engineeringkit.aerodynamics.utils.gas import check_gamma

__all__ = ["fanno_flow", "fanno_choking_length", "rayleigh_flow", "rayleigh_choking_heat"]


def fanno_flow(M, gamma=1.4) -> DuctFlowState:
    r"""Fanno-line ratios of adiabatic flow with wall friction in a constant-area duct.

    .. math::

        \frac{T}{T^*} = \frac{\gamma+1}{2 + (\gamma-1)M^2},\qquad
        \frac{p}{p^*} = \frac1M\sqrt{\frac{T}{T^*}},\qquad
        \frac{f L^*}{D} = \frac{1 - M^2}{\gamma M^2}
        + \frac{\gamma+1}{2\gamma}\ln\frac{(\gamma+1)M^2}{2 + (\gamma-1)M^2}

    Parameters
    ----------
    M : float or array_like
        Mach number [-], positive.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    DuctFlowState
        Ratios to the sonic state, including the friction parameter
        :math:`f L^*/D` with the Darcy friction factor :math:`f`
        (four times the Fanning factor).

    Notes
    -----
    Friction drives subsonic flow up and supersonic flow down towards
    :math:`M = 1`, while :math:`T_0` stays constant and :math:`p_0` falls.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    References
    ----------
    J. D. Anderson, *Modern Compressible Flow*, 3rd ed. (McGraw-Hill,
    2003), sec. 3.9, after G. Fanno's 1904 diploma thesis at the ETH Zürich.

    Examples
    --------
    >>> state = fanno_flow(2.0)
    >>> round(float(state.friction_parameter), 5), round(float(state.pressure_ratio), 5)
    (1.21999, 0.40825)
    """
    _positive(M=M)
    check_gamma(gamma)
    M = np.asarray(M, float)
    m2 = M**2
    t = (gamma + 1) / (2 + (gamma - 1) * m2)
    friction = 4 * ((1 - m2) / (gamma * m2) + (gamma + 1) / (2 * gamma) * np.log((gamma + 1) * m2 / (2 + (gamma - 1) * m2)))
    return DuctFlowState(
        kind="fanno",
        mach=M,
        temperature_ratio=t,
        pressure_ratio=np.sqrt(t) / M,
        density_ratio=1 / (M * np.sqrt(t)),
        stagnation_pressure_ratio=(1 / t) ** ((gamma + 1) / (2 * (gamma - 1))) / M,
        stagnation_temperature_ratio=np.ones_like(M),
        friction_parameter=friction,
    )


def fanno_choking_length(M, D, f, gamma=1.4):
    r"""Duct length :math:`L^*` over which friction brings the flow to :math:`M = 1`.

    Parameters
    ----------
    M : float or array_like
        Inlet Mach number [-], positive.
    D : float or array_like
        Hydraulic diameter [m], positive.
    f : float or array_like
        Darcy friction factor [-], positive (for example from
        :func:`engineeringkit.fluid_engineering.colebrook_friction_factor`).
    gamma : float, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Choking length [m]. A longer duct cannot pass the same mass flow:
        subsonic inlet flow is throttled and supersonic flow shocks.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(fanno_choking_length(0.5, 0.05, 0.02)), 4)
    10.6906
    """
    _positive(D=D, f=f)
    return fanno_flow(M, gamma).friction_parameter * np.asarray(D) / f


def rayleigh_flow(M, gamma=1.4) -> DuctFlowState:
    r"""Rayleigh-line ratios of frictionless flow with heat transfer in a constant-area duct.

    .. math::

        \frac{p}{p^*} = \frac{1+\gamma}{1+\gamma M^2},\qquad
        \frac{T}{T^*} = M^2\left(\frac{p}{p^*}\right)^2,\qquad
        \frac{T_0}{T_0^*} = \frac{(\gamma+1)M^2\,[2 + (\gamma-1)M^2]}{(1+\gamma M^2)^2}

    Parameters
    ----------
    M : float or array_like
        Mach number [-], positive.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    DuctFlowState
        Ratios to the sonic state; ``friction_parameter`` is zero.

    Notes
    -----
    Heating drives the Mach number towards one from either side and
    cooling drives it away. The stagnation temperature peaks at
    :math:`M = 1`, so a given inlet state can absorb only a finite heat.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    References
    ----------
    Lord Rayleigh, "Aerial plane waves of finite amplitude," *Proc. R. Soc.
    Lond. A* 84 (1910) 247-284.

    Examples
    --------
    >>> state = rayleigh_flow(2.0)
    >>> round(float(state.stagnation_temperature_ratio), 4), round(float(state.pressure_ratio), 4)
    (0.7934, 0.3636)
    """
    _positive(M=M)
    check_gamma(gamma)
    M = np.asarray(M, float)
    m2 = M**2
    p = (1 + gamma) / (1 + gamma * m2)
    k = (2 + (gamma - 1) * m2) / (gamma + 1)
    return DuctFlowState(
        kind="rayleigh",
        mach=M,
        temperature_ratio=m2 * p**2,
        pressure_ratio=p,
        density_ratio=1 / (m2 * p),
        stagnation_pressure_ratio=p * k ** (gamma / (gamma - 1)),
        stagnation_temperature_ratio=m2 * p**2 * k,
        friction_parameter=np.zeros_like(M),
    )


def rayleigh_choking_heat(M, T0, gamma=1.4, R=R_AIR):
    r"""Heat per unit mass :math:`q^* = c_p (T_0^* - T_0)` that brings Rayleigh flow to :math:`M = 1`.

    Parameters
    ----------
    M : float or array_like
        Inlet Mach number [-], positive.
    T0 : float or array_like
        Inlet stagnation temperature [K], positive.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.
    R : float, optional
        Specific gas constant [J/(kg K)], positive.

    Returns
    -------
    float or ndarray
        Heat addition to choke the duct [J/kg]; adding more forces the inlet
        state to change (thermal choking).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(rayleigh_choking_heat(0.2, 300.0)) / 1e3, 1)
    1435.3
    """
    _positive(T0=T0, R=R)
    ratio = rayleigh_flow(M, gamma).stagnation_temperature_ratio
    cp = gamma * R / (gamma - 1)
    return cp * np.asarray(T0) * (1 / ratio - 1)
