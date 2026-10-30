r"""Refrigerators and heat pumps: coefficients of performance and their Carnot limits.

Perkins's vapour-compression machine pumps heat from a cold space; Kelvin
saw that the same machine run against the weather heats a building with
less work than the heat it delivers. Both are bounded by the reversed
Carnot cycle, whose two COPs differ by exactly one. The vapour-compression
cycle itself is run on a Peng-Robinson refrigerant.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.thermodynamics.core.base import VaporCompressionResult
from engineeringkit.thermodynamics.systems.peng_robinson import R134A, refrigerant_properties
from engineeringkit.thermodynamics.systems.steady_flow import compression_outlet_enthalpy

__all__ = ["carnot_cop_refrigerator", "carnot_cop_heat_pump", "refrigeration_cop", "heat_pump_cop", "vapor_compression_cycle"]


def carnot_cop_refrigerator(T_hot, T_cold):
    r"""Coefficient of performance of a reversible refrigerator, :math:`T_c/(T_h - T_c)`.

    Parameters
    ----------
    T_hot, T_cold : float
        Absolute temperatures of the heat sink and the cooled space [K], with
        :math:`T_h > T_c > 0`.

    Returns
    -------
    float
        Maximum heat removed per unit work [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.
    """
    _finite(T_hot=T_hot, T_cold=T_cold)
    _positive(T_hot=T_hot, T_cold=T_cold)
    if np.any(np.asarray(T_hot) <= T_cold):
        raise ValueError("hot reservoir must be hotter")
    return T_cold / (T_hot - T_cold)


def carnot_cop_heat_pump(T_hot, T_cold):
    r"""Coefficient of performance of a reversible heat pump, :math:`T_h/(T_h - T_c)`.

    Parameters
    ----------
    T_hot, T_cold : float
        Absolute temperatures of the heated space and the heat source [K],
        with :math:`T_h > T_c > 0`.

    Returns
    -------
    float
        Maximum heat delivered per unit work [-], one more than the
        refrigerator COP between the same temperatures.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    W. Thomson, "On the economy of the heating or cooling of buildings by
    means of currents of air," *Proc. Phil. Soc. Glasgow* 3 (1852) 269-272.
    """
    _finite(T_hot=T_hot, T_cold=T_cold)
    _positive(T_hot=T_hot, T_cold=T_cold)
    if np.any(np.asarray(T_hot) <= T_cold):
        raise ValueError("hot reservoir must be hotter")
    return T_hot / (T_hot - T_cold)


def refrigeration_cop(q_cold, work_input):
    r"""Refrigerator coefficient of performance, :math:`Q_c/W`.

    Parameters
    ----------
    q_cold : float
        Heat removed from the cold space [J or W], nonnegative.
    work_input : float
        Work (or power) supplied [J or W, as ``q_cold``], positive.

    Returns
    -------
    float
        COP [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    J. Perkins, "Apparatus and means for producing ice, and in cooling
    fluids," British Patent 6662 (1834).
    """
    _finite(q_cold=q_cold, work_input=work_input)
    _nonnegative(q_cold=q_cold)
    _positive(work_input=work_input)
    return q_cold / work_input


def heat_pump_cop(q_hot, work_input):
    r"""Heat-pump coefficient of performance, :math:`Q_h/W`.

    Parameters
    ----------
    q_hot : float
        Heat delivered to the warm space [J or W], at least ``work_input``.
    work_input : float
        Work (or power) supplied [J or W, as ``q_hot``], positive.

    Returns
    -------
    float
        COP [-], at least one.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.
    """
    _finite(q_hot=q_hot, work_input=work_input)
    _positive(work_input=work_input)
    if np.any(np.asarray(q_hot) < work_input):
        raise ValueError("heat-pump heat output must be >= work input")
    return q_hot / work_input


def vapor_compression_cycle(p_evaporator, p_condenser, *, fluid=R134A, superheat=0.0, subcooling=0.0, compressor_efficiency=1.0):
    r"""Vapour-compression refrigeration cycle on a Peng-Robinson refrigerant.

    Parameters
    ----------
    p_evaporator : float
        Evaporator pressure [Pa], positive and below ``p_condenser``.
    p_condenser : float
        Condenser pressure [Pa], below the fluid's critical pressure.
    fluid : CubicFluid, optional
        Refrigerant; :data:`R134A` by default.
    superheat : float, optional
        Superheat of the vapour leaving the evaporator [K], nonnegative.
    subcooling : float, optional
        Subcooling of the liquid leaving the condenser [K], nonnegative.
    compressor_efficiency : float, optional
        Isentropic efficiency of the compressor [-] in (0, 1].

    Returns
    -------
    VaporCompressionResult
        The four states (compressor inlet and outlet, condenser outlet,
        evaporator inlet), the heat absorbed and rejected and the
        compressor work per kilogram, and both COPs.

    Notes
    -----
    Perkins's cycle: the vapour is compressed, condensed, throttled through
    an expansion valve at constant enthalpy and evaporated. Heat exchangers
    have no pressure drop. Properties come from
    :func:`refrigerant_properties`. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/thermodynamics`.

    References
    ----------
    Y. A. Çengel and M. A. Boles, *Thermodynamics: An Engineering
    Approach*, 9th ed., New York: McGraw-Hill (2019), Ch. 11.

    Examples
    --------
    >>> res = vapor_compression_cycle(0.14e6, 0.8e6)
    >>> round(res.cop_refrigerator, 2)
    3.94
    """
    _positive(p_evaporator=p_evaporator, p_condenser=p_condenser)
    if p_evaporator >= p_condenser:
        raise ValueError("p_evaporator must be below p_condenser")
    _nonnegative(superheat=superheat, subcooling=subcooling)
    _interval("compressor_efficiency", compressor_efficiency, 0, 1, closed="right")
    suction = refrigerant_properties(fluid, p_evaporator, x=1.0)
    if superheat:
        suction = refrigerant_properties(fluid, p_evaporator, suction.T + superheat)
    h_s = refrigerant_properties(fluid, p_condenser, s=suction.s).h
    discharge = refrigerant_properties(fluid, p_condenser, h=compression_outlet_enthalpy(suction.h, h_s, compressor_efficiency))
    liquid = refrigerant_properties(fluid, p_condenser, x=0.0)
    if subcooling:
        liquid = refrigerant_properties(fluid, p_condenser, liquid.T - subcooling)
    valve_exit = refrigerant_properties(fluid, p_evaporator, h=liquid.h)
    return VaporCompressionResult(
        fluid=fluid,
        labels=("compressor inlet", "compressor outlet", "condenser outlet", "evaporator inlet"),
        states=(suction, discharge, liquid, valve_exit),
        q_evaporator=suction.h - valve_exit.h,
        q_condenser=discharge.h - liquid.h,
        compressor_work=discharge.h - suction.h,
    )
