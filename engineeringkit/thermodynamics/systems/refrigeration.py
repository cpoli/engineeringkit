r"""Refrigerators and heat pumps: coefficients of performance and their Carnot limits.

Perkins's vapour-compression machine pumps heat from a cold space; Kelvin
saw that the same machine run against the weather heats a building with
less work than the heat it delivers. Both are bounded by the reversed
Carnot cycle, whose two COPs differ by exactly one.
"""

from __future__ import annotations

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["carnot_cop_refrigerator", "carnot_cop_heat_pump", "refrigeration_cop", "heat_pump_cop"]


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
    if T_hot <= T_cold:
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
    if T_hot <= T_cold:
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
    if q_hot < work_input:
        raise ValueError("heat-pump heat output must be >= work input")
    return q_hot / work_input
