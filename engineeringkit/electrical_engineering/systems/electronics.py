r"""Electronic circuits: ideal operational amplifiers and the junction diode.

Ragazzini, Randall and Russell named the operational amplifier: a
high-gain amplifier whose negative feedback through resistors and
capacitors makes it add, scale and integrate voltages. Shockley's theory of
the p-n junction gave the diode its exponential law, which a load line
turns into an operating point.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive

__all__ = [
    "inverting_amplifier_gain",
    "noninverting_amplifier_gain",
    "op_amp_integrator",
    "thermal_voltage",
    "shockley_diode_current",
    "diode_operating_point",
]

_BOLTZMANN = 1.380649e-23
_ELEMENTARY_CHARGE = 1.602176634e-19


def inverting_amplifier_gain(r_feedback, r_input):
    r"""Voltage gain :math:`v_o/v_i = -R_f/R_{in}` of an ideal inverting amplifier.

    The ideal op amp draws no input current and, with negative feedback,
    holds its inverting input at the grounded non-inverting one (a virtual
    ground), so the input current :math:`v_i/R_{in}` all flows through
    :math:`R_f`.

    Parameters
    ----------
    r_feedback : float or array_like
        Feedback resistance :math:`R_f` [ohm], positive.
    r_input : float or array_like
        Input resistance :math:`R_{in}` [ohm], positive.

    Returns
    -------
    float or ndarray
        Voltage gain [-], negative.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    J. R. Ragazzini, R. H. Randall and F. A. Russell, "Analysis of problems
    in dynamics by electronic circuits," *Proc. IRE* 35 (1947) 444-452.

    Examples
    --------
    >>> inverting_amplifier_gain(10e3, 1e3)
    -10.0
    """
    _finite(r_feedback=r_feedback, r_input=r_input)
    _positive(r_feedback=r_feedback, r_input=r_input)
    return -r_feedback / r_input


def noninverting_amplifier_gain(r_feedback, r_ground):
    r"""Voltage gain :math:`v_o/v_i = 1 + R_f/R_g` of an ideal non-inverting amplifier.

    Parameters
    ----------
    r_feedback : float or array_like
        Feedback resistance :math:`R_f` from the output to the inverting
        input [ohm], nonnegative (zero gives the unity-gain follower).
    r_ground : float or array_like
        Resistance :math:`R_g` from the inverting input to ground [ohm],
        positive.

    Returns
    -------
    float or ndarray
        Voltage gain [-], at least one.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> noninverting_amplifier_gain(9e3, 1e3)
    10.0
    """
    _finite(r_feedback=r_feedback, r_ground=r_ground)
    _positive(r_ground=r_ground)
    if np.any(np.asarray(r_feedback) < 0):
        raise ValueError("r_feedback must be nonnegative")
    return 1 + r_feedback / r_ground


def op_amp_integrator(v_in, t, R, C, v0=0.0):
    r"""Output :math:`v_o(t) = v_0 - \frac{1}{RC}\int_0^t v_i\,dt'` of an ideal inverting integrator.

    The virtual ground makes the input current :math:`v_i/R`, and it can
    only go into the feedback capacitor. The integral of the sampled input
    is taken by the trapezoidal rule.

    Parameters
    ----------
    v_in : array_like
        Input voltage at each sample [V].
    t : array_like
        Sample times [s], strictly increasing, same length as ``v_in``
        (at least two).
    R : float
        Input resistance [ohm], positive.
    C : float
        Feedback capacitance [F], positive.
    v0 : float, optional
        Output voltage at ``t[0]`` [V].

    Returns
    -------
    ndarray
        Output voltage [V] at each sample.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`. A real integrator saturates at
    its supply rails and drifts on its input offset; neither is modelled.

    References
    ----------
    J. R. Ragazzini, R. H. Randall and F. A. Russell, "Analysis of problems
    in dynamics by electronic circuits," *Proc. IRE* 35 (1947) 444-452.

    Examples
    --------
    >>> t = np.linspace(0, 1e-3, 11)
    >>> round(float(op_amp_integrator(np.ones(11), t, 1e3, 1e-6)[-1]), 6)
    -1.0
    """
    v_in, t = np.asarray(v_in, dtype=float), np.asarray(t, dtype=float)
    _finite(v_in=v_in, t=t, R=R, C=C, v0=v0)
    _positive(R=R, C=C)
    if t.ndim != 1 or t.size < 2 or np.any(np.diff(t) <= 0):
        raise ValueError("t must be a strictly increasing time vector of at least two samples")
    if v_in.shape != t.shape:
        raise ValueError("v_in must have one sample per time")
    return v0 - cumulative_trapezoid(v_in, t, initial=0.0) / (R * C)


def thermal_voltage(temperature=300.0):
    r"""Thermal voltage :math:`V_T = k_B T/q`.

    Parameters
    ----------
    temperature : float or array_like, optional
        Junction temperature [K], positive.

    Returns
    -------
    float or ndarray
        Thermal voltage [V]; about 25.85 mV at 300 K.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> round(thermal_voltage(300.0) * 1e3, 3)
    25.852
    """
    _finite(temperature=temperature)
    _positive(temperature=temperature)
    return _BOLTZMANN * temperature / _ELEMENTARY_CHARGE


def shockley_diode_current(v, saturation_current, ideality=1.0, temperature=300.0):
    r"""Shockley's diode equation :math:`i = I_S\left(e^{v/nV_T} - 1\right)`.

    Parameters
    ----------
    v : float or array_like
        Voltage across the diode, anode to cathode [V], at most
        :math:`40\,nV_T` (beyond that the exponential overflows any real
        device's limits).
    saturation_current : float
        Reverse saturation current :math:`I_S` [A], positive.
    ideality : float, optional
        Ideality factor :math:`n` [-], positive (1 for Shockley's ideal
        junction, up to about 2 for real silicon diodes).
    temperature : float, optional
        Junction temperature [K], positive.

    Returns
    -------
    float or ndarray
        Diode current [A]; tends to :math:`-I_S` under reverse bias.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`. Reverse breakdown, series
    resistance and high injection are not modelled.

    References
    ----------
    W. Shockley, "The theory of p-n junctions in semiconductors and p-n
    junction transistors," *Bell Syst. Tech. J.* 28 (1949) 435-489.

    Examples
    --------
    >>> round(float(shockley_diode_current(0.7, 1e-14)) * 1e3, 4)
    5.7475
    """
    _finite(v=v, saturation_current=saturation_current, ideality=ideality, temperature=temperature)
    _positive(saturation_current=saturation_current, ideality=ideality)
    scale = ideality * thermal_voltage(temperature)
    v = np.asarray(v, dtype=float)
    if np.any(v > 40 * scale):
        raise ValueError("v must be at most 40 n V_T; limit the current with a series resistance")
    return saturation_current * np.expm1(v / scale)


def diode_operating_point(v_source, resistance, saturation_current, ideality=1.0, temperature=300.0):
    r"""Operating point of a diode in series with a resistor: the load-line solution.

    The resistor fixes the load line :math:`i = (V_s - v)/R`, and the
    operating point is where it crosses the diode curve of
    :func:`shockley_diode_current`. The crossing is found by bracketed root
    finding (Brent's method), on the current for forward bias, where the
    diode voltage :math:`nV_T\ln(1 + i/I_S)` stays bounded, and on the
    voltage for reverse bias.

    Parameters
    ----------
    v_source : float
        Source voltage driving the anode through the resistor [V].
    resistance : float
        Series resistance [ohm], positive.
    saturation_current : float
        Reverse saturation current [A], positive.
    ideality : float, optional
        Ideality factor [-], positive.
    temperature : float, optional
        Junction temperature [K], positive.

    Returns
    -------
    tuple of float
        ``(v_diode, i_diode)``: diode voltage [V] and current [A].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`. All parameters are scalar-only:
    each operating point is one root find.

    Examples
    --------
    >>> v, i = diode_operating_point(5.0, 1e3, 1e-14)
    >>> round(v, 4), round(i * 1e3, 4)
    (0.6925, 4.3075)
    """
    _finite(v_source=v_source)
    _positive(resistance=resistance, saturation_current=saturation_current, ideality=ideality, temperature=temperature)
    scale = ideality * thermal_voltage(temperature)
    if v_source > 0:
        current = brentq(
            lambda i: v_source - i * resistance - scale * np.log1p(i / saturation_current), 0.0, v_source / resistance, xtol=1e-18
        )
        return float(scale * np.log1p(current / saturation_current)), float(current)
    if v_source < 0:
        voltage = brentq(lambda v: saturation_current * np.expm1(v / scale) - (v_source - v) / resistance, v_source, 0.0, xtol=1e-15)
        return float(voltage), float((v_source - voltage) / resistance)
    return 0.0, 0.0
