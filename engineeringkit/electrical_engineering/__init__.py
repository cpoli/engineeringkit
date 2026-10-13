"""Electrical engineering fundamentals from Ohm/Kirchhoff through AC and transmission lines."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive


def ohms_law(voltage=None, current=None, resistance=None):
    """Ohms law.

    Parameters
    ----------
    voltage : V
        Input in the units stated above.
    current : A
        Input in the units stated above.
    resistance : ohm
        Input in the units stated above.

    Returns
    -------
    result
        missing voltage [V], current [A], or resistance [ohm].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    vals = [voltage, current, resistance]
    if sum(v is None for v in vals) != 1:
        raise ValueError("provide exactly two quantities")
    for name, value in zip(("voltage", "current", "resistance"), vals):
        if value is not None:
            _finite(**{name: value})
    if resistance is not None:
        _positive(resistance=resistance)
    if resistance is None and current == 0:
        raise ValueError("nonzero current is needed to infer resistance")
    if voltage is None:
        return current * resistance
    if current is None:
        return voltage / resistance
    result = voltage / current
    _nonnegative(resistance=result)
    return result


def series_resistance(*r):
    """Series resistance.

    Parameters
    ----------
    r : ohm
        Input in the units stated above.

    Returns
    -------
    result
        ohm.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _nonnegative(r=r)
    return sum(r)


def parallel_resistance(*r):
    """Parallel resistance.

    Parameters
    ----------
    r : ohm
        Input in the units stated above.

    Returns
    -------
    result
        ohm.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _positive(r=r)
    return 1 / sum(1 / x for x in r)


def capacitor_impedance(C, omega):
    """Capacitor impedance.

    Parameters
    ----------
    C : F
        Input in the units stated above.
    omega : rad/s
        Input in the units stated above.

    Returns
    -------
    result
        complex ohm.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(C=C, omega=omega)
    _positive(C=C, omega=omega)
    return 1 / (1j * omega * C)


def inductor_impedance(L, omega):
    """Inductor impedance.

    Parameters
    ----------
    L : H
        Input in the units stated above.
    omega : rad/s
        Input in the units stated above.

    Returns
    -------
    result
        complex ohm.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(L=L, omega=omega)
    _positive(L=L)
    _nonnegative(omega=omega)
    return 1j * omega * L


def series_rlc_impedance(R, L, C, omega):
    """Series rlc impedance.

    Parameters
    ----------
    R : ohm
        Input in the units stated above.
    L : H
        Input in the units stated above.
    C : F
        Input in the units stated above.
    omega : rad/s
        Input in the units stated above.

    Returns
    -------
    result
        complex ohm.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(R=R, L=L, C=C, omega=omega)
    _nonnegative(R=R)
    return R + inductor_impedance(L, omega) + capacitor_impedance(C, omega)


def resonance_frequency(L, C):
    """Resonance frequency.

    Parameters
    ----------
    L : H
        Input in the units stated above.
    C : F
        Input in the units stated above.

    Returns
    -------
    result
        Hz.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(L=L, C=C)
    _positive(L=L, C=C)
    return 1 / (2 * np.pi * np.sqrt(L * C))


def rms_sine(amplitude):
    """Rms sine.

    Parameters
    ----------
    amplitude : V
        Input in the units stated above.

    Returns
    -------
    result
        same amplitude units as input.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(amplitude=amplitude)
    _nonnegative(amplitude=amplitude)
    return amplitude / np.sqrt(2)


def real_power(v_rms, i_rms, pf=1):
    """Real power.

    Parameters
    ----------
    v_rms : V
        Input in the units stated above.
    i_rms : A
        Input in the units stated above.
    pf : 1
        Input in the units stated above.

    Returns
    -------
    result
        W.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(v_rms=v_rms, i_rms=i_rms, pf=pf)
    _nonnegative(v_rms=v_rms, i_rms=i_rms)
    _interval("pf", pf, -1, 1)
    return v_rms * i_rms * pf


def three_phase_power(v_line, i_line, pf=1):
    """Three phase power.

    Parameters
    ----------
    v_line : V
        Input in the units stated above.
    i_line : A
        Input in the units stated above.
    pf : 1
        Input in the units stated above.

    Returns
    -------
    result
        W.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(v_line=v_line, i_line=i_line, pf=pf)
    _nonnegative(v_line=v_line, i_line=i_line)
    _interval("pf", pf, -1, 1)
    return np.sqrt(3) * v_line * i_line * pf


def transformer_secondary(v_primary, n_primary, n_secondary):
    """Transformer secondary.

    Parameters
    ----------
    v_primary : V
        Input in the units stated above.
    n_primary : 1
        Input in the units stated above.
    n_secondary : 1
        Input in the units stated above.

    Returns
    -------
    result
        V.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(v_primary=v_primary, n_primary=n_primary, n_secondary=n_secondary)
    _positive(n_primary=n_primary, n_secondary=n_secondary)
    return v_primary * n_secondary / n_primary


def thevenin_from_source(vs, rs):
    """Thevenin from source.

    Parameters
    ----------
    vs : V
        Input in the units stated above.
    rs : ohm
        Input in the units stated above.

    Returns
    -------
    result
        (voltage [V], resistance [ohm]).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(vs=vs, rs=rs)
    _nonnegative(rs=rs)
    return float(vs), float(rs)


def norton_from_thevenin(vth, rth):
    """Norton from thevenin.

    Parameters
    ----------
    vth : V
        Input in the units stated above.
    rth : ohm
        Input in the units stated above.

    Returns
    -------
    result
        (current [A], resistance [ohm]).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(vth=vth, rth=rth)
    _positive(rth=rth)
    return vth / rth, rth


def wheatstone_balanced(r1, r2, r3, r4, rtol=1e-9):
    """Wheatstone balanced.

    Parameters
    ----------
    r1 : ohm
        Input in the units stated above.
    r2 : ohm
        Input in the units stated above.
    r3 : ohm
        Input in the units stated above.
    r4 : ohm
        Input in the units stated above.
    rtol : 1
        Input in the units stated above.

    Returns
    -------
    result
        bool.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(r1=r1, r2=r2, r3=r3, r4=r4, rtol=rtol)
    _positive(r1=r1, r2=r2, r3=r3, r4=r4)
    _nonnegative(rtol=rtol)
    return bool(np.isclose(r1 / r2, r3 / r4, rtol=rtol))


def transmission_line_z0(L_per_m, C_per_m):
    """Transmission line z0.

    Parameters
    ----------
    L_per_m : H/m
        Input in the units stated above.
    C_per_m : F/m
        Input in the units stated above.

    Returns
    -------
    result
        ohm.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(L_per_m=L_per_m, C_per_m=C_per_m)
    _positive(L_per_m=L_per_m, C_per_m=C_per_m)
    return np.sqrt(L_per_m / C_per_m)


def transmission_line_velocity(L_per_m, C_per_m):
    """Transmission line velocity.

    Parameters
    ----------
    L_per_m : H/m
        Input in the units stated above.
    C_per_m : F/m
        Input in the units stated above.

    Returns
    -------
    result
        m/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(L_per_m=L_per_m, C_per_m=C_per_m)
    _positive(L_per_m=L_per_m, C_per_m=C_per_m)
    return 1 / np.sqrt(L_per_m * C_per_m)


@dataclass(frozen=True)
class SeriesRLC:
    resistance: float
    inductance: float
    capacitance: float

    def __post_init__(self):
        _nonnegative(resistance=self.resistance)
        _positive(inductance=self.inductance, capacitance=self.capacitance)

    def transfer_function(self):
        """Current/voltage transfer I(s)/V(s); SI output A/V, zero initial energy."""
        return signal.TransferFunction([1, 0], [self.inductance, self.resistance, 1 / self.capacitance])

    def step_response(self, T=None):
        """Current [A] for a 1 V step, initially uncharged/unenergized."""
        if T is not None:
            _nonnegative(T=T)
            if np.asarray(T).ndim != 1 or len(T) < 2 or np.any(np.diff(T) <= 0):
                raise ValueError("T must be a strictly increasing time vector")
        t, i = signal.step(self.transfer_function(), T=T)
        return CircuitResult(t, i)


@dataclass(frozen=True)
class CircuitResult:
    t: np.ndarray
    response: np.ndarray

    def summary(self):
        return {"peak": float(np.max(self.response)), "final": float(self.response[-1])}

    def plot(self):
        fig, ax = plt.subplots()
        ax.plot(self.t, self.response)
        ax.set(xlabel="time [s]", ylabel="current [A]")
        return fig, ax


def capacitors_parallel(*c):
    """Capacitors parallel.

    Parameters
    ----------
    c : F
        Input in the units stated above.

    Returns
    -------
    result
        F.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _positive(c=c)
    return sum(c)


def capacitors_series(*c):
    """Capacitors series.

    Parameters
    ----------
    c : F
        Input in the units stated above.

    Returns
    -------
    result
        F.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _positive(c=c)
    return 1 / sum(1 / x for x in c)


def rc_time_constant(R, C):
    """Rc time constant.

    Parameters
    ----------
    R : ohm
        Input in the units stated above.
    C : F
        Input in the units stated above.

    Returns
    -------
    result
        s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(R=R, C=C)
    _positive(R=R, C=C)
    return R * C


def rl_time_constant(L, R):
    """Rl time constant.

    Parameters
    ----------
    L : H
        Input in the units stated above.
    R : ohm
        Input in the units stated above.

    Returns
    -------
    result
        s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(L=L, R=R)
    _positive(L=L, R=R)
    return L / R


def rc_charging_voltage(v_source, R, C, t, v0=0.0):
    """Rc charging voltage.

    Parameters
    ----------
    v_source : V
        Input in the units stated above.
    R : ohm
        Input in the units stated above.
    C : F
        Input in the units stated above.
    t : s
        Input in the units stated above.
    v0 : V
        Input in the units stated above.

    Returns
    -------
    result
        V.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(v_source=v_source, R=R, C=C, t=t, v0=v0)
    _nonnegative(t=t)
    t = np.asarray(t, float)
    tau = rc_time_constant(R, C)
    return v_source + (v0 - v_source) * np.exp(-t / tau)


def maximum_power_transfer(vth, rth):
    """Matched resistive load and maximum delivered power for a Thévenin source.

    Parameters
    ----------
    vth : V
        Input in the units stated above.
    rth : ohm
        Input in the units stated above.

    Returns
    -------
    result
        (matched resistance [ohm], power [W]).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(vth=vth, rth=rth)
    _positive(rth=rth)
    return rth, vth * vth / (4 * rth)


def power_factor(real_power_value, apparent_power):
    """Power factor.

    Parameters
    ----------
    real_power_value : W
        Input in the units stated above.
    apparent_power : VA
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(real_power_value=real_power_value, apparent_power=apparent_power)
    _positive(apparent_power=apparent_power)
    if abs(real_power_value) > apparent_power:
        raise ValueError("real power magnitude cannot exceed apparent power")
    return real_power_value / apparent_power


def reactive_power(v_rms, i_rms, phase_rad):
    """Reactive power.

    Parameters
    ----------
    v_rms : V
        Input in the units stated above.
    i_rms : A
        Input in the units stated above.
    phase_rad : rad
        Input in the units stated above.

    Returns
    -------
    result
        var.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/electrical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(v_rms=v_rms, i_rms=i_rms, phase_rad=phase_rad)
    _nonnegative(v_rms=v_rms, i_rms=i_rms)
    return v_rms * i_rms * np.sin(phase_rad)
