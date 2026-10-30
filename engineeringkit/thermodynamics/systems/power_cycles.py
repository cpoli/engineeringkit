r"""Power plant cycles state by state: Rankine with reheat and feedwater heating, Brayton with regeneration.

The efficiency formulas in :mod:`~engineeringkit.thermodynamics.systems.cycles`
describe ideal cycles. A real plant is analysed device by device with the
steady-flow energy balance. Steam is superheated and reheated to keep the
turbine exhaust dry, part of it is bled to heat the feedwater, and pumps and
turbines have isentropic efficiencies. A gas turbine recovers exhaust heat in
a regenerator and divides its compression and expansion with intercoolers and
reheaters. These changes move the cycles towards Carnot's bound.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive
from engineeringkit.thermodynamics.core.base import BraytonResult, RankineResult
from engineeringkit.thermodynamics.systems.steady_flow import compression_outlet_enthalpy, expansion_outlet_enthalpy
from engineeringkit.thermodynamics.systems.steam import saturation_temperature, steam_properties
from engineeringkit.thermodynamics.utils import if97

__all__ = ["rankine_cycle", "brayton_cycle"]

_P_TRIPLE = 611.213


def _superheated(p, T, name):
    """Steam at ``p`` and ``T``, which must lie above the saturation (or critical) temperature."""
    T_min = float(saturation_temperature(p)) if p < if97.P_CRIT else if97.T_CRIT
    if T <= T_min:
        raise ValueError(f"{name} must be above {T_min:.2f} K at this pressure")
    return steam_properties(p, T)


def rankine_cycle(
    p_boiler,
    p_condenser,
    T_turbine_inlet=None,
    *,
    reheat_pressure=None,
    T_reheat=None,
    feedwater_pressure=None,
    turbine_efficiency=1.0,
    pump_efficiency=1.0,
):
    r"""Rankine steam cycle with optional superheat, reheat and an open feedwater heater.

    Parameters
    ----------
    p_boiler : float
        Boiler (turbine inlet) pressure [Pa], at most 100 MPa.
    p_condenser : float
        Condenser pressure [Pa], from 611.213 Pa to below ``p_boiler``.
    T_turbine_inlet : float, optional
        Turbine inlet temperature [K], above saturation (or above the
        critical temperature for a supercritical boiler). Omitted, the
        turbine takes saturated vapour, which needs a subcritical boiler.
    reheat_pressure : float, optional
        Pressure [Pa] at which the steam leaves the high-pressure turbine
        for the reheater, between the condenser and boiler pressures.
    T_reheat : float, optional
        Reheater outlet temperature [K]; given with ``reheat_pressure``.
    feedwater_pressure : float, optional
        Pressure [Pa] of an open (direct-contact) feedwater heater, between
        the condenser and boiler pressures and not equal to
        ``reheat_pressure``. Steam is bled from the turbine at this
        pressure; the heater delivers saturated liquid.
    turbine_efficiency, pump_efficiency : float, optional
        Isentropic efficiencies [-] in (0, 1], applied to every turbine
        section and every pump.

    Returns
    -------
    RankineResult
        State points, efficiency, turbine, pump, boiler and condenser
        energies per kilogram of boiler steam, and the bleed fraction.

    Notes
    -----
    The condenser delivers saturated liquid. The bleed fraction :math:`y`
    follows from the heater's energy balance,
    :math:`y h_{bleed} + (1 - y) h_{pump\,1} = h_f(p_{fwh})`. Each turbine
    section between two of the pressures (boiler, reheat, bleed,
    condenser) is expanded with ``turbine_efficiency`` from its own inlet
    state. Properties come from :func:`steam_properties` (IAPWS-IF97).
    Every parameter is scalar-only: one plant is solved state by state.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    Y. A. Çengel and M. A. Boles, *Thermodynamics: An Engineering
    Approach*, 9th ed., New York: McGraw-Hill (2019), Ch. 10; M. J. Moran
    et al., *Fundamentals of Engineering Thermodynamics*, 9th ed., Wiley
    (2018), Ch. 8.

    Examples
    --------
    >>> res = rankine_cycle(15e6, 10e3, 873.15, reheat_pressure=4e6, T_reheat=873.15)
    >>> round(res.efficiency, 3), round(res.turbine_exit_quality, 3)
    (0.45, 0.896)
    """
    _positive(p_boiler=p_boiler)
    _interval("p_condenser", p_condenser, _P_TRIPLE, p_boiler, closed="left")
    _interval("p_boiler", p_boiler, 0, 100e6, closed="right")
    _interval("turbine_efficiency", turbine_efficiency, 0, 1, closed="right")
    _interval("pump_efficiency", pump_efficiency, 0, 1, closed="right")
    if (reheat_pressure is None) != (T_reheat is None):
        raise ValueError("give reheat_pressure and T_reheat together")
    for name, value in (("reheat_pressure", reheat_pressure), ("feedwater_pressure", feedwater_pressure)):
        if value is not None:
            _interval(name, value, p_condenser, p_boiler, closed=False)
    if reheat_pressure is not None and reheat_pressure == feedwater_pressure:
        raise ValueError("reheat_pressure and feedwater_pressure must differ")
    if T_turbine_inlet is None:
        if p_boiler >= if97.P_CRIT:
            raise ValueError("a supercritical boiler needs T_turbine_inlet")
        inlet = steam_properties(p_boiler, x=1.0)
    else:
        inlet = _superheated(p_boiler, T_turbine_inlet, "T_turbine_inlet")

    def pump(state, p_out):
        h_s = steam_properties(p_out, s=state.s).h
        return steam_properties(p_out, h=compression_outlet_enthalpy(state.h, h_s, pump_efficiency))

    def expand(state, p_out):
        h_s = steam_properties(p_out, s=state.s).h
        return steam_properties(p_out, h=expansion_outlet_enthalpy(state.h, h_s, turbine_efficiency))

    # feed train
    labels, states = ["pump inlet"], [steam_properties(p_condenser, x=0.0)]
    if feedwater_pressure is None:
        labels.append("pump outlet")
        states.append(pump(states[0], p_boiler))
    else:
        heater = steam_properties(feedwater_pressure, x=0.0)
        labels += ["pump outlet", "heater outlet", "feed pump outlet"]
        states += [pump(states[0], feedwater_pressure), heater, pump(heater, p_boiler)]
    boiler_inlet = states[-1]
    labels.append("turbine inlet")
    states.append(inlet)

    # turbine sections, recording (flow fraction placeholder, inlet, outlet)
    stops = sorted({p for p in (reheat_pressure, feedwater_pressure) if p is not None}, reverse=True) + [p_condenser]
    state, sections, reheats = inlet, [], []
    for p in stops:
        out = expand(state, p)
        sections.append((state, out))
        if p == reheat_pressure:
            labels.append("HP turbine exit")
            states.append(out)
            state = _superheated(p, T_reheat, "T_reheat")
            reheats.append((len(sections), out, state))
            labels.append("reheat outlet")
        else:
            state = out
            labels.append("bleed" if p == feedwater_pressure else "turbine exit")
        states.append(state)

    y = 0.0
    if feedwater_pressure is not None:
        bleed = states[labels.index("bleed")]
        pumped = states[labels.index("pump outlet")]
        y = (heater.h - pumped.h) / (bleed.h - pumped.h)

    def fraction(i):
        """Flow fraction through turbine section ``i``, which starts at ``stops[i - 1]`` (the boiler for 0)."""
        return 1.0 - y if i and feedwater_pressure is not None and stops[i - 1] <= feedwater_pressure else 1.0

    turbine_work = sum(fraction(i) * (a.h - b.h) for i, (a, b) in enumerate(sections))
    heat_in = (inlet.h - boiler_inlet.h) + sum(fraction(i) * (hot.h - cold.h) for i, cold, hot in reheats)
    exit_state = states[labels.index("turbine exit")]
    heat_out = (1 - y) * (exit_state.h - states[0].h)
    if feedwater_pressure is None:
        pump_work = states[1].h - states[0].h
    else:
        pump_work = (1 - y) * (states[1].h - states[0].h) + (states[3].h - states[2].h)
    return RankineResult(
        labels=tuple(labels),
        states=tuple(states),
        efficiency=(turbine_work - pump_work) / heat_in,
        turbine_work=turbine_work,
        pump_work=pump_work,
        heat_in=heat_in,
        heat_out=heat_out,
        bleed_fraction=y,
    )


def brayton_cycle(
    pressure_ratio,
    T_min,
    T_max,
    *,
    compressor_efficiency=1.0,
    turbine_efficiency=1.0,
    regenerator_effectiveness=0.0,
    compressor_stages=1,
    turbine_stages=1,
    cp=1005.0,
    gamma=1.4,
    p_min=101325.0,
):
    r"""Air-standard gas-turbine cycle with regeneration, intercooling and reheat.

    Parameters
    ----------
    pressure_ratio : float
        Overall compressor pressure ratio [-], greater than one.
    T_min : float
        Compressor inlet temperature [K], positive; every intercooler
        returns the air to it.
    T_max : float
        Turbine inlet temperature [K], above the compressor outlet
        temperature; every reheater returns the gas to it.
    compressor_efficiency, turbine_efficiency : float, optional
        Isentropic efficiency of each stage [-] in (0, 1].
    regenerator_effectiveness : float, optional
        Regenerator effectiveness [-] in [0, 1]; zero for no regenerator.
        It heats the compressed air by
        :math:`\varepsilon(T_{exhaust} - T_{compressor})` when the exhaust
        is the hotter stream and does nothing otherwise.
    compressor_stages, turbine_stages : int, optional
        Number of compressor stages (with an intercooler between each pair)
        and turbine stages (with a reheater between each pair), each with
        the same pressure ratio :math:`r_p^{1/N}`.
    cp : float, optional
        Isobaric specific heat [J/(kg K)], positive.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.
    p_min : float, optional
        Compressor inlet pressure [Pa], positive; it only sets the pressure
        labels of the states.

    Returns
    -------
    BraytonResult
        State temperatures, pressures and entropies; efficiency; compressor
        and turbine work, heat added and heat regenerated per kilogram.

    Notes
    -----
    With ideal components and a perfect regenerator the efficiency is
    :math:`1 - (T_{min}/T_{max})\,r_p^{(\gamma-1)/\gamma}`, which falls with
    pressure ratio; with many intercooled and reheated stages it approaches
    the Ericsson cycle and Carnot's bound. Pressure losses in the heat
    exchangers are neglected. Every parameter is scalar-only: one cycle is
    built with its list of states. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/thermodynamics`.

    References
    ----------
    Y. A. Çengel and M. A. Boles, *Thermodynamics: An Engineering
    Approach*, 9th ed., New York: McGraw-Hill (2019), Ch. 9.

    Examples
    --------
    >>> res = brayton_cycle(8.0, 300.0, 1300.0, compressor_efficiency=0.8, turbine_efficiency=0.85)
    >>> round(res.efficiency, 3), round(res.back_work_ratio, 3)
    (0.274, 0.615)
    """
    _positive(T_min=T_min, T_max=T_max, cp=cp, p_min=p_min)
    _finite(pressure_ratio=pressure_ratio, gamma=gamma)
    if pressure_ratio <= 1 or gamma <= 1:
        raise ValueError("require pressure_ratio > 1 and gamma > 1")
    for name, value in (("compressor_efficiency", compressor_efficiency), ("turbine_efficiency", turbine_efficiency)):
        _interval(name, value, 0, 1, closed="right")
    _interval("regenerator_effectiveness", regenerator_effectiveness, 0, 1)
    _integer("compressor_stages", compressor_stages, 1)
    _integer("turbine_stages", turbine_stages, 1)
    k = (gamma - 1) / gamma
    T_comp = T_min + T_min * ((pressure_ratio ** (1 / compressor_stages)) ** k - 1) / compressor_efficiency
    T_turb = T_max - turbine_efficiency * T_max * (1 - (pressure_ratio ** (-1 / turbine_stages)) ** k)
    if T_max <= T_comp:
        raise ValueError("T_max must exceed the compressor outlet temperature")
    regenerated = regenerator_effectiveness * max(T_turb - T_comp, 0.0)
    T_regen = T_comp + regenerated

    labels: list[str] = []
    temperatures: list[float] = []
    pressures: list[float] = []

    def add(label, temperature, pressure):
        labels.append(label)
        temperatures.append(temperature)
        pressures.append(pressure)

    p_high = p_min * pressure_ratio
    for i in range(compressor_stages):
        p_in = p_min * pressure_ratio ** (i / compressor_stages)
        add("compressor inlet" if i == 0 else f"intercooler {i} outlet", T_min, p_in)
        add(f"compressor {i + 1} outlet", T_comp, p_in * pressure_ratio ** (1 / compressor_stages))
    if regenerated:
        add("regenerator outlet (cold side)", T_regen, p_high)
    for i in range(turbine_stages):
        p_in = p_high * pressure_ratio ** (-i / turbine_stages)
        add("turbine inlet" if i == 0 else f"reheater {i} outlet", T_max, p_in)
        add(f"turbine {i + 1} outlet", T_turb, p_in * pressure_ratio ** (-1 / turbine_stages))
    if regenerated:
        add("regenerator outlet (hot side)", T_turb - regenerated, p_min)
    T, p = np.array(temperatures), np.array(pressures)
    s = cp * np.log(T / T_min) - cp * k * np.log(p / p_min)

    compressor_work = compressor_stages * cp * (T_comp - T_min)
    turbine_work = turbine_stages * cp * (T_max - T_turb)
    heat_in = cp * (T_max - T_regen) + (turbine_stages - 1) * cp * (T_max - T_turb)
    return BraytonResult(
        labels=tuple(labels),
        T=T,
        p=p,
        s=s,
        efficiency=(turbine_work - compressor_work) / heat_in,
        compressor_work=compressor_work,
        turbine_work=turbine_work,
        heat_in=heat_in,
        heat_regenerated=cp * regenerated,
        cp=cp,
    )
