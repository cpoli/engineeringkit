"""Classical thermodynamics and idealized power/refrigeration cycles."""

from dataclasses import dataclass

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

R_UNIVERSAL = 8.31446261815324


def ideal_gas_pressure(n, T, V):
    """Ideal gas pressure.

    Parameters
    ----------
    n : mol
        Input in the units stated above.
    T : K (absolute)
        Input in the units stated above.
    V : m^3
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(n=n, T=T, V=V)
    _positive(n=n, T=T, V=V)
    return n * R_UNIVERSAL * T / V


def ideal_gas_temperature(P, V, n):
    """Ideal gas temperature.

    Parameters
    ----------
    P : Pa
        Input in the units stated above.
    V : m^3
        Input in the units stated above.
    n : mol
        Input in the units stated above.

    Returns
    -------
    result
        K.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(P=P, V=V, n=n)
    _positive(P=P, V=V, n=n)
    return P * V / (n * R_UNIVERSAL)


def first_law_heat(delta_u, work_by_system):
    """First law heat.

    Parameters
    ----------
    delta_u : J
        Input in the units stated above.
    work_by_system : J
        Input in the units stated above.

    Returns
    -------
    result
        J.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return delta_u + work_by_system


def entropy_change_isothermal_ideal_gas(n, T, V1, V2):
    """Entropy change isothermal ideal gas.

    Parameters
    ----------
    n : mol
        Input in the units stated above.
    T : K (absolute)
        Input in the units stated above.
    V1 : m^3
        Input in the units stated above.
    V2 : m^3
        Input in the units stated above.

    Returns
    -------
    result
        J/K.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(n=n, T=T, V1=V1, V2=V2)
    _positive(n=n, T=T, V1=V1, V2=V2)
    return n * R_UNIVERSAL * np.log(V2 / V1)


def carnot_efficiency(T_hot, T_cold):
    """Carnot efficiency.

    Parameters
    ----------
    T_hot : K (absolute)
        Input in the units stated above.
    T_cold : K (absolute)
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(T_hot=T_hot, T_cold=T_cold)
    _positive(T_hot=T_hot, T_cold=T_cold)
    if T_hot <= T_cold or T_cold <= 0:
        raise ValueError("require T_hot > T_cold > 0")
    return 1 - T_cold / T_hot


def carnot_cop_refrigerator(T_hot, T_cold):
    """Carnot cop refrigerator.

    Parameters
    ----------
    T_hot : K (absolute)
        Input in the units stated above.
    T_cold : K (absolute)
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(T_hot=T_hot, T_cold=T_cold)
    _positive(T_hot=T_hot, T_cold=T_cold)
    if T_hot <= T_cold:
        raise ValueError("hot reservoir must be hotter")
    return T_cold / (T_hot - T_cold)


def carnot_cop_heat_pump(T_hot, T_cold):
    """Carnot cop heat pump.

    Parameters
    ----------
    T_hot : K (absolute)
        Input in the units stated above.
    T_cold : K (absolute)
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(T_hot=T_hot, T_cold=T_cold)
    _positive(T_hot=T_hot, T_cold=T_cold)
    if T_hot <= T_cold:
        raise ValueError("hot reservoir must be hotter")
    return T_hot / (T_hot - T_cold)


def otto_efficiency(compression_ratio, gamma=1.4):
    """Otto efficiency.

    Parameters
    ----------
    compression_ratio : 1
        Input in the units stated above.
    gamma : 1
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(compression_ratio=compression_ratio, gamma=gamma)
    _positive(compression_ratio=compression_ratio)
    if compression_ratio < 1 or gamma <= 1:
        raise ValueError("require r >= 1 and gamma > 1")
    return 1 - compression_ratio ** (1 - gamma)


def diesel_efficiency(r, cutoff, gamma=1.4):
    """Diesel efficiency.

    Parameters
    ----------
    r : 1
        Input in the units stated above.
    cutoff : 1
        Input in the units stated above.
    gamma : 1
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(r=r, cutoff=cutoff, gamma=gamma)
    if not 1 <= cutoff <= r or gamma <= 1:
        raise ValueError("require 1 <= cutoff <= r and gamma > 1")
    if cutoff == 1:
        return otto_efficiency(r, gamma)
    return 1 - np.expm1(gamma * np.log(cutoff)) / (gamma * r ** (gamma - 1) * (cutoff - 1))


def brayton_efficiency(pressure_ratio, gamma=1.4):
    """Brayton efficiency.

    Parameters
    ----------
    pressure_ratio : 1
        Input in the units stated above.
    gamma : 1
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(pressure_ratio=pressure_ratio, gamma=gamma)
    if pressure_ratio < 1 or gamma <= 1:
        raise ValueError("require pressure ratio >= 1 and gamma > 1")
    return 1 - pressure_ratio ** ((1 - gamma) / gamma)


def isentropic_temperature_ratio(pressure_ratio, gamma=1.4):
    """Isentropic temperature ratio.

    Parameters
    ----------
    pressure_ratio : 1
        Input in the units stated above.
    gamma : 1
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(pressure_ratio=pressure_ratio, gamma=gamma)
    _positive(pressure_ratio=pressure_ratio)
    if gamma <= 1:
        raise ValueError("gamma must exceed 1")
    return pressure_ratio ** ((gamma - 1) / gamma)


def specific_internal_energy_ideal_gas(T, cv):
    """Specific internal energy ideal gas.

    Parameters
    ----------
    T : K (absolute)
        Input in the units stated above.
    cv : J/(kg K)
        Input in the units stated above.

    Returns
    -------
    result
        J/kg.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(T=T, cv=cv)
    _nonnegative(T=T)
    _positive(cv=cv)
    return cv * T


def specific_enthalpy_ideal_gas(T, cp):
    """Specific enthalpy ideal gas.

    Parameters
    ----------
    T : K (absolute)
        Input in the units stated above.
    cp : J/(kg K)
        Input in the units stated above.

    Returns
    -------
    result
        J/kg.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(T=T, cp=cp)
    _nonnegative(T=T)
    _positive(cp=cp)
    return cp * T


def exergy_heat(Q, T_source, T0):
    """Exergy heat.

    Parameters
    ----------
    Q : J
        Input in the units stated above.
    T_source : K (absolute)
        Input in the units stated above.
    T0 : K (absolute)
        Input in the units stated above.

    Returns
    -------
    result
        J.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(Q=Q, T_source=T_source, T0=T0)
    _positive(T_source=T_source, T0=T0)
    return Q * (1 - T0 / T_source)


def heat_engine_work(Qin, Qout):
    """Heat engine work.

    Parameters
    ----------
    Qin : J
        Input in the units stated above.
    Qout : J
        Input in the units stated above.

    Returns
    -------
    result
        J.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(Qin=Qin, Qout=Qout)
    _nonnegative(Qin=Qin, Qout=Qout)
    if Qout > Qin:
        raise ValueError("heat engine requires Qin >= Qout")
    return Qin - Qout


def thermal_efficiency(Qin, Qout):
    """Thermal efficiency.

    Parameters
    ----------
    Qin : J
        Input in the units stated above.
    Qout : J
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(Qin=Qin, Qout=Qout)
    _positive(Qin=Qin)
    _nonnegative(Qout=Qout)
    if Qout > Qin:
        raise ValueError("heat engine requires Qin >= Qout")
    return 1 - Qout / Qin


@dataclass(frozen=True)
class OttoCycle:
    compression_ratio: float
    gamma: float = 1.4

    def __post_init__(self):
        otto_efficiency(self.compression_ratio, self.gamma)

    @property
    def efficiency(self):
        return otto_efficiency(self.compression_ratio, self.gamma)

    def summary(self):
        return {"compression_ratio": self.compression_ratio, "gamma": self.gamma, "efficiency": self.efficiency}


def clausius_inequality(cyclic_heat_over_temperature):
    """Return whether a discretized cyclic Clausius integral satisfies sum(dQ/T) <= 0.

    Parameters
    ----------
    cyclic_heat_over_temperature : J/K
        Input in the units stated above.

    Returns
    -------
    result
        bool.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return float(np.sum(cyclic_heat_over_temperature)) <= 1e-12


def reversible_isothermal_work_ideal_gas(n, T, V1, V2):
    """Reversible isothermal work ideal gas.

    Parameters
    ----------
    n : mol
        Input in the units stated above.
    T : K (absolute)
        Input in the units stated above.
    V1 : m^3
        Input in the units stated above.
    V2 : m^3
        Input in the units stated above.

    Returns
    -------
    result
        J.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(n=n, T=T, V1=V1, V2=V2)
    _positive(n=n, T=T, V1=V1, V2=V2)
    return n * R_UNIVERSAL * T * np.log(V2 / V1)


def isentropic_pressure_volume(p1, v1, v2, gamma=1.4):
    """Isentropic pressure volume.

    Parameters
    ----------
    p1 : Pa
        Input in the units stated above.
    v1 : consistent volume units (e.g. m^3/kg)
        Input in the units stated above.
    v2 : consistent volume units (e.g. m^3/kg)
        Input in the units stated above.
    gamma : 1
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(p1=p1, v1=v1, v2=v2, gamma=gamma)
    _positive(p1=p1, v1=v1, v2=v2)
    if gamma <= 1:
        raise ValueError("gamma must exceed 1")
    return p1 * (v1 / v2) ** gamma


def rankine_ideal_efficiency(h1, h2, h3, h4):
    """Ideal Rankine efficiency from pump-in, boiler-out, turbine-out, condenser-out enthalpies.

    Parameters
    ----------
    h1 : J/kg
        Input in the units stated above.
    h2 : J/kg
        Input in the units stated above.
    h3 : J/kg
        Input in the units stated above.
    h4 : J/kg
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(h1=h1, h2=h2, h3=h3, h4=h4)
    if not h2 > h1 >= h4 or h2 < h3 or h3 < h4 or h2 - h3 < h1 - h4:
        raise ValueError("enthalpies must describe a power-producing cycle")
    wt = h2 - h3
    wp = h1 - h4
    qin = h2 - h1
    return (wt - wp) / qin


def refrigeration_cop(q_cold, work_input):
    """Refrigeration cop.

    Parameters
    ----------
    q_cold : J
        Input in the units stated above.
    work_input : J
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(q_cold=q_cold, work_input=work_input)
    _nonnegative(q_cold=q_cold)
    _positive(work_input=work_input)
    return q_cold / work_input


def heat_pump_cop(q_hot, work_input):
    """Heat pump cop.

    Parameters
    ----------
    q_hot : J
        Input in the units stated above.
    work_input : J
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(q_hot=q_hot, work_input=work_input)
    _positive(work_input=work_input)
    if q_hot < work_input:
        raise ValueError("heat-pump heat output must be >= work input")
    return q_hot / work_input


def exergy_destruction(T0, entropy_generation):
    """Exergy destruction.

    Parameters
    ----------
    T0 : K (absolute)
        Input in the units stated above.
    entropy_generation : J/K
        Input in the units stated above.

    Returns
    -------
    result
        J.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/thermodynamics`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(T0=T0, entropy_generation=entropy_generation)
    _positive(T0=T0)
    _nonnegative(entropy_generation=entropy_generation)
    return T0 * entropy_generation
