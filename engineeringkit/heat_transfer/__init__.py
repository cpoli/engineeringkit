"""Classical heat-transfer models: conduction, convection, radiation and heat exchangers."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from scipy import special

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive


def fourier_heat_rate(k, A, dT, L):
    """Fourier heat rate.

    Parameters
    ----------
    k : W/(m K)
        Input in the units stated above.
    A : m^2
        Input in the units stated above.
    dT : K (difference)
        Input in the units stated above.
    L : m
        Input in the units stated above.

    Returns
    -------
    result
        W.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(k=k, A=A, dT=dT, L=L)
    _positive(k=k, A=A, L=L)
    return k * A * dT / L


def plane_wall_resistance(L, k, A=1.0):
    """Plane wall resistance.

    Parameters
    ----------
    L : m
        Input in the units stated above.
    k : W/(m K)
        Input in the units stated above.
    A : m^2
        Input in the units stated above.

    Returns
    -------
    result
        K/W.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(L=L, k=k, A=A)
    _nonnegative(L=L)
    _positive(k=k, A=A)
    return L / (k * A)


def cylinder_resistance(r1, r2, k, L=1.0):
    """Cylinder resistance.

    Parameters
    ----------
    r1 : m
        Input in the units stated above.
    r2 : m
        Input in the units stated above.
    k : W/(m K)
        Input in the units stated above.
    L : m
        Input in the units stated above.

    Returns
    -------
    result
        K/W.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(r1=r1, r2=r2, k=k, L=L)
    _positive(r1=r1, r2=r2, k=k, L=L)
    if r2 <= r1:
        raise ValueError("outer radius must exceed inner radius")
    return np.log(r2 / r1) / (2 * np.pi * k * L)


def sphere_resistance(r1, r2, k):
    """Sphere resistance.

    Parameters
    ----------
    r1 : m
        Input in the units stated above.
    r2 : m
        Input in the units stated above.
    k : W/(m K)
        Input in the units stated above.

    Returns
    -------
    result
        K/W.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(r1=r1, r2=r2, k=k)
    _positive(r1=r1, r2=r2, k=k)
    if r2 <= r1:
        raise ValueError("outer radius must exceed inner radius")
    return (1 / r1 - 1 / r2) / (4 * np.pi * k)


def convection_heat_rate(h, A, Ts, Tinf):
    """Convection heat rate.

    Parameters
    ----------
    h : W/(m^2 K)
        Input in the units stated above.
    A : m^2
        Input in the units stated above.
    Ts : K (absolute)
        Input in the units stated above.
    Tinf : K (absolute)
        Input in the units stated above.

    Returns
    -------
    result
        W.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(h=h, A=A, Ts=Ts, Tinf=Tinf)
    _nonnegative(h=h, Ts=Ts, Tinf=Tinf)
    _positive(A=A)
    return h * A * (Ts - Tinf)


def convection_resistance(h, A=1.0):
    """Convection resistance.

    Parameters
    ----------
    h : W/(m^2 K)
        Input in the units stated above.
    A : m^2
        Input in the units stated above.

    Returns
    -------
    result
        K/W.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(h=h, A=A)
    _positive(h=h, A=A)
    return 1 / (h * A)


def radiation_heat_rate(epsilon, A, Ts, Tsur, sigma=5.670374419e-8):
    """Radiation heat rate.

    Parameters
    ----------
    epsilon : 1
        Input in the units stated above.
    A : m^2
        Input in the units stated above.
    Ts : K (absolute)
        Input in the units stated above.
    Tsur : K (absolute)
        Input in the units stated above.
    sigma : W/(m^2 K^4)
        Input in the units stated above.

    Returns
    -------
    result
        W.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(epsilon=epsilon, A=A, Ts=Ts, Tsur=Tsur, sigma=sigma)
    _interval("epsilon", epsilon, 0, 1)
    _positive(A=A, sigma=sigma)
    _nonnegative(Ts=Ts, Tsur=Tsur)
    return epsilon * sigma * A * (Ts**4 - Tsur**4)


def radiation_coefficient(epsilon, Ts, Tsur, sigma=5.670374419e-8):
    """Radiation coefficient.

    Parameters
    ----------
    epsilon : 1
        Input in the units stated above.
    Ts : K (absolute)
        Input in the units stated above.
    Tsur : K (absolute)
        Input in the units stated above.
    sigma : W/(m^2 K^4)
        Input in the units stated above.

    Returns
    -------
    result
        W/(m^2 K).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(epsilon=epsilon, Ts=Ts, Tsur=Tsur, sigma=sigma)
    _interval("epsilon", epsilon, 0, 1)
    _positive(sigma=sigma)
    _nonnegative(Ts=Ts, Tsur=Tsur)
    return epsilon * sigma * (Ts + Tsur) * (Ts**2 + Tsur**2)


def overall_u(*resistances):
    """Combine full thermal resistances into conductance (legacy name).

    Parameters
    ----------
    resistances : K/W (each full resistance, including area)
        Input in the units stated above.

    Returns
    -------
    result
        W/K (conductance, not area-normalized U).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _positive(resistances=resistances)
    return 1 / sum(resistances)


def critical_insulation_radius_cylinder(k, h):
    """Critical insulation radius cylinder.

    Parameters
    ----------
    k : W/(m K)
        Input in the units stated above.
    h : W/(m^2 K)
        Input in the units stated above.

    Returns
    -------
    result
        m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(k=k, h=h)
    _positive(k=k, h=h)
    return k / h


def critical_insulation_radius_sphere(k, h):
    """Critical insulation radius sphere.

    Parameters
    ----------
    k : W/(m K)
        Input in the units stated above.
    h : W/(m^2 K)
        Input in the units stated above.

    Returns
    -------
    result
        m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(k=k, h=h)
    _positive(k=k, h=h)
    return 2 * k / h


def biot_number(h, Lc, k):
    """Biot number.

    Parameters
    ----------
    h : W/(m^2 K)
        Input in the units stated above.
    Lc : m
        Input in the units stated above.
    k : W/(m K)
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(h=h, Lc=Lc, k=k)
    _nonnegative(h=h)
    _positive(Lc=Lc, k=k)
    return h * Lc / k


def fourier_number(alpha, t, Lc):
    """Fourier number.

    Parameters
    ----------
    alpha : m^2/s
        Input in the units stated above.
    t : s
        Input in the units stated above.
    Lc : m
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(alpha=alpha, t=t, Lc=Lc)
    _positive(alpha=alpha, Lc=Lc)
    _nonnegative(t=t)
    return alpha * t / Lc**2


def lumped_temperature(Ti, Tinf, h, A, rho, cp, V, t):
    """Lumped temperature.

    Parameters
    ----------
    Ti : K (absolute)
        Input in the units stated above.
    Tinf : K (absolute)
        Input in the units stated above.
    h : W/(m^2 K)
        Input in the units stated above.
    A : m^2
        Input in the units stated above.
    rho : kg/m^3
        Input in the units stated above.
    cp : J/(kg K)
        Input in the units stated above.
    V : m^3
        Input in the units stated above.
    t : s
        Input in the units stated above.

    Returns
    -------
    result
        K.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(Ti=Ti, Tinf=Tinf, h=h, A=A, rho=rho, cp=cp, V=V, t=t)
    _nonnegative(Ti=Ti, Tinf=Tinf, h=h, t=t)
    _positive(A=A, rho=rho, cp=cp, V=V)
    return Tinf + (Ti - Tinf) * np.exp(-h * A * t / (rho * cp * V))


def semi_infinite_temperature(Ti, Ts, x, alpha, t):
    """Semi infinite temperature.

    Parameters
    ----------
    Ti : K (absolute)
        Input in the units stated above.
    Ts : K (absolute)
        Input in the units stated above.
    x : m
        Input in the units stated above.
    alpha : m^2/s
        Input in the units stated above.
    t : s
        Input in the units stated above.

    Returns
    -------
    result
        K.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(Ti=Ti, Ts=Ts, x=x, alpha=alpha, t=t)
    _nonnegative(Ti=Ti, Ts=Ts, x=x)
    _positive(alpha=alpha, t=t)
    return Ts + (Ti - Ts) * special.erf(x / (2 * np.sqrt(alpha * t)))


def fin_parameter(h, P, k, Ac):
    """Fin parameter.

    Parameters
    ----------
    h : W/(m^2 K)
        Input in the units stated above.
    P : m
        Input in the units stated above.
    k : W/(m K)
        Input in the units stated above.
    Ac : m^2
        Input in the units stated above.

    Returns
    -------
    result
        1/m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(h=h, P=P, k=k, Ac=Ac)
    _nonnegative(h=h)
    _positive(P=P, k=k, Ac=Ac)
    return np.sqrt(h * P / (k * Ac))


def fin_efficiency_adiabatic(m, L):
    """Adiabatic-tip straight-fin efficiency; tends to one for m*L = 0.

    Parameters
    ----------
    m : 1/m
        Input in the units stated above.
    L : m
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    if not np.all(np.isfinite([m, L])) or m < 0 or L < 0:
        raise ValueError("m and L must be finite and nonnegative")
    z = m * L
    return 1.0 if z == 0 else np.tanh(z) / z


def fin_heat_rate_adiabatic(h, P, k, Ac, L, Tb, Tinf):
    """Fin heat rate adiabatic.

    Parameters
    ----------
    h : W/(m^2 K)
        Input in the units stated above.
    P : m
        Input in the units stated above.
    k : W/(m K)
        Input in the units stated above.
    Ac : m^2
        Input in the units stated above.
    L : m
        Input in the units stated above.
    Tb : K (absolute)
        Input in the units stated above.
    Tinf : K (absolute)
        Input in the units stated above.

    Returns
    -------
    result
        W.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(h=h, P=P, k=k, Ac=Ac, L=L, Tb=Tb, Tinf=Tinf)
    _nonnegative(L=L, Tb=Tb, Tinf=Tinf)
    m = fin_parameter(h, P, k, Ac)
    return np.sqrt(h * P * k * Ac) * (Tb - Tinf) * np.tanh(m * L)


def log_mean_temperature_difference(dt1, dt2):
    """Log mean temperature difference.

    Parameters
    ----------
    dt1 : K (difference)
        Input in the units stated above.
    dt2 : K (difference)
        Input in the units stated above.

    Returns
    -------
    result
        K.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(dt1=dt1, dt2=dt2)
    _positive(dt1=dt1, dt2=dt2)
    if dt1 == dt2:
        return dt1
    relative = (dt1 - dt2) / dt2
    return dt2 * relative / np.log1p(relative)


def heat_exchanger_rate(U, A, dt_lm, F=1.0):
    """Heat exchanger rate.

    Parameters
    ----------
    U : W/(m^2 K)
        Input in the units stated above.
    A : m^2
        Input in the units stated above.
    dt_lm : K (difference)
        Input in the units stated above.
    F : 1
        Input in the units stated above.

    Returns
    -------
    result
        W.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(U=U, A=A, dt_lm=dt_lm, F=F)
    _positive(U=U, A=A, dt_lm=dt_lm)
    _interval("F", F, 0, 1)
    return U * A * F * dt_lm


def exchanger_effectiveness_parallel(NTU, Cr):
    """Exchanger effectiveness parallel.

    Parameters
    ----------
    NTU : 1
        Input in the units stated above.
    Cr : 1
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(NTU=NTU, Cr=Cr)
    _nonnegative(NTU=NTU)
    _interval("Cr", Cr, 0, 1)
    return -np.expm1(-NTU * (1 + Cr)) / (1 + Cr)


def exchanger_effectiveness_counterflow(NTU, Cr):
    """Exchanger effectiveness counterflow.

    Parameters
    ----------
    NTU : 1
        Input in the units stated above.
    Cr : 1
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(NTU=NTU, Cr=Cr)
    _nonnegative(NTU=NTU)
    _interval("Cr", Cr, 0, 1)
    if Cr == 1:
        return NTU / (1 + NTU)
    numerator = -np.expm1(-NTU * (1 - Cr))
    return numerator / (1 - Cr + Cr * numerator)


def blackbody_peak_wavelength(T):
    """Blackbody peak wavelength.

    Parameters
    ----------
    T : K (absolute)
        Input in the units stated above.

    Returns
    -------
    result
        m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(T=T)
    _positive(T=T)
    return 2.897771955e-3 / T


def stefan_boltzmann_emissive_power(T, sigma=5.670374419e-8):
    """Stefan boltzmann emissive power.

    Parameters
    ----------
    T : K (absolute)
        Input in the units stated above.
    sigma : W/(m^2 K^4)
        Input in the units stated above.

    Returns
    -------
    result
        W/m^2.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(T=T, sigma=sigma)
    _nonnegative(T=T)
    _positive(sigma=sigma)
    return sigma * T**4


@dataclass(frozen=True)
class WallResult:
    heat_rate: float
    interface_temperatures: np.ndarray

    def summary(self):
        return {"heat_rate_W": self.heat_rate, "interfaces_K": self.interface_temperatures.tolist()}

    def plot(self):
        fig, ax = plt.subplots()
        ax.plot(np.arange(len(self.interface_temperatures)), self.interface_temperatures, marker="o")
        ax.set(xlabel="interface", ylabel="temperature [K]")
        return fig, ax


def composite_wall(T_hot, T_cold, layers, A=1.0, h_hot=None, h_cold=None):
    """Composite wall.

    Parameters
    ----------
    T_hot : K (absolute)
        Input in the units stated above.
    T_cold : K (absolute)
        Input in the units stated above.
    layers : sequence of (thickness [m], conductivity [W/(m K)])
        Input in the units stated above.
    A : m^2
        Input in the units stated above.
    h_hot : W/(m^2 K)
        Input in the units stated above.
    h_cold : W/(m^2 K)
        Input in the units stated above.

    Returns
    -------
    result
        WallResult: heat rate [W], interface temperatures [K].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/heat_transfer`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(T_hot=T_hot, T_cold=T_cold)
    _nonnegative(T_hot=T_hot, T_cold=T_cold)
    _positive(A=A)
    layers = list(layers)
    if not layers:
        raise ValueError("at least one wall layer is required")
    for L, k in layers:
        _positive(L=L, k=k)
    if h_hot is not None:
        _positive(h_hot=h_hot)
    if h_cold is not None:
        _positive(h_cold=h_cold)
    rs = []
    if h_hot is not None:
        rs.append(1 / (h_hot * A))
    rs += [L / (k * A) for L, k in layers]
    if h_cold is not None:
        rs.append(1 / (h_cold * A))
    q = (T_hot - T_cold) / sum(rs)
    temps = [T_hot]
    for r in rs:
        temps.append(temps[-1] - q * r)
    return WallResult(float(q), np.array(temps))
