r"""Heat engines: Carnot's bound, the air-standard cycles and the Rankine steam cycle.

Carnot showed that no engine between two reservoirs beats the reversible
one, :math:`\eta = 1 - T_c/T_h`. The Otto, Diesel, dual, Brayton and Stirling
cycles idealize real engines with a calorically perfect working gas; Stirling's
regenerator brings his cycle up to Carnot's bound. Rankine's cycle uses
the liquid and vapour of water; Curzon and Ahlborn asked what efficiency an
engine has when it is run for maximum power instead.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "carnot_efficiency",
    "curzon_ahlborn_efficiency",
    "heat_engine_work",
    "thermal_efficiency",
    "otto_efficiency",
    "OttoCycle",
    "diesel_efficiency",
    "dual_efficiency",
    "brayton_efficiency",
    "stirling_efficiency",
    "rankine_ideal_efficiency",
]


def carnot_efficiency(T_hot, T_cold):
    r"""Efficiency of a reversible engine between two reservoirs, :math:`1 - T_c/T_h`.

    Parameters
    ----------
    T_hot, T_cold : float
        Absolute reservoir temperatures [K], with :math:`T_h > T_c > 0`.

    Returns
    -------
    float
        Thermal efficiency [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    S. Carnot, *Réflexions sur la puissance motrice du feu*, Paris:
    Bachelier (1824).

    Examples
    --------
    >>> carnot_efficiency(600.0, 300.0)
    0.5
    """
    _finite(T_hot=T_hot, T_cold=T_cold)
    _positive(T_hot=T_hot, T_cold=T_cold)
    if np.any(np.asarray(T_hot) <= T_cold):
        raise ValueError("require T_hot > T_cold > 0")
    return 1 - T_cold / T_hot


def curzon_ahlborn_efficiency(T_hot, T_cold):
    r"""Efficiency of an endoreversible engine at maximum power, :math:`1 - \sqrt{T_c/T_h}`.

    Parameters
    ----------
    T_hot, T_cold : float or array_like
        Absolute reservoir temperatures [K], with :math:`T_h > T_c > 0`.

    Returns
    -------
    float or ndarray
        Thermal efficiency at maximum power output [-], independent of the
        two heat-transfer conductances.

    Notes
    -----
    The engine is a Carnot cycle that exchanges heat with its reservoirs
    through finite conductances, :math:`\dot Q = K\,\Delta T`. Run reversibly
    it gives Carnot efficiency at zero power; run fast it gives no work at
    all. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    F. L. Curzon and B. Ahlborn, "Efficiency of a Carnot engine at maximum
    power output," *Am. J. Phys.* 43 (1975) 22-24.

    Examples
    --------
    >>> round(float(curzon_ahlborn_efficiency(838.15, 298.15)), 3)
    0.404
    """
    _positive(T_hot=T_hot, T_cold=T_cold)
    if np.any(np.asarray(T_hot) <= np.asarray(T_cold)):
        raise ValueError("require T_hot > T_cold > 0")
    return 1 - np.sqrt(np.asarray(T_cold) / np.asarray(T_hot))


def heat_engine_work(Qin, Qout):
    r"""Net work of a heat engine, :math:`W = Q_{in} - Q_{out}`.

    Parameters
    ----------
    Qin, Qout : float
        Heat absorbed and heat rejected per cycle [J], nonnegative, with
        :math:`Q_{in} \ge Q_{out}`.

    Returns
    -------
    float
        Net work [J].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.
    """
    _finite(Qin=Qin, Qout=Qout)
    _nonnegative(Qin=Qin, Qout=Qout)
    if np.any(np.asarray(Qout) > Qin):
        raise ValueError("heat engine requires Qin >= Qout")
    return Qin - Qout


def thermal_efficiency(Qin, Qout):
    r"""Thermal efficiency from the heat flows, :math:`1 - Q_{out}/Q_{in}`.

    Parameters
    ----------
    Qin : float
        Heat absorbed [J], positive.
    Qout : float
        Heat rejected [J], nonnegative and at most ``Qin``.

    Returns
    -------
    float
        Thermal efficiency [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.
    """
    _finite(Qin=Qin, Qout=Qout)
    _positive(Qin=Qin)
    _nonnegative(Qout=Qout)
    if np.any(np.asarray(Qout) > Qin):
        raise ValueError("heat engine requires Qin >= Qout")
    return 1 - Qout / Qin


def otto_efficiency(compression_ratio, gamma=1.4):
    r"""Air-standard Otto cycle efficiency, :math:`1 - r^{1-\gamma}`.

    Parameters
    ----------
    compression_ratio : float or array_like
        Volumetric compression ratio :math:`r = V_1/V_2` [-], at least one.
    gamma : float or array_like, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Thermal efficiency [-].

    Notes
    -----
    Heat is added and rejected at constant volume. Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/thermodynamics`.

    References
    ----------
    N. A. Otto, "Gas motor engine," U.S. Patent 194,047 (1877).
    """
    _finite(compression_ratio=compression_ratio, gamma=gamma)
    _positive(compression_ratio=compression_ratio)
    if np.any(np.asarray(compression_ratio) < 1) or np.any(np.asarray(gamma) <= 1):
        raise ValueError("require r >= 1 and gamma > 1")
    return (1 - np.asarray(compression_ratio, dtype=float) ** (1 - np.asarray(gamma, dtype=float)))[()]


@dataclass(frozen=True)
class OttoCycle:
    """Air-standard Otto cycle with a given compression ratio.

    Parameters
    ----------
    compression_ratio : float
        Volumetric compression ratio [-], at least one.
    gamma : float, optional
        Ratio of specific heats [-], greater than one.
    """

    compression_ratio: float
    gamma: float = 1.4

    def __post_init__(self):
        otto_efficiency(self.compression_ratio, self.gamma)

    @property
    def efficiency(self):
        """float: Air-standard thermal efficiency [-]."""
        return otto_efficiency(self.compression_ratio, self.gamma)

    def summary(self):
        """Compression ratio, gamma and efficiency."""
        return {"compression_ratio": self.compression_ratio, "gamma": self.gamma, "efficiency": self.efficiency}


def diesel_efficiency(r, cutoff, gamma=1.4):
    r"""Air-standard Diesel cycle efficiency.

    .. math:: \eta = 1 - \frac{r_c^\gamma - 1}{\gamma\, r^{\gamma-1}(r_c - 1)}

    Parameters
    ----------
    r : float or array_like
        Compression ratio [-].
    cutoff : float or array_like
        Cutoff ratio :math:`r_c = V_3/V_2` [-], with :math:`1 \le r_c \le r`.
    gamma : float or array_like, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Thermal efficiency [-]; the Otto value at :math:`r_c = 1`.

    Notes
    -----
    Heat is added at constant pressure. ``expm1`` keeps cutoffs near one
    accurate. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    R. Diesel, *Theorie und Konstruktion eines rationellen Wärmemotors*,
    Berlin: Springer (1893).
    """
    _finite(r=r, cutoff=cutoff, gamma=gamma)
    r, cutoff, gamma = np.broadcast_arrays(*(np.asarray(v, dtype=float) for v in (r, cutoff, gamma)))
    if np.any(cutoff < 1) or np.any(cutoff > r) or np.any(gamma <= 1):
        raise ValueError("require 1 <= cutoff <= r and gamma > 1")
    # (r_c^gamma - 1) / (gamma (r_c - 1)), which tends to one (the Otto cycle) as r_c -> 1
    excess = cutoff - 1
    safe = np.where(excess == 0, 1.0, excess)
    load = np.where(excess == 0, 1.0, np.expm1(gamma * np.log1p(excess)) / (gamma * safe))
    return (1 - load * r ** (1 - gamma))[()]


def dual_efficiency(r, pressure_ratio, cutoff, gamma=1.4):
    r"""Air-standard dual (limited-pressure) cycle efficiency.

    .. math:: \eta = 1 - \frac{\alpha r_c^\gamma - 1}{r^{\gamma-1}[(\alpha - 1) + \gamma\alpha(r_c - 1)]}

    Parameters
    ----------
    r : float or array_like
        Compression ratio [-], at least one.
    pressure_ratio : float or array_like
        Pressure ratio of the constant-volume heat addition
        :math:`\alpha = p_3/p_2` [-], at least one.
    cutoff : float or array_like
        Cutoff ratio of the constant-pressure heat addition
        :math:`r_c = V_4/V_3` [-], with :math:`1 \le r_c \le r`.
    gamma : float or array_like, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Thermal efficiency [-]: the Otto value at :math:`r_c = 1` and the
        Diesel value at :math:`\alpha = 1`.

    Notes
    -----
    Heat is added first at constant volume, then at constant pressure, as
    in a fast diesel engine whose peak pressure is limited. At
    :math:`\alpha = r_c = 1` no heat is added and the Otto limit is
    returned. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    Y. A. Çengel and M. A. Boles, *Thermodynamics: An Engineering
    Approach*, 9th ed., New York: McGraw-Hill (2019), Ch. 9.

    Examples
    --------
    >>> round(float(dual_efficiency(18.0, 1.5, 1.2)), 4)
    0.6798
    """
    _finite(r=r, pressure_ratio=pressure_ratio, cutoff=cutoff, gamma=gamma)
    r, alpha, cutoff, gamma = np.broadcast_arrays(*(np.asarray(v, dtype=float) for v in (r, pressure_ratio, cutoff, gamma)))
    if np.any(alpha < 1) or np.any(cutoff < 1) or np.any(cutoff > r) or np.any(gamma <= 1):
        raise ValueError("require pressure_ratio >= 1, 1 <= cutoff <= r and gamma > 1")
    heat = (alpha - 1) + gamma * alpha * (cutoff - 1)
    none = heat == 0
    load = np.where(none, 1.0, (alpha * cutoff**gamma - 1) / np.where(none, 1.0, heat))
    return (1 - load * r ** (1 - gamma))[()]


def brayton_efficiency(pressure_ratio, gamma=1.4):
    r"""Ideal Brayton (Joule) cycle efficiency, :math:`1 - r_p^{(1-\gamma)/\gamma}`.

    Parameters
    ----------
    pressure_ratio : float or array_like
        Compressor pressure ratio :math:`r_p` [-], at least one.
    gamma : float or array_like, optional
        Ratio of specific heats [-], greater than one.

    Returns
    -------
    float or ndarray
        Thermal efficiency [-].

    Notes
    -----
    Heat is added and rejected at constant pressure; compression and
    expansion are isentropic. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/thermodynamics`.

    References
    ----------
    G. B. Brayton, "Improvement in gas-engines," U.S. Patent 125,166 (1872).
    """
    _finite(pressure_ratio=pressure_ratio, gamma=gamma)
    if np.any(np.asarray(pressure_ratio) < 1) or np.any(np.asarray(gamma) <= 1):
        raise ValueError("require pressure ratio >= 1 and gamma > 1")
    gamma = np.asarray(gamma, dtype=float)
    return (1 - np.asarray(pressure_ratio, dtype=float) ** ((1 - gamma) / gamma))[()]


def stirling_efficiency(T_hot, T_cold, compression_ratio, gamma=1.4, regenerator_effectiveness=1.0):
    r"""Ideal-gas Stirling cycle efficiency with an imperfect regenerator.

    .. math:: \eta = \frac{(T_h - T_c)\ln r}{T_h \ln r + (1 - \varepsilon)(T_h - T_c)/(\gamma - 1)}

    Parameters
    ----------
    T_hot, T_cold : float or array_like
        Temperatures of the isothermal expansion and compression [K], with
        :math:`T_h > T_c > 0`.
    compression_ratio : float or array_like
        Volume ratio :math:`r = V_{max}/V_{min}` [-], greater than one.
    gamma : float or array_like, optional
        Ratio of specific heats [-], greater than one.
    regenerator_effectiveness : float or array_like, optional
        Fraction :math:`\varepsilon` [-] in [0, 1] of the constant-volume
        heat rejected on cooling that the regenerator returns on heating.

    Returns
    -------
    float or ndarray
        Thermal efficiency [-]; Carnot's :math:`1 - T_c/T_h` for a perfect
        regenerator, whatever the compression ratio.

    Notes
    -----
    Two isotherms joined by two constant-volume processes. The work is
    :math:`R(T_h - T_c)\ln r` per unit mass; the external heat is the
    isothermal heat :math:`RT_h\ln r` plus the part of the constant-volume
    heating :math:`c_v(T_h - T_c)` that the regenerator does not supply.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    References
    ----------
    R. Stirling, "Improvements for diminishing the consumption of fuel,"
    British Patent 4081 (1816).

    Examples
    --------
    >>> round(float(stirling_efficiency(900.0, 300.0, 3.0)), 4)  # Carnot's 1 - 300/900
    0.6667
    >>> round(float(stirling_efficiency(900.0, 300.0, 3.0, regenerator_effectiveness=0.0)), 4)
    0.2649
    """
    _positive(T_hot=T_hot, T_cold=T_cold)
    _finite(compression_ratio=compression_ratio, gamma=gamma)
    _interval("regenerator_effectiveness", regenerator_effectiveness, 0, 1)
    T_hot, T_cold, r, gamma, eps = np.broadcast_arrays(
        *(np.asarray(v, dtype=float) for v in (T_hot, T_cold, compression_ratio, gamma, regenerator_effectiveness))
    )
    if np.any(T_hot <= T_cold) or np.any(r <= 1) or np.any(gamma <= 1):
        raise ValueError("require T_hot > T_cold, compression_ratio > 1 and gamma > 1")
    log_r = np.log(r)
    return ((T_hot - T_cold) * log_r / (T_hot * log_r + (1 - eps) * (T_hot - T_cold) / (gamma - 1)))[()]


def rankine_ideal_efficiency(h1, h2, h3, h4):
    r"""Ideal Rankine cycle efficiency from the four state enthalpies.

    Parameters
    ----------
    h1 : float
        Pump outlet (boiler inlet) enthalpy [J/kg].
    h2 : float
        Boiler outlet (turbine inlet) enthalpy [J/kg].
    h3 : float
        Turbine outlet (condenser inlet) enthalpy [J/kg].
    h4 : float
        Condenser outlet (pump inlet) enthalpy [J/kg].

    Returns
    -------
    float
        :math:`[(h_2 - h_3) - (h_1 - h_4)]/(h_2 - h_1)` [-].

    Notes
    -----
    The enthalpies can come from :func:`steam_properties`. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/thermodynamics`.

    References
    ----------
    W. J. M. Rankine, *A Manual of the Steam Engine and Other Prime Movers*,
    London: Griffin (1859).
    """
    _finite(h1=h1, h2=h2, h3=h3, h4=h4)
    h1, h2, h3, h4 = np.broadcast_arrays(*(np.asarray(h, dtype=float) for h in (h1, h2, h3, h4)))
    if np.any((h2 <= h1) | (h1 < h4) | (h2 < h3) | (h3 < h4) | (h2 - h3 < h1 - h4)):
        raise ValueError("enthalpies must describe a power-producing cycle")
    wt = h2 - h3
    wp = h1 - h4
    qin = h2 - h1
    return (wt - wp) / qin
