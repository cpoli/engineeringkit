r"""Heat exchangers: the log-mean temperature difference and effectiveness-NTU.

The LMTD method sizes an exchanger when all four terminal temperatures are
known. Kays and London (1955) organised compact-exchanger data around the
effectiveness :math:`\varepsilon = \dot Q/\dot Q_{max}` as a function of the
number of transfer units :math:`NTU = UA/C_{min}` and the capacity ratio
:math:`C_r = C_{min}/C_{max}`, which rates an exchanger when only the inlet
temperatures are known. Bowman, Mueller and Nagle (1940) collected the
correction factor :math:`F` that adapts the counterflow LMTD to shell-and-tube
and crossflow arrangements.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "log_mean_temperature_difference",
    "heat_exchanger_rate",
    "exchanger_effectiveness_parallel",
    "exchanger_effectiveness_counterflow",
    "exchanger_effectiveness_crossflow",
    "exchanger_effectiveness_shell_tube",
    "lmtd_correction_factor",
]


def log_mean_temperature_difference(dt1, dt2):
    r"""Log-mean temperature difference :math:`(\Delta T_1 - \Delta T_2)/\ln(\Delta T_1/\Delta T_2)`.

    Parameters
    ----------
    dt1, dt2 : float
        Hot-minus-cold temperature differences at the two ends [K], positive.

    Returns
    -------
    float
        LMTD [K]; equal to ``dt1`` when the two ends match.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Evaluated with ``log1p`` so that nearly
    equal differences keep full precision.

    Examples
    --------
    >>> round(float(log_mean_temperature_difference(60.0, 20.0)), 3)
    36.41
    """
    _finite(dt1=dt1, dt2=dt2)
    _positive(dt1=dt1, dt2=dt2)
    dt1, dt2 = np.asarray(dt1, dtype=float), np.asarray(dt2, dtype=float)
    equal = dt1 == dt2
    relative = np.where(equal, 1.0, (dt1 - dt2) / dt2)
    return np.where(equal, dt1, dt2 * relative / np.log1p(relative))[()]


def heat_exchanger_rate(U, A, dt_lm, F=1.0):
    r"""Heat-exchanger duty from the LMTD method, :math:`\dot Q = UAF\,\Delta T_{lm}`.

    Parameters
    ----------
    U : float
        Overall heat-transfer coefficient [W/(m^2 K)], positive.
    A : float
        Heat-transfer area [m^2], positive.
    dt_lm : float
        Counterflow log-mean temperature difference [K], positive.
    F : float, optional
        Correction factor for multipass or crossflow arrangements [-], in [0, 1].

    Returns
    -------
    float
        Heat rate [W].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> heat_exchanger_rate(500.0, 2.0, 36.41)
    36410.0
    """
    _finite(U=U, A=A, dt_lm=dt_lm, F=F)
    _positive(U=U, A=A, dt_lm=dt_lm)
    _interval("F", F, 0, 1)
    return U * A * F * dt_lm


def exchanger_effectiveness_parallel(NTU, Cr):
    r"""Effectiveness of a parallel-flow exchanger, :math:`\varepsilon = [1 - e^{-NTU(1+C_r)}]/(1 + C_r)`.

    Parameters
    ----------
    NTU : float or array_like
        Number of transfer units :math:`UA/C_{min}` [-], nonnegative.
    Cr : float
        Capacity ratio :math:`C_{min}/C_{max}` [-], in [0, 1].

    Returns
    -------
    float or ndarray
        Effectiveness :math:`\dot Q/[C_{min}(T_{h,in} - T_{c,in})]` [-];
        it never exceeds :math:`1/(1 + C_r)`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    References
    ----------
    W. M. Kays and A. L. London, *Compact Heat Exchangers* (1955).

    Examples
    --------
    >>> round(float(exchanger_effectiveness_parallel(2.0, 0.5)), 4)
    0.6335
    """
    _finite(NTU=NTU, Cr=Cr)
    _nonnegative(NTU=NTU)
    _interval("Cr", Cr, 0, 1)
    return -np.expm1(-np.asarray(NTU, dtype=float) * (1 + Cr)) / (1 + Cr)


def exchanger_effectiveness_counterflow(NTU, Cr):
    r"""Effectiveness of a counterflow exchanger.

    .. math::

        \varepsilon = \frac{1 - e^{-NTU(1-C_r)}}{1 - C_r e^{-NTU(1-C_r)}},
        \qquad \varepsilon = \frac{NTU}{1 + NTU} \text{ at } C_r = 1.

    Parameters
    ----------
    NTU : float
        Number of transfer units :math:`UA/C_{min}` [-], nonnegative.
    Cr : float
        Capacity ratio :math:`C_{min}/C_{max}` [-], in [0, 1].

    Returns
    -------
    float
        Effectiveness [-]; tends to one as NTU grows, for every :math:`C_r`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    References
    ----------
    W. M. Kays and A. L. London, *Compact Heat Exchangers* (1955).

    Examples
    --------
    >>> round(float(exchanger_effectiveness_counterflow(2.0, 0.5)), 4)
    0.7746
    """
    _finite(NTU=NTU, Cr=Cr)
    _nonnegative(NTU=NTU)
    _interval("Cr", Cr, 0, 1)
    NTU, Cr = np.asarray(NTU, dtype=float), np.asarray(Cr, dtype=float)
    balanced = Cr == 1
    c = np.where(balanced, 0.0, Cr)
    numerator = -np.expm1(-NTU * (1 - c))
    return np.where(balanced, NTU / (1 + NTU), numerator / (1 - c + c * numerator))[()]


_MIXED = ("none", "cmax", "cmin")


def exchanger_effectiveness_crossflow(NTU, Cr, mixed="none"):
    r"""Effectiveness of a single-pass crossflow exchanger.

    .. math::

        \varepsilon &= 1 - \exp\left[\frac{NTU^{0.22}}{C_r}\left(e^{-C_r NTU^{0.78}} - 1\right)\right]
        && \text{both fluids unmixed (approximate)},\\
        \varepsilon &= \frac{1}{C_r}\left\{1 - \exp\left[-C_r\left(1 - e^{-NTU}\right)\right]\right\}
        && C_{max} \text{ mixed}, C_{min} \text{ unmixed},\\
        \varepsilon &= 1 - \exp\left[-\frac{1 - e^{-C_r NTU}}{C_r}\right]
        && C_{min} \text{ mixed}, C_{max} \text{ unmixed}.

    Parameters
    ----------
    NTU : float or array_like
        Number of transfer units :math:`UA/C_{min}` [-], nonnegative.
    Cr : float
        Capacity ratio :math:`C_{min}/C_{max}` [-], in [0, 1].
    mixed : {"none", "cmax", "cmin"}, optional
        Which stream is mixed across the flow passage: neither (finned
        tubes, plate-fin cores), the one with the larger capacity rate, or
        the one with the smaller.

    Returns
    -------
    float or ndarray
        Effectiveness [-]; every arrangement reduces to :math:`1 - e^{-NTU}`
        at :math:`C_r = 0`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. The unmixed-unmixed expression is the usual
    fit to the exact series solution; it is within 4% of it for
    :math:`NTU \le 5`, and within about 1% for :math:`2 \le NTU \le 5`.

    References
    ----------
    W. M. Kays and A. L. London, *Compact Heat Exchangers* (1955); F. P.
    Incropera et al., *Fundamentals of Heat and Mass Transfer*, Table 11.3.

    Examples
    --------
    >>> round(float(exchanger_effectiveness_crossflow(2.0, 0.5)), 4)
    0.7388
    """
    _finite(NTU=NTU, Cr=Cr)
    _nonnegative(NTU=NTU)
    _interval("Cr", Cr, 0, 1)
    if mixed not in _MIXED:
        raise ValueError(f"mixed must be one of {list(_MIXED)}")
    NTU = np.asarray(NTU, dtype=float)
    if Cr == 0:
        return -np.expm1(-NTU)
    if mixed == "none":
        return -np.expm1(NTU**0.22 / Cr * np.expm1(-Cr * NTU**0.78))
    if mixed == "cmax":
        return -np.expm1(Cr * np.expm1(-NTU)) / Cr
    return -np.expm1(np.expm1(-Cr * NTU) / Cr)


def exchanger_effectiveness_shell_tube(NTU, Cr, shells=1):
    r"""Effectiveness of a shell-and-tube exchanger with ``shells`` shell passes and 2, 4, ... tube passes each.

    One shell pass, with :math:`\Gamma = \sqrt{1 + C_r^2}` and :math:`N_1 = NTU/n`:

    .. math::

        \varepsilon_1 = 2\left\{1 + C_r + \Gamma\,\frac{1 + e^{-N_1\Gamma}}{1 - e^{-N_1\Gamma}}\right\}^{-1},

    and :math:`n` shells in counterflow series combine as

    .. math::

        \varepsilon = \frac{[(1 - \varepsilon_1 C_r)/(1 - \varepsilon_1)]^n - 1}
                            {[(1 - \varepsilon_1 C_r)/(1 - \varepsilon_1)]^n - C_r},
        \qquad \varepsilon = \frac{n\varepsilon_1}{1 + (n - 1)\varepsilon_1} \text{ at } C_r = 1.

    Parameters
    ----------
    NTU : float or array_like
        Total number of transfer units :math:`UA/C_{min}` [-], nonnegative.
    Cr : float
        Capacity ratio :math:`C_{min}/C_{max}` [-], in [0, 1].
    shells : int, optional
        Number of shell passes :math:`n`, at least one.

    Returns
    -------
    float or ndarray
        Effectiveness [-]; between parallel flow and counterflow, tending to
        counterflow as ``shells`` grows.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. The number of tube passes per shell does
    not matter once it is even; mixing across each pass is assumed.

    References
    ----------
    W. M. Kays and A. L. London, *Compact Heat Exchangers* (1955); F. P.
    Incropera et al., *Fundamentals of Heat and Mass Transfer*, Table 11.3.

    Examples
    --------
    >>> round(float(exchanger_effectiveness_shell_tube(2.0, 0.5)), 4)
    0.6931
    """
    _finite(NTU=NTU, Cr=Cr)
    _nonnegative(NTU=NTU)
    _interval("Cr", Cr, 0, 1)
    _integer("shells", shells, 1)
    NTU = np.asarray(NTU, dtype=float)
    gamma = np.sqrt(1 + Cr**2)
    decay = -np.expm1(-NTU / shells * gamma)
    safe = np.where(decay == 0, 1.0, decay)
    eps1 = np.where(decay == 0, 0.0, 2 * safe / ((1 + Cr) * safe + gamma * (2 - safe)))
    if shells == 1:
        return eps1[()]
    if Cr == 1:
        return (shells * eps1 / (1 + (shells - 1) * eps1))[()]
    ratio = ((1 - eps1 * Cr) / (1 - eps1)) ** shells
    return ((ratio - 1) / (ratio - Cr))[()]


def lmtd_correction_factor(P, R, shells=1):
    r"""LMTD correction factor :math:`F` of a shell-and-tube exchanger (Bowman, Mueller and Nagle).

    With the temperature effectiveness of the tube fluid :math:`P = (t_o - t_i)/(T_i - t_i)`
    and the capacity ratio :math:`R = (T_i - T_o)/(t_o - t_i) = C_t/C_s`, one
    shell pass gives

    .. math::

        F = \frac{\sqrt{R^2 + 1}}{R - 1}\,
            \frac{\ln[(1 - P)/(1 - PR)]}
                 {\ln\dfrac{2 - P(R + 1 - \sqrt{R^2 + 1})}{2 - P(R + 1 + \sqrt{R^2 + 1})}}.

    For :math:`n` shells, :math:`P` is first replaced by the per-shell value
    :math:`P_1 = (1 - X)/(R - X)` with :math:`X = [(1 - PR)/(1 - P)]^{1/n}`
    (:math:`P_1 = P/[n - P(n - 1)]` at :math:`R = 1`).

    Parameters
    ----------
    P : float
        Temperature effectiveness of the tube-side fluid [-], in [0, 1).
    R : float
        Ratio of tube-side to shell-side capacity rates [-], positive.
    shells : int, optional
        Number of shell passes, at least one; each with an even number of tube passes.

    Returns
    -------
    float
        Correction factor [-] in (0, 1], so that :math:`\dot Q = UAF\,\Delta T_{lm,cf}`
        with the counterflow LMTD. Designs below about 0.75 are usually rejected.

    Raises
    ------
    ValueError
        If ``P`` is at or beyond the largest value the arrangement can reach
        for this ``R``, where :math:`F` falls to zero (a temperature cross).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. :math:`F` is symmetric in the two
    streams: :math:`F(P, R) = F(PR, 1/R)`.

    References
    ----------
    R. A. Bowman, A. C. Mueller and W. M. Nagle, "Mean temperature
    difference in design," *Trans. ASME* 62 (1940) 283-294.

    Examples
    --------
    >>> round(float(lmtd_correction_factor(0.4, 1.5)), 4)
    0.8033
    """
    _finite(P=P, R=R)
    _interval("P", P, 0, 1, closed="left")
    _positive(R=R)
    _integer("shells", shells, 1)
    if P == 0:
        return 1.0
    cross = ValueError("P is beyond the reach of this arrangement for this R (temperature cross)")
    if P * R >= 1:
        raise cross
    if shells > 1:
        if R == 1:
            P = P / (shells - P * (shells - 1))
        else:
            X = ((1 - P * R) / (1 - P)) ** (1 / shells)
            P = (1 - X) / (R - X)
    root = np.sqrt(R**2 + 1)
    denominator = 2 - P * (R + 1 + root)
    if denominator <= 0:
        raise cross
    outer = np.log((2 - P * (R + 1 - root)) / denominator)
    if R == 1:
        return float(P * root / (1 - P) / outer)
    return float(root / (R - 1) * np.log((1 - P) / (1 - P * R)) / outer)
