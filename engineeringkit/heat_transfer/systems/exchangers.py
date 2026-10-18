r"""Heat exchangers: the log-mean temperature difference and effectiveness-NTU.

The LMTD method sizes an exchanger when all four terminal temperatures are
known. Kays and London (1955) organised compact-exchanger data around the
effectiveness :math:`\varepsilon = \dot Q/\dot Q_{max}` as a function of the
number of transfer units :math:`NTU = UA/C_{min}` and the capacity ratio
:math:`C_r = C_{min}/C_{max}`, which rates an exchanger when only the inlet
temperatures are known.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "log_mean_temperature_difference",
    "heat_exchanger_rate",
    "exchanger_effectiveness_parallel",
    "exchanger_effectiveness_counterflow",
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
    if dt1 == dt2:
        return dt1
    relative = (dt1 - dt2) / dt2
    return dt2 * relative / np.log1p(relative)


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
    if Cr == 1:
        return NTU / (1 + NTU)
    numerator = -np.expm1(-NTU * (1 - Cr))
    return numerator / (1 - Cr + Cr * numerator)
