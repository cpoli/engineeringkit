r"""Thermal noise and the noise figure of cascaded amplifiers.

Johnson measured, and Nyquist explained, the noise power :math:`kTB`
delivered by any resistor; Friis showed that in a chain of amplifiers the
first stage dominates the total noise factor.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.communications.utils.decibels import db

__all__ = ["BOLTZMANN", "thermal_noise_power", "noise_figure", "cascaded_noise_factor"]

BOLTZMANN = 1.380649e-23
"""Boltzmann constant [J/K] (exact)."""


def thermal_noise_power(B, T=290.0, k=BOLTZMANN):
    r"""Available thermal (Johnson-Nyquist) noise power :math:`N = kTB`.

    Parameters
    ----------
    B : float or array_like
        Noise bandwidth [Hz], nonnegative.
    T : float or array_like, optional
        Absolute noise temperature [K], nonnegative; 290 K is the standard
        reference.
    k : float, optional
        Boltzmann constant [J/K].

    Returns
    -------
    float or ndarray
        Noise power [W] (-174 dBm in 1 Hz at 290 K).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    J. B. Johnson, "Thermal agitation of electricity in conductors," *Phys.
    Rev.* 32 (1928) 97-109; H. Nyquist, *Phys. Rev.* 32 (1928) 110-113.

    Examples
    --------
    >>> round(float(10 * np.log10(thermal_noise_power(1.0) / 1e-3)), 2)
    -173.98
    """
    _nonnegative(B=B, T=T)
    _positive(k=k)
    return k * np.asarray(T) * B


def noise_figure(F):
    r"""Noise figure :math:`NF = 10\log_{10} F` from a linear noise factor.

    Parameters
    ----------
    F : float or array_like
        Noise factor [-], at least one (one for a noiseless device).

    Returns
    -------
    float or ndarray
        Noise figure [dB].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(noise_figure(2.0)) > 3.0
    True
    """
    _interval("noise factor", F, 1, np.inf, closed="left")
    return db(F)


def cascaded_noise_factor(gains, noise_factors):
    r"""Friis's formula for the noise factor of amplifiers in cascade.

    .. math::

        F = F_1 + \frac{F_2 - 1}{G_1} + \frac{F_3 - 1}{G_1 G_2} + \cdots

    Parameters
    ----------
    gains : array_like
        Linear available power gain of each stage [-], positive. The last
        stage's gain does not affect the result.
    noise_factors : array_like
        Linear noise factor of each stage [-], at least one.

    Returns
    -------
    float
        Overall noise factor [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    H. T. Friis, "Noise figures of radio receivers," *Proc. IRE* 32 (1944)
    419-422.

    Examples
    --------
    >>> round(cascaded_noise_factor([10.0, 10.0], [2.0, 3.0]), 6)
    2.2
    """
    gains = _vector("gains", gains)
    noise_factors = _vector("noise_factors", noise_factors)
    if gains.shape != noise_factors.shape:
        raise ValueError("one gain and noise factor are required per stage")
    _positive(gains=gains)
    _interval("noise factors", noise_factors, 1, np.inf, closed="left")
    upstream_gain = np.concatenate(([1.0], np.cumprod(gains[:-1])))
    return float(noise_factors[0] + np.sum((noise_factors[1:] - 1) / upstream_gain[1:]))
