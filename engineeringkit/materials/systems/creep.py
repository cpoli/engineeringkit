r"""Creep: Norton's power law and the Larson-Miller parameter.

Norton measured the steady creep rate of steels at temperature and fitted it
with a power of stress, :math:`\dot\varepsilon = A\sigma^n`; the temperature
dependence is an Arrhenius factor :math:`e^{-Q/RT}`. Larson and Miller
collapsed rupture data at many temperatures onto one curve of stress against
:math:`P = T(C + \log_{10} t_r)`, so short hot tests predict long cool lives.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["arrhenius_creep_rate", "larson_miller_parameter", "larson_miller_rupture_time"]

_HOUR = 3600.0


def arrhenius_creep_rate(A, Q, R, T, sigma, n=1):
    r"""Steady-state (secondary) creep rate, Norton's law with an Arrhenius factor.

    .. math:: \dot\varepsilon = A\,\sigma^n e^{-Q/RT}

    Parameters
    ----------
    A : float or array_like
        Material constant [1/(s Pa^n)], nonnegative.
    Q : float or array_like
        Activation energy for creep [J/mol], nonnegative.
    R : float or array_like
        Gas constant [J/(mol K)], positive (8.314).
    T : float or array_like
        Absolute temperature [K], positive.
    sigma : float or array_like
        Applied stress [Pa], nonnegative.
    n : float or array_like
        Stress exponent [-], positive (1 for diffusion creep, 3 to 8 for
        dislocation creep).

    Returns
    -------
    float or ndarray
        Creep strain rate [1/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    F. H. Norton, *The Creep of Steel at High Temperatures* (McGraw-Hill,
    New York, 1929).
    """
    _finite(T=T)
    _nonnegative(A=A, Q=Q, sigma=sigma)
    _positive(R=R, T=T, n=n)
    return A * np.asarray(sigma) ** n * np.exp(-np.asarray(Q) / (R * np.asarray(T)))


def larson_miller_parameter(T, rupture_time, C=20.0):
    r"""Larson-Miller parameter :math:`P = T\,(C + \log_{10} t_r)` with :math:`t_r` in hours.

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K], positive.
    rupture_time : float or array_like
        Time to rupture [s], positive; converted to hours inside, the unit
        the constant ``C`` is defined for.
    C : float
        Larson-Miller constant [-], positive (about 20 for steels).

    Returns
    -------
    float or ndarray
        Parameter [K], in the conventional K log10(h) units.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    F. R. Larson and J. Miller, "A time-temperature relationship for rupture
    and creep stresses," *Trans. ASME* 74 (1952) 765-775.

    Examples
    --------
    >>> float(larson_miller_parameter(1000.0, 3600.0 * 1000))
    23000.0
    """
    _positive(T=T, rupture_time=rupture_time, C=C)
    return np.asarray(T) * (C + np.log10(np.asarray(rupture_time) / _HOUR))


def larson_miller_rupture_time(P, T, C=20.0):
    r"""Time to rupture at temperature :math:`T` for a given Larson-Miller parameter.

    .. math:: t_r = 10^{P/T - C}\ \text{hours}

    Parameters
    ----------
    P : float or array_like
        Larson-Miller parameter [K], positive, read from a master curve at
        the working stress.
    T : float or array_like
        Absolute temperature [K], positive.
    C : float
        Larson-Miller constant [-], positive.

    Returns
    -------
    float or ndarray
        Time to rupture [s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    _positive(P=P, T=T, C=C)
    return _HOUR * 10.0 ** (np.asarray(P) / T - C)
