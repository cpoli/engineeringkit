r"""Brittle strength: Weibull's weakest-link statistics.

A ceramic or glass part breaks at its worst flaw, so its strength scatters
and a bigger part is weaker. Weibull (1939) treated the part as a chain of
volume elements each of which can fail, which gives a survival probability
:math:`\exp[-(V/V_0)(\sigma/\sigma_0)^m]`: the Weibull distribution of
:mod:`engineeringkit.reliability` with the stress in place of time and a
scale that shrinks with volume.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.reliability import weibull_cdf, weibull_mttf

__all__ = ["weibull_failure_probability", "weibull_size_effect", "weibull_mean_strength"]


def _scale(sigma_0, m, volume, reference_volume):
    _positive(sigma_0=sigma_0, m=m, volume=volume, reference_volume=reference_volume)
    return np.asarray(sigma_0) * (np.asarray(reference_volume) / volume) ** (1 / np.asarray(m))


def weibull_failure_probability(stress, sigma_0, m, volume=1.0, reference_volume=1.0):
    r"""Probability that a uniformly stressed brittle part fails.

    .. math:: P_f = 1 - \exp\left[-\frac{V}{V_0}\left(\frac{\sigma}{\sigma_0}\right)^m\right]

    Parameters
    ----------
    stress : float or array_like
        Uniform tensile stress [Pa], nonnegative.
    sigma_0 : float or array_like
        Characteristic strength of the reference volume (63.2% fail) [Pa],
        positive.
    m : float or array_like
        Weibull modulus [-], positive (5 to 20 for ceramics, above 30 for
        metals).
    volume : float or array_like, optional
        Stressed volume :math:`V` [m^3], positive.
    reference_volume : float or array_like, optional
        Volume :math:`V_0` at which ``sigma_0`` was measured [m^3], positive.

    Returns
    -------
    float or ndarray
        Failure probability [-].

    Notes
    -----
    Evaluated as :func:`engineeringkit.reliability.weibull_cdf` with scale
    :math:`\sigma_0(V_0/V)^{1/m}`. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/materials`.

    References
    ----------
    W. Weibull, "A statistical theory of the strength of materials,"
    *Ingeniörsvetenskapsakademiens Handlingar* 151 (Stockholm, 1939).

    Examples
    --------
    >>> round(float(weibull_failure_probability(300e6, 300e6, 10)), 4)
    0.6321
    """
    _nonnegative(stress=stress)
    return weibull_cdf(stress, _scale(sigma_0, m, volume, reference_volume), m)


def weibull_size_effect(strength, volume, new_volume, m):
    r"""Strength of a different volume at the same failure probability.

    .. math:: \sigma_2 = \sigma_1\left(\frac{V_1}{V_2}\right)^{1/m}

    Parameters
    ----------
    strength : float or array_like
        Strength measured on volume ``volume`` [Pa], positive.
    volume : float or array_like
        Volume of the test pieces [m^3], positive.
    new_volume : float or array_like
        Volume of the part [m^3], positive.
    m : float or array_like
        Weibull modulus [-], positive.

    Returns
    -------
    float or ndarray
        Strength of the part [Pa]; lower when it is bigger.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    Examples
    --------
    >>> round(float(weibull_size_effect(400e6, 1e-6, 1e-3, 10)) / 1e6, 1)
    200.5
    """
    _positive(strength=strength, volume=volume, new_volume=new_volume, m=m)
    return np.asarray(strength) * (np.asarray(volume) / new_volume) ** (1 / np.asarray(m))


def weibull_mean_strength(sigma_0, m, volume=1.0, reference_volume=1.0):
    r"""Mean strength :math:`\sigma_0(V_0/V)^{1/m}\,\Gamma(1 + 1/m)`.

    Parameters
    ----------
    sigma_0 : float or array_like
        Characteristic strength of the reference volume [Pa], positive.
    m : float or array_like
        Weibull modulus [-], positive.
    volume, reference_volume : float or array_like, optional
        As in :func:`weibull_failure_probability` [m^3], positive.

    Returns
    -------
    float or ndarray
        Mean strength [Pa], from :func:`engineeringkit.reliability.weibull_mttf`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    return weibull_mttf(_scale(sigma_0, m, volume, reference_volume), m)
