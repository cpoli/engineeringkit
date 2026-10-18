"""Signal measures shared by the sampling, spectral and filter models: RMS, energy, power, decibels and SNR."""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite

__all__ = ["db20", "rms", "energy", "power", "snr_db"]


def db20(x):
    r"""Amplitude ratio in decibels, :math:`20\log_{10}|x|`.

    Parameters
    ----------
    x : float or array_like
        Amplitude ratio or gain [-]; complex values use their magnitude.

    Returns
    -------
    float or ndarray
        Ratio [dB]; zero gives ``-inf``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> float(db20(10.0))
    20.0
    """
    return 20 * np.log10(np.abs(x))


def rms(x):
    r"""Root-mean-square value :math:`\sqrt{\overline{|x|^2}}`.

    Parameters
    ----------
    x : array_like
        Samples, real or complex [signal units].

    Returns
    -------
    float
        RMS value [signal units].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> float(rms([3.0, -3.0]))
    3.0
    """
    x = np.asarray(x)
    return np.sqrt(np.mean(np.abs(x) ** 2))


def energy(x):
    r"""Discrete signal energy :math:`\sum_n |x_n|^2`.

    Parameters
    ----------
    x : array_like
        Samples, real or complex [signal units].

    Returns
    -------
    float
        Energy [signal units squared], without sample-period weighting;
        multiply by :math:`1/f_s` to approximate a continuous integral.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> energy([1.0, 2.0])
    5.0
    """
    x = np.asarray(x)
    return float(np.sum(np.abs(x) ** 2))


def power(x):
    r"""Mean signal power :math:`\overline{|x|^2}`.

    Parameters
    ----------
    x : array_like
        Samples, real or complex [signal units].

    Returns
    -------
    float
        Power [signal units squared].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> power([1.0, -1.0])
    1.0
    """
    x = np.asarray(x)
    return float(np.mean(np.abs(x) ** 2))


def snr_db(signal_values, noise):
    r"""Signal-to-noise ratio :math:`10\log_{10}(P_s/P_n)` from sample arrays.

    Parameters
    ----------
    signal_values : array_like
        Signal samples, real or complex [signal units].
    noise : array_like
        Noise samples of the same shape [signal units].

    Returns
    -------
    float
        SNR [dB]; ``inf`` for zero noise and ``-inf`` for zero signal.

    Raises
    ------
    ValueError
        If the shapes differ or both powers are zero.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> float(snr_db([1.0, -1.0], [0.1, -0.1]))
    20.0
    """
    values, noise = np.asarray(signal_values), np.asarray(noise)
    _finite(signal_values=values, noise=noise)
    if values.shape != noise.shape:
        raise ValueError("signal and noise must have matching shapes")
    sp, npower = np.mean(abs(values) ** 2), np.mean(abs(noise) ** 2)
    if sp == 0 and npower == 0:
        raise ValueError("SNR is undefined for zero signal and zero noise")
    if npower == 0:
        return float("inf")
    if sp == 0:
        return float("-inf")
    return 10 * np.log10(sp / npower)
