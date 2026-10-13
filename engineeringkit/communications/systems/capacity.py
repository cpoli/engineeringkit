r"""Signalling rate, information and channel capacity.

Nyquist bounded the pulse rate through a bandwidth :math:`B` at :math:`2B`;
Hartley counted :math:`\log_2 M` bits per distinguishable level; Shannon
defined information as entropy and proved that every noisy channel has a
capacity below which error-free communication is possible.
"""

from __future__ import annotations

import numpy as np
from scipy import special

from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "nyquist_symbol_rate",
    "symbol_rate",
    "hartley_capacity",
    "spectral_efficiency",
    "shannon_capacity",
    "entropy_binary",
    "binary_symmetric_capacity",
]


def nyquist_symbol_rate(B):
    r"""Maximum intersymbol-interference-free pulse rate through an ideal low-pass channel, :math:`2B`.

    Parameters
    ----------
    B : float or array_like
        Channel bandwidth [Hz], nonnegative.

    Returns
    -------
    float or ndarray
        Symbol rate [baud].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    H. Nyquist, "Certain topics in telegraph transmission theory," *Trans.
    AIEE* 47 (1928) 617-644.

    Examples
    --------
    >>> float(nyquist_symbol_rate(3000.0))
    6000.0
    """
    _nonnegative(B=B)
    return 2 * np.asarray(B, dtype=float)


def symbol_rate(bit_rate, M):
    r"""Symbol (baud) rate for an ``M``-ary alphabet, :math:`R_b / \log_2 M`.

    Parameters
    ----------
    bit_rate : float or array_like
        Bit rate [bit/s], nonnegative.
    M : int
        Alphabet size, at least 2.

    Returns
    -------
    float or ndarray
        Symbol rate [baud].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(symbol_rate(9600.0, 16))
    2400.0
    """
    _nonnegative(bit_rate=bit_rate)
    _integer("M", M, 2)
    return np.asarray(bit_rate) / np.log2(M)


def hartley_capacity(B, M):
    r"""Hartley's line rate, :math:`2B\log_2 M`, for ``M`` distinguishable levels.

    Parameters
    ----------
    B : float or array_like
        Bandwidth [Hz], nonnegative.
    M : int
        Number of distinguishable signal levels, at least 2.

    Returns
    -------
    float or ndarray
        Bit rate [bit/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    R. V. L. Hartley, "Transmission of information," *Bell Syst. Tech. J.* 7
    (1928) 535-563.

    Examples
    --------
    >>> float(hartley_capacity(3000.0, 4))
    12000.0
    """
    _nonnegative(B=B)
    _integer("M", M, 2)
    return 2 * np.asarray(B) * np.log2(M)


def spectral_efficiency(bit_rate, B):
    r"""Spectral efficiency :math:`R_b/B`.

    Parameters
    ----------
    bit_rate : float or array_like
        Bit rate [bit/s], nonnegative.
    B : float or array_like
        Occupied bandwidth [Hz], positive.

    Returns
    -------
    float or ndarray
        Efficiency [bit/s/Hz].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(spectral_efficiency(1e6, 2.5e5))
    4.0
    """
    _nonnegative(bit_rate=bit_rate)
    _positive(B=B)
    return np.asarray(bit_rate) / B


def shannon_capacity(B, snr_linear):
    r"""Shannon-Hartley capacity of a band-limited AWGN channel, :math:`C = B\log_2(1 + S/N)`.

    Parameters
    ----------
    B : float or array_like
        Bandwidth [Hz], nonnegative.
    snr_linear : float or array_like
        Signal-to-noise power ratio [-] (linear), nonnegative.

    Returns
    -------
    float or ndarray
        Capacity [bit/s]: the highest rate at which arbitrarily reliable
        communication is possible with ideal coding.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    C. E. Shannon, "Communication in the presence of noise," *Proc. IRE* 37
    (1949) 10-21.

    Examples
    --------
    >>> round(float(shannon_capacity(3000.0, 1000.0)))
    29902
    """
    _nonnegative(B=B, snr_linear=snr_linear)
    return np.asarray(B) * np.log2(1 + np.asarray(snr_linear))


def entropy_binary(p):
    r"""Binary entropy :math:`H(p) = -p\log_2 p - (1-p)\log_2(1-p)`.

    Parameters
    ----------
    p : float or array_like
        Probability [-] in ``[0, 1]``.

    Returns
    -------
    float or ndarray
        Entropy [bit]; zero at ``p = 0`` and ``p = 1``, one at ``p = 0.5``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    C. E. Shannon, "A mathematical theory of communication," *Bell Syst.
    Tech. J.* 27 (1948) 379-423, 623-656.

    Examples
    --------
    >>> float(entropy_binary(0.5))
    1.0
    """
    _interval("p", p, 0, 1)
    p = np.asarray(p, dtype=float)
    h = -(special.xlogy(p, p) + special.xlogy(1 - p, 1 - p)) / np.log(2)
    return float(h) if h.ndim == 0 else h


def binary_symmetric_capacity(p):
    r"""Capacity of the binary symmetric channel, :math:`C = 1 - H(p)`.

    Parameters
    ----------
    p : float or array_like
        Crossover (bit-flip) probability [-] in ``[0, 1]``.

    Returns
    -------
    float or ndarray
        Capacity [bit per channel use].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> round(binary_symmetric_capacity(0.11), 3)
    0.5
    """
    return 1 - entropy_binary(p)
