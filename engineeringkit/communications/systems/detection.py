r"""Optimum detection of antipodal signals in additive white Gaussian noise.

Kotel'nikov's theory of optimum noise immunity gives the bit-error rate of
coherently detected BPSK as :math:`Q(\sqrt{2E_b/N_0})`; Gray-labelled QPSK
achieves the same per-bit error rate at twice the spectral efficiency.
"""

from __future__ import annotations

import numpy as np
from scipy import special

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative

__all__ = ["q_function", "ber_bpsk_awgn", "ber_qpsk_awgn"]


def q_function(x):
    r"""Gaussian tail probability :math:`Q(x) = \tfrac12\operatorname{erfc}(x/\sqrt2)`.

    Parameters
    ----------
    x : float or array_like
        Threshold in standard deviations [-].

    Returns
    -------
    float or ndarray
        :math:`P(Z > x)` for a standard normal :math:`Z`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(q_function(0.0))
    0.5
    """
    _finite(x=x)
    return 0.5 * special.erfc(np.asarray(x) / np.sqrt(2))


def ber_bpsk_awgn(ebn0):
    r"""Bit-error rate of coherent BPSK in AWGN, :math:`Q(\sqrt{2E_b/N_0}) = \tfrac12\operatorname{erfc}\sqrt{E_b/N_0}`.

    Parameters
    ----------
    ebn0 : float or array_like
        Energy per bit over noise spectral density, :math:`E_b/N_0` [-]
        (linear, not dB), nonnegative.

    Returns
    -------
    float or ndarray
        Bit-error probability [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    V. A. Kotel'nikov, *The Theory of Optimum Noise Immunity* (1947; English
    translation, McGraw-Hill, 1959).

    Examples
    --------
    >>> f"{float(ber_bpsk_awgn(10 ** (9.6 / 10))):.1e}"
    '9.7e-06'
    """
    _nonnegative(ebn0=ebn0)
    return 0.5 * special.erfc(np.sqrt(ebn0))


def ber_qpsk_awgn(ebn0):
    r"""Bit-error rate of coherent Gray-labelled QPSK in AWGN; equal to BPSK per bit.

    Parameters
    ----------
    ebn0 : float or array_like
        :math:`E_b/N_0` [-] (linear), nonnegative.

    Returns
    -------
    float or ndarray
        Bit-error probability [-].

    Notes
    -----
    The in-phase and quadrature rails are independent BPSK channels. See
    :doc:`/guides/communications`.

    Examples
    --------
    >>> bool(ber_qpsk_awgn(4.0) == ber_bpsk_awgn(4.0))
    True
    """
    return ber_bpsk_awgn(ebn0)
