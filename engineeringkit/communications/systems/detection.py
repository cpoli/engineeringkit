r"""Optimum detection of known signals in additive white Gaussian noise.

North's matched filter maximises the output signal-to-noise ratio at
:math:`2E/N_0`; Kotel'nikov's theory of optimum noise immunity gives the
bit-error rate of coherently detected BPSK as :math:`Q(\sqrt{2E_b/N_0})`;
Gray-labelled QPSK achieves the same per-bit error rate at twice the
spectral efficiency; and minimum-distance detection gives the exact symbol
error rates of M-PSK and square M-QAM.
"""

from __future__ import annotations

import numpy as np
from scipy import integrate, special

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "q_function",
    "ber_bpsk_awgn",
    "ber_qpsk_awgn",
    "matched_filter",
    "matched_filter_snr",
    "minimum_distance_detect",
    "ser_mpsk_awgn",
    "ser_mqam_awgn",
]


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


def _complex_vector(name, value):
    a = np.asarray(value)
    if a.ndim != 1 or not a.size:
        raise ValueError(f"{name} must be a nonempty 1-D array")
    _finite(**{name: a})
    return a


def matched_filter(pulse):
    r"""Impulse response of the filter matched to a pulse, :math:`h[n] = s^*[N - 1 - n]`.

    Parameters
    ----------
    pulse : array_like
        Samples of the known pulse :math:`s` (real or complex).

    Returns
    -------
    ndarray
        The time-reversed, conjugated pulse. Convolving the received samples
        with it and sampling at the pulse's end, sample :math:`N - 1`, gives
        the correlation :math:`\sum r[n] s^*[n]`.

    Notes
    -----
    In white noise this filter maximises the output signal-to-noise ratio
    over all linear filters; see :func:`matched_filter_snr` and
    :doc:`/guides/communications`.

    References
    ----------
    D. O. North, "An analysis of the factors which determine signal/noise
    discrimination in pulsed-carrier systems," RCA Laboratories Report
    PTR-6C (25 June 1943); reprinted in *Proc. IEEE* 51 (1963) 1016-1027.

    Examples
    --------
    >>> matched_filter([1.0, 2.0, 3.0]).tolist()
    [3.0, 2.0, 1.0]
    """
    return np.conj(_complex_vector("pulse", pulse)[::-1])


def matched_filter_snr(pulse_energy, N0):
    r"""Peak output signal-to-noise ratio of a matched filter in white noise, :math:`2E/N_0`.

    Parameters
    ----------
    pulse_energy : float or array_like
        Received pulse energy :math:`E` [J], nonnegative.
    N0 : float or array_like
        One-sided noise power spectral density [W/Hz], positive (the
        two-sided density is :math:`N_0/2`).

    Returns
    -------
    float or ndarray
        Output SNR [-] (linear). It depends only on the pulse energy, not on
        its shape.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(matched_filter_snr(2e-9, 1e-9))
    4.0
    """
    _nonnegative(pulse_energy=pulse_energy)
    _positive(N0=N0)
    return 2 * np.asarray(pulse_energy, dtype=float) / N0


def minimum_distance_detect(received, constellation):
    """Index of the constellation point nearest each received sample (maximum likelihood in AWGN).

    Parameters
    ----------
    received : array_like
        Received complex (or real) samples.
    constellation : array_like
        Candidate points, for example from :func:`psk_constellation` or
        :func:`qam_constellation`.

    Returns
    -------
    ndarray of int
        The detected label of each sample.

    Notes
    -----
    See :doc:`/guides/communications`.

    Examples
    --------
    >>> from engineeringkit.communications import psk_constellation
    >>> minimum_distance_detect([0.9 + 0.2j, -0.1 - 1.2j], psk_constellation(4)).tolist()
    [0, 2]
    """
    r = np.asarray(received)
    if not r.size:
        raise ValueError("received must be nonempty")
    _finite(received=r)
    points = _complex_vector("constellation", constellation)
    return np.argmin(np.abs(r.reshape(-1, 1) - points) ** 2, axis=1).reshape(r.shape)


def ser_mpsk_awgn(esn0, M):
    r"""Exact symbol-error rate of coherent M-PSK in AWGN, by Craig's integral.

    .. math::

        P_s = \frac{1}{\pi}\int_0^{(M-1)\pi/M}
        \exp\left(-\frac{E_s}{N_0}\frac{\sin^2(\pi/M)}{\sin^2\theta}\right)d\theta

    Parameters
    ----------
    esn0 : float or array_like
        Energy per symbol over noise spectral density, :math:`E_s/N_0` [-]
        (linear), nonnegative; :math:`E_s = E_b\log_2 M`.
    M : int
        Number of phases, at least 2.

    Returns
    -------
    float or ndarray
        Symbol-error probability [-]. At high SNR it approaches
        :math:`2Q(\sqrt{2E_s/N_0}\sin(\pi/M))`; with Gray labelling the
        bit-error rate is close to :math:`P_s/\log_2 M`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    J. W. Craig, "A new, simple and exact result for calculating the
    probability of error for two-dimensional signal constellations," *IEEE
    MILCOM* (1991) 571-575.

    Examples
    --------
    >>> bool(np.isclose(ser_mpsk_awgn(4.0, 2), ber_bpsk_awgn(4.0)))
    True
    """
    _nonnegative(esn0=esn0)
    _integer("M", M, 2)
    g = np.asarray(esn0, dtype=float)
    s2 = np.sin(np.pi / M) ** 2
    upper = (M - 1) * np.pi / M
    values = [integrate.quad(lambda th, x=x: np.exp(-x * s2 / max(np.sin(th) ** 2, 1e-300)), 0, upper)[0] / np.pi for x in g.ravel()]
    p = np.array(values).reshape(g.shape)
    return float(p) if p.ndim == 0 else p


def ser_mqam_awgn(esn0, M):
    r"""Exact symbol-error rate of coherent square M-QAM in AWGN.

    .. math::

        P_s = 1 - \left[1 - 2\left(1 - \frac{1}{\sqrt M}\right)Q\left(\sqrt{\frac{3E_s}{(M-1)N_0}}\right)\right]^2

    Parameters
    ----------
    esn0 : float or array_like
        :math:`E_s/N_0` [-] (linear), nonnegative.
    M : int
        Number of points, a power of four, at least 4.

    Returns
    -------
    float or ndarray
        Symbol-error probability [-]; each rail is an independent
        :math:`\sqrt M`-level amplitude decision.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> bool(np.isclose(ser_mqam_awgn(10.0, 4), ser_mpsk_awgn(10.0, 4)))
    True
    """
    _nonnegative(esn0=esn0)
    _integer("M", M, 4)
    k = int(M).bit_length() - 1
    if M & (M - 1) or k % 2:
        raise ValueError("M must be a power of four")
    rail = 2 * (1 - 1 / np.sqrt(M)) * q_function(np.sqrt(3 * np.asarray(esn0, dtype=float) / (M - 1)))
    return 1 - (1 - rail) ** 2
