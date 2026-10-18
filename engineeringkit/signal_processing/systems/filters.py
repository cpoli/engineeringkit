r"""Digital filters: Butterworth IIR design, Kaiser-window FIR design, filtering and frequency response.

Hurewicz (1947) described a sampled system by its z-transform
:math:`H(z) = B(z)/A(z)`; its frequency response is :math:`H` on the unit
circle, :math:`z = e^{j2\pi f/f_s}`. Butterworth's (1930) maximally flat
low-pass, carried into discrete time by the bilinear transform, has
:math:`|H|^2 = 1/(1 + (\tan(\pi f/f_s)/\tan(\pi f_c/f_s))^{2n})`. Kaiser
(1966) designed FIR filters by windowing the ideal impulse response, with
formulas for the length and window shape that meet a given attenuation.
Jackson (1970) showed why high-order recursive filters are built as
cascades of second-order sections: the roots of one high-degree polynomial
are too sensitive to rounding in its coefficients.
"""

from __future__ import annotations

import numpy as np
from scipy import signal

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive

__all__ = ["butterworth_lowpass", "butterworth_highpass", "fir_lowpass", "kaiser_order", "filter_signal", "frequency_response"]


def _butterworth(order, cutoff, fs, btype, output):
    _integer("order", order, 1)
    _positive(cutoff=cutoff, fs=fs)
    if cutoff >= fs / 2:
        raise ValueError("cutoff must lie below the Nyquist frequency fs/2")
    if output not in ("ba", "sos"):
        raise ValueError('output must be "ba" or "sos"')
    return signal.butter(order, cutoff, btype=btype, fs=fs, output=output)


def butterworth_lowpass(order, cutoff, fs, output="ba"):
    r"""Digital Butterworth low-pass filter by the bilinear transform.

    Parameters
    ----------
    order : int
        Filter order :math:`n`, at least 1.
    cutoff : float
        -3 dB frequency :math:`f_c` [Hz], between 0 and :math:`f_s/2`.
    fs : float
        Sampling rate [Hz], positive.
    output : {"ba", "sos"}, optional
        ``"ba"`` for numerator and denominator polynomials in
        :math:`z^{-1}`; ``"sos"`` for a cascade of second-order sections.

    Returns
    -------
    tuple of ndarray or ndarray
        ``(b, a)`` [-], or an ``(n_sections, 6)`` array of
        ``[b0, b1, b2, 1, a1, a2]`` rows [-].

    Notes
    -----
    :math:`|H(f)|^2 = 1/(1 + (\tan(\pi f/f_s)/\tan(\pi f_c/f_s))^{2n})`: the
    first :math:`2n-1` derivatives vanish at :math:`f = 0`, and the gain is
    :math:`1/\sqrt2` at :math:`f_c` for every order. Above about order 8,
    or for a cutoff far below :math:`f_s/2`, the ``"ba"`` polynomials lose
    precision and the filter can come out wrong or unstable; use ``"sos"``.
    See :doc:`/guides/signal_processing`.

    References
    ----------
    S. Butterworth, "On the theory of filter amplifiers," *Wireless
    Engineer* 7 (1930) 536-541; L. B. Jackson, "On the interaction of
    roundoff noise and dynamic range in digital filters," *Bell Syst. Tech.
    J.* 49 (1970) 159-184.

    Examples
    --------
    >>> b, a = butterworth_lowpass(1, 250.0, 1000.0)
    >>> np.round(b, 12).tolist()
    [0.5, 0.5]
    """
    return _butterworth(order, cutoff, fs, "lowpass", output)


def butterworth_highpass(order, cutoff, fs, output="ba"):
    r"""Digital Butterworth high-pass filter by the bilinear transform.

    Parameters
    ----------
    order : int
        Filter order :math:`n`, at least 1.
    cutoff : float
        -3 dB frequency :math:`f_c` [Hz], between 0 and :math:`f_s/2`.
    fs : float
        Sampling rate [Hz], positive.
    output : {"ba", "sos"}, optional
        Polynomials or second-order sections, as in
        :func:`butterworth_lowpass`.

    Returns
    -------
    tuple of ndarray or ndarray
        ``(b, a)`` or second-order sections [-].

    Notes
    -----
    :math:`|H(f)|^2 = 1/(1 + (\tan(\pi f_c/f_s)/\tan(\pi f/f_s))^{2n})`. See
    :func:`butterworth_lowpass` on precision, and
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> b, a = butterworth_highpass(1, 250.0, 1000.0)
    >>> b.tolist()
    [0.5, -0.5]
    """
    return _butterworth(order, cutoff, fs, "highpass", output)


def fir_lowpass(numtaps, cutoff, fs, window="hamming"):
    r"""Linear-phase FIR low-pass filter by the window method (SciPy's ``firwin``).

    Parameters
    ----------
    numtaps : int
        Number of coefficients :math:`N`, at least 1.
    cutoff : float
        Cutoff frequency [Hz], between 0 and :math:`f_s/2`; the gain there
        is about one half (-6 dB).
    fs : float
        Sampling rate [Hz], positive.
    window : str or tuple, optional
        Window, e.g. ``"hamming"`` or ``("kaiser", beta)``; see
        :func:`kaiser_order`.

    Returns
    -------
    ndarray
        Symmetric impulse response [-], scaled to unit gain at DC.

    Notes
    -----
    The ideal low-pass impulse response
    :math:`2(f_c/f_s)\operatorname{sinc}(2 f_c (n - (N-1)/2)/f_s)` is
    truncated to :math:`N` samples and tapered by the window. Symmetry gives
    a constant group delay of :math:`(N-1)/2` samples. See
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> round(float(fir_lowpass(31, 100.0, 1000.0).sum()), 12)
    1.0
    """
    _integer("numtaps", numtaps, 1)
    _positive(cutoff=cutoff, fs=fs)
    return signal.firwin(numtaps, cutoff, fs=fs, window=window)


def kaiser_order(attenuation_db, transition_width, fs):
    r"""Kaiser's formulas for the FIR length and window shape that meet a stopband attenuation.

    Parameters
    ----------
    attenuation_db : float
        Required stopband attenuation (and passband ripple)
        :math:`A = -20\log_{10}\delta` [dB], positive.
    transition_width : float
        Width of the transition band :math:`\Delta f` [Hz], positive and
        below :math:`f_s/2`.
    fs : float
        Sampling rate [Hz], positive.

    Returns
    -------
    tuple
        ``(numtaps, beta)``: number of coefficients (int) and Kaiser window
        parameter :math:`\beta` [-], for ``window=("kaiser", beta)``.

    Notes
    -----
    With :math:`\Delta\omega = 2\pi\Delta f/f_s`,
    :math:`N - 1 = \lceil (A - 7.95)/(2.285\,\Delta\omega) \rceil` and
    :math:`\beta = 0.1102(A - 8.7)` for :math:`A > 50`,
    :math:`0.5842(A-21)^{0.4} + 0.07886(A-21)` for :math:`21 \le A \le 50`,
    and 0 below. The formulas are empirical and usually land within about
    a decibel of :math:`A`. See :doc:`/guides/signal_processing`.

    References
    ----------
    J. F. Kaiser, "Digital filters," ch. 7 in F. F. Kuo and J. F. Kaiser
    (eds.), *System Analysis by Digital Computer*, Wiley (1966); J. F.
    Kaiser, "Nonrecursive digital filter design using the I0-sinh window
    function," *Proc. IEEE Int. Symp. Circuits Syst.* (1974) 20-23.

    Examples
    --------
    >>> numtaps, beta = kaiser_order(60.0, 50.0, 1000.0)
    >>> numtaps, round(beta, 4)
    (74, 5.6533)
    """
    _positive(attenuation_db=attenuation_db, transition_width=transition_width, fs=fs)
    if transition_width >= fs / 2:
        raise ValueError("transition_width must lie below fs/2")
    a = float(attenuation_db)
    if a > 50:
        beta = 0.1102 * (a - 8.7)
    elif a >= 21:
        beta = 0.5842 * (a - 21) ** 0.4 + 0.07886 * (a - 21)
    else:
        beta = 0.0
    numtaps = int(np.ceil((a - 7.95) / (2.285 * 2 * np.pi * transition_width / fs))) + 1
    return max(numtaps, 1), beta


def _sections(sos):
    sos = np.asarray(sos, dtype=float)
    if sos.ndim != 2 or sos.shape[1] != 6 or not sos.shape[0]:
        raise ValueError("second-order sections must be an (n_sections, 6) array; pass a=None only with sections")
    _finite(sos=sos)
    return sos


def filter_signal(b, a, x):
    r"""Apply a digital filter, given as polynomials or as second-order sections.

    Parameters
    ----------
    b : array_like
        Numerator coefficients in :math:`z^{-1}` [-], or, when ``a`` is
        ``None``, an ``(n_sections, 6)`` array of second-order sections.
    a : array_like or None
        Denominator coefficients [-], or ``None`` for sections.
    x : array_like
        Input samples [signal units].

    Returns
    -------
    ndarray
        Output :math:`y`, with :math:`\sum_k a_k y_{n-k} = \sum_k b_k x_{n-k}`
        and zero initial state [signal units].

    Notes
    -----
    The filter runs forward once, so its start-up transient and phase delay
    remain. Sections are applied one after another (``scipy.signal.sosfilt``);
    see :func:`butterworth_lowpass` and :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> filter_signal([1, 1], [1], [1, 0, 0]).tolist()
    [1.0, 1.0, 0.0]
    """
    if a is None:
        return signal.sosfilt(_sections(b), np.asarray(x, dtype=np.result_type(x, np.float64)))
    b, a, x = np.asarray(b), np.asarray(a), np.asarray(x)
    dtype = np.result_type(b, a, x, np.float64)
    return signal.lfilter(b.astype(dtype), a.astype(dtype), x.astype(dtype))


def frequency_response(b, a=1, worN=512, fs=2 * np.pi):
    r"""Frequency response :math:`H(e^{j2\pi f/f_s})` of a digital filter.

    Parameters
    ----------
    b : array_like
        Numerator coefficients [-], or second-order sections when ``a`` is
        ``None``.
    a : array_like or None, optional
        Denominator coefficients [-]; 1 for an FIR filter, ``None`` for
        sections.
    worN : int or array_like, optional
        Number of frequencies from 0 to :math:`f_s/2`, or the frequencies
        themselves [units of ``fs``].
    fs : float, optional
        Sampling rate [Hz]; the default :math:`2\pi` gives frequencies in
        rad/sample.

    Returns
    -------
    tuple of ndarray
        ``(frequency, H)``: frequencies [units of ``fs``] and complex gain [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> f, H = frequency_response([0.5, 0.5], 1, worN=[0.0, 250.0], fs=1000.0)
    >>> np.round(np.abs(H), 6).tolist()
    [1.0, 0.707107]
    """
    if a is None:
        return signal.sosfreqz(_sections(b), worN=worN, fs=fs)
    return signal.freqz(b, a, worN=worN, fs=fs)
