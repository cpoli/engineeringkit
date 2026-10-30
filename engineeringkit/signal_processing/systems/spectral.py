r"""Spectral estimation: the periodogram, windows, Welch averaging, the STFT, Goertzel's single bin and the cepstrum.

Schuster (1898) searched records for hidden periodicities with the
periodogram :math:`|X_k|^2/N`. Its variance does not fall as the record
grows; Welch (1967) averaged the periodograms of overlapping, windowed
segments to trade resolution for variance. Windows (Blackman and Tukey,
1958) taper each record to reduce leakage, and Gabor's (1946) short-time
transform slides a window along the signal to show how its spectrum
changes. Goertzel (1958) computes a single DFT bin with one real
multiplication per sample. Bogert, Healy and Tukey (1963) took the spectrum
of the log spectrum, the cepstrum, in which an echo becomes a single peak at
its delay.
"""

from __future__ import annotations

import numpy as np
from scipy import signal

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive
from engineeringkit.signal_processing.core.base import STFTResult

__all__ = ["power_spectral_density", "windowed_spectrum", "welch_psd", "stft", "goertzel", "real_cepstrum"]


def _real_signal(x, minimum=1):
    x = np.asarray(x)
    _finite(x=x)
    if x.ndim != 1 or np.iscomplexobj(x) or x.size < minimum:
        raise ValueError(f"x must be a real 1D signal with at least {minimum} samples")
    return x.astype(float)


def power_spectral_density(x, fs):
    r"""One-sided periodogram :math:`P_k = \frac{2|X_k|^2}{f_s N}` (SciPy's ``periodogram``).

    Parameters
    ----------
    x : array_like
        Samples [signal units].
    fs : float
        Sampling rate [Hz], positive.

    Returns
    -------
    tuple of ndarray
        ``(frequency, psd)``: bin frequencies [Hz] and power spectral density
        [signal units squared per Hz]; its sum times the bin spacing is the
        signal variance.

    Notes
    -----
    The periodogram's standard deviation equals its mean at every bin, so a
    longer record sharpens it without smoothing it; see :func:`welch_psd`.
    SciPy removes the mean before transforming. See
    :doc:`/guides/signal_processing`.

    References
    ----------
    A. Schuster, "On the investigation of hidden periodicities with
    application to a supposed 26 day period of meteorological phenomena,"
    *Terr. Magn.* 3 (1898) 13-41.
    """
    _positive(fs=fs)
    return signal.periodogram(x, fs=fs)


def windowed_spectrum(x, fs, window="hann", n_fft=None):
    r"""One-sided magnitude :math:`|\mathrm{FFT}(w\,x)|` of a windowed record.

    Parameters
    ----------
    x : array_like
        Real samples [signal units].
    fs : float
        Sampling rate [Hz], positive.
    window : str or tuple, optional
        Window name accepted by :func:`scipy.signal.get_window`, e.g.
        ``"boxcar"``, ``"hann"``, ``"hamming"``, ``"blackman"`` or
        ``("kaiser", beta)``.
    n_fft : int, optional
        Transform length, at least ``len(x)``; the windowed record is
        zero-padded to interpolate the spectrum. Defaults to ``len(x)``.

    Returns
    -------
    tuple of ndarray
        ``(frequency, magnitude)``: bin frequencies [Hz] and unnormalized
        magnitude [signal units].

    Notes
    -----
    A tone between bins leaks into its neighbours through the window's
    sidelobes. The window does not change the bin spacing :math:`f_s/N`;
    zero-padding only samples the same spectrum more finely. See
    :doc:`/guides/signal_processing`.

    References
    ----------
    R. B. Blackman and J. W. Tukey, "The measurement of power spectra from
    the point of view of communications engineering," *Bell Syst. Tech. J.*
    37 (1958) 185-282 and 485-569.

    Examples
    --------
    >>> f, magnitude = windowed_spectrum([1.0, 1.0, 1.0, 1.0], 4.0, window="boxcar")
    >>> magnitude.tolist()
    [4.0, 0.0, 0.0]
    """
    _positive(fs=fs)
    x = np.asarray(x)
    n = len(x) if n_fft is None else n_fft
    _integer("n_fft", n, len(x))
    w = signal.get_window(window, len(x))
    return np.fft.rfftfreq(n, 1 / fs), np.abs(np.fft.rfft(x * w, n))


def welch_psd(x, fs, segment_length=256, overlap=None, window="hann"):
    r"""Welch's averaged, modified periodogram, written out segment by segment.

    Parameters
    ----------
    x : array_like
        Real samples [signal units].
    fs : float
        Sampling rate [Hz], positive.
    segment_length : int, optional
        Samples per segment :math:`L`, between 2 and ``len(x)``.
    overlap : int, optional
        Samples shared by consecutive segments, from 0 to :math:`L-1`.
        Defaults to :math:`L/2` (rounded down).
    window : str or tuple, optional
        Window applied to each segment, as in :func:`windowed_spectrum`.

    Returns
    -------
    tuple of ndarray
        ``(frequency, psd)``: one-sided bin frequencies [Hz] and power
        spectral density [signal units squared per Hz].

    Notes
    -----
    Each segment's periodogram is :math:`|\mathrm{FFT}(w\,x)|^2/(f_s\sum w^2)`,
    doubled for every bin but DC and Nyquist, and the segments are averaged.
    Averaging :math:`K` nearly independent segments divides the variance by
    about :math:`K` and widens the resolution to about :math:`f_s/L`. No
    detrending is applied, so the result matches
    ``scipy.signal.welch(..., detrend=False)``. See
    :doc:`/guides/signal_processing`.

    References
    ----------
    P. D. Welch, "The use of fast Fourier transform for the estimation of
    power spectra: a method based on time averaging over short, modified
    periodograms," *IEEE Trans. Audio Electroacoust.* 15 (1967) 70-73.

    Examples
    --------
    >>> f, psd = welch_psd(np.ones(8), 8.0, segment_length=4, window="boxcar")
    >>> psd.tolist()
    [0.5, 0.0, 0.0]
    """
    _positive(fs=fs)
    x = _real_signal(x, 2)
    _integer("segment_length", segment_length, 2)
    if segment_length > x.size:
        raise ValueError("segment_length must not exceed the signal length")
    overlap = segment_length // 2 if overlap is None else overlap
    _integer("overlap", overlap, 0)
    if overlap >= segment_length:
        raise ValueError("overlap must be smaller than segment_length")
    w = signal.get_window(window, segment_length)
    step = segment_length - overlap
    starts = range(0, x.size - segment_length + 1, step)
    spectra = [np.abs(np.fft.rfft(w * x[s : s + segment_length])) ** 2 for s in starts]
    psd = np.mean(spectra, axis=0) / (fs * np.sum(w**2))
    psd[1 : (segment_length + 1) // 2] *= 2
    return np.fft.rfftfreq(segment_length, 1 / fs), psd


def stft(x, fs, window_length, hop=None, window="hann"):
    r"""Short-time Fourier transform :math:`X_m(f_k) = \sum_n x_{n+mH}\,w_n\,e^{-j2\pi kn/L}`.

    Parameters
    ----------
    x : array_like
        Real samples [signal units].
    fs : float
        Sampling rate [Hz], positive.
    window_length : int
        Frame length :math:`L` [samples], between 2 and ``len(x)``.
    hop : int, optional
        Frame step :math:`H` [samples], at least 1. Defaults to :math:`L/4`
        (rounded down, at least 1).
    window : str or tuple, optional
        Analysis window, as in :func:`windowed_spectrum`; Gabor's choice is
        ``("gaussian", std)``.

    Returns
    -------
    STFTResult
        Frame centre times [s], one-sided bin frequencies [Hz] and the
        complex coefficients [signal units].

    Notes
    -----
    A window of duration :math:`\Delta t` resolves frequencies only to about
    :math:`1/\Delta t`. For RMS widths of :math:`|w|^2` in time and of
    :math:`|W|^2` in frequency, :math:`\sigma_t\sigma_f \ge 1/(4\pi)`, with
    equality only for a Gaussian. Frames are not padded, so the first and
    last centres lie :math:`L/2` samples inside the record. See
    :doc:`/guides/signal_processing`.

    References
    ----------
    D. Gabor, "Theory of communication," *J. IEE* 93, part III (1946)
    429-457.

    Examples
    --------
    >>> result = stft(np.ones(16), 16.0, 8, hop=4, window="boxcar")
    >>> result.summary()["frames"], float(abs(result.values[0, 0]))
    (3, 8.0)
    """
    _positive(fs=fs)
    x = _real_signal(x, 2)
    _integer("window_length", window_length, 2)
    if window_length > x.size:
        raise ValueError("window_length must not exceed the signal length")
    hop = max(1, window_length // 4) if hop is None else hop
    _integer("hop", hop, 1)
    w = signal.get_window(window, window_length)
    starts = np.arange(0, x.size - window_length + 1, hop)
    frames = x[starts[:, None] + np.arange(window_length)] * w
    return STFTResult((starts + window_length / 2) / fs, np.fft.rfftfreq(window_length, 1 / fs), np.fft.rfft(frames, axis=1).T)


def goertzel(x, frequency, fs):
    r"""Goertzel's algorithm: the DFT of ``x`` at one frequency from a second-order recursion.

    Parameters
    ----------
    x : array_like
        Real samples [signal units].
    frequency : float or array_like
        Frequency or frequencies to evaluate [Hz]; need not fall on a bin.
    fs : float
        Sampling rate [Hz], positive.

    Returns
    -------
    complex or ndarray of complex
        :math:`X(f) = \sum_n x_n e^{-j2\pi fn/f_s}` [signal units]; at
        :math:`f = kf_s/N` this is DFT bin :math:`k`.

    Notes
    -----
    The recursion :math:`s_n = x_n + 2\cos\omega\,s_{n-1} - s_{n-2}`, with
    :math:`\omega = 2\pi f/f_s`, costs one real multiplication per sample.
    One extra step with zero input and
    :math:`X = e^{-j\omega N}(s_N - e^{-j\omega}s_{N-1})` finish it. It
    beats an FFT when fewer than about :math:`\log_2 N` bins are wanted, as
    in a touch-tone (DTMF) decoder. See :doc:`/guides/signal_processing`.

    References
    ----------
    G. Goertzel, "An algorithm for the evaluation of finite trigonometric
    series," *Amer. Math. Monthly* 65 (1958) 34-35.

    Examples
    --------
    >>> X = goertzel([1.0, 0.0, -1.0, 0.0], 1.0, 4.0)
    >>> round(float(X.real), 12), round(float(X.imag), 12)
    (2.0, 0.0)
    """
    _positive(fs=fs)
    _finite(frequency=frequency)
    x = _real_signal(x)
    omega = 2 * np.pi * np.asarray(frequency, dtype=float) / fs
    coefficient = 2 * np.cos(omega)
    s1 = np.zeros_like(omega)
    s2 = np.zeros_like(omega)
    for sample in np.append(x, 0.0):
        s1, s2 = sample + coefficient * s1 - s2, s1
    return (np.exp(-1j * omega * x.size) * (s1 - np.exp(-1j * omega) * s2))[()]


def real_cepstrum(x, n_fft=None):
    r"""Real cepstrum :math:`c_n = \mathrm{IDFT}(\ln|\mathrm{DFT}(x)|)`, indexed by quefrency :math:`n/f_s`.

    Parameters
    ----------
    x : array_like
        Real 1D samples [signal units].
    n_fft : int, optional
        Transform length, at least ``len(x)``; zero-padding reduces the
        wrap-around of the cepstrum. Defaults to ``len(x)``.

    Returns
    -------
    ndarray
        Real, even cepstrum of length ``n_fft`` [-]; sample :math:`n` is
        quefrency :math:`n/f_s` [s].

    Notes
    -----
    An echo :math:`x_n + \alpha x_{n-d}` multiplies the spectrum by
    :math:`1 + \alpha e^{-j\omega d}`, a ripple in :math:`\ln|X|` with
    "frequency" :math:`d`; its log expands as
    :math:`\sum_m (-1)^{m+1}(\alpha^m/m)\cos(m\omega d)`, so for
    :math:`|\alpha| < 1` the cepstrum adds :math:`\alpha/2` at :math:`n = d`
    (and smaller peaks at multiples of :math:`d`) to that of :math:`x`.
    Convolution becomes addition, which is how speech is split into pitch
    and vocal tract. A zero in the spectrum has no logarithm and raises
    :class:`ValueError`. See :doc:`/guides/signal_processing`.

    References
    ----------
    B. P. Bogert, M. J. R. Healy and J. W. Tukey, "The quefrency alanysis
    of time series for echoes: cepstrum, pseudo-autocovariance,
    cross-cepstrum and saphe cracking," in M. Rosenblatt (ed.),
    *Proceedings of the Symposium on Time Series Analysis*, Wiley (1963)
    209-243.

    Examples
    --------
    >>> x = np.zeros(64)
    >>> x[0], x[10] = 1.0, 0.5
    >>> round(float(real_cepstrum(x)[10]), 6)
    0.25
    """
    x = _real_signal(x)
    n = x.size if n_fft is None else n_fft
    _integer("n_fft", n, x.size)
    magnitude = np.abs(np.fft.rfft(x, n))
    if np.any(magnitude <= np.finfo(float).tiny):
        raise ValueError("the spectrum of x has a zero, so its logarithm is undefined")
    return np.fft.irfft(np.log(magnitude), n)
