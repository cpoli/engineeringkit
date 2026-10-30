r"""Sampling, quantization and sample-rate conversion: aliasing, sinc reconstruction, the uniform ADC, decimation and interpolation.

Shannon (1949) proved that a signal with no energy above :math:`B` is fixed
by samples taken at :math:`f_s > 2B`, and rebuilt by a sum of sinc pulses.
A tone above :math:`f_s/2` is indistinguishable from its alias in
:math:`[0, f_s/2]`. Bennett (1948) showed that a fine uniform quantizer of
step :math:`q` adds noise of power :math:`q^2/12`, which gives the
:math:`6.02N + 1.76` dB rule for an :math:`N`-bit converter. Schafer and
Rabiner (1973) treated changing the rate of a sampled signal as digital
filtering: insert zeros and low-pass filter to interpolate, low-pass filter
and keep every :math:`M`-th sample to decimate.
"""

from __future__ import annotations

import numpy as np
from scipy import signal

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive

__all__ = [
    "sampling_period",
    "nyquist_frequency",
    "alias_frequency",
    "sinc_reconstruct",
    "quantization_step",
    "quantize",
    "downsample",
    "upsample",
    "decimate",
    "interpolate",
]


def sampling_period(fs):
    r"""Sampling period :math:`T = 1/f_s`.

    Parameters
    ----------
    fs : float or array_like
        Sampling rate [Hz], positive.

    Returns
    -------
    float or ndarray
        Sampling period [s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> float(sampling_period(10))
    0.1
    """
    _positive(fs=fs)
    return 1 / np.asarray(fs, dtype=float)[()]


def nyquist_frequency(fs):
    r"""Nyquist frequency :math:`f_s/2`, the highest frequency the samples represent unambiguously.

    Parameters
    ----------
    fs : float or array_like
        Sampling rate [Hz], positive.

    Returns
    -------
    float or ndarray
        Nyquist frequency [Hz].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> float(nyquist_frequency(100))
    50.0
    """
    _positive(fs=fs)
    return np.asarray(fs, dtype=float)[()] / 2


def alias_frequency(f, fs):
    r"""Apparent frequency in :math:`[0, f_s/2]` of a tone at ``f`` sampled at ``fs``.

    Parameters
    ----------
    f : float or array_like
        Tone frequency [Hz]; negative frequencies fold like positive ones.
    fs : float
        Sampling rate [Hz], positive.

    Returns
    -------
    float or ndarray
        Alias frequency :math:`|((f + f_s/2) \bmod f_s) - f_s/2|` [Hz].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> float(alias_frequency(70, 100))
    30.0
    """
    _finite(f=f)
    _positive(fs=fs)
    return np.abs(((np.asarray(f, dtype=float)[()] + fs / 2) % fs) - fs / 2)


def sinc_reconstruct(samples, fs, t):
    r"""Whittaker-Shannon interpolation :math:`x(t) = \sum_n x_n \operatorname{sinc}(f_s t - n)`.

    Parameters
    ----------
    samples : array_like
        Samples :math:`x_n` taken at :math:`t_n = n/f_s` [signal units].
    fs : float
        Sampling rate [Hz], positive.
    t : array_like
        Times at which to evaluate the reconstruction [s].

    Returns
    -------
    ndarray
        Reconstructed signal at ``t`` [signal units].

    Notes
    -----
    The sum is truncated to the samples supplied, so it is exact at the
    sample instants but only approximate between them near the ends of the
    record. See :doc:`/guides/signal_processing`.

    References
    ----------
    C. E. Shannon, "Communication in the presence of noise," *Proc. IRE* 37
    (1949) 10-21.

    Examples
    --------
    >>> sinc_reconstruct([0.0, 1.0, 0.0], 10.0, [0.1]).tolist()
    [1.0]
    """
    _positive(fs=fs)
    samples = np.asarray(samples)
    tn = np.arange(len(samples)) / fs
    return np.sinc(fs * (np.asarray(t)[:, None] - tn)) @ samples


def quantization_step(vmin, vmax, bits):
    r"""Step :math:`q = (v_{max} - v_{min})/2^N` of a uniform mid-rise ADC.

    Parameters
    ----------
    vmin, vmax : float
        Input range [V or any amplitude unit], ``vmax > vmin``.
    bits : int
        Resolution :math:`N`, between 1 and 52.

    Returns
    -------
    float
        Quantization step [same unit as ``vmin``].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> quantization_step(-1.0, 1.0, 3)
    0.25
    """
    if not isinstance(bits, (int, np.integer)) or isinstance(bits, bool) or not 1 <= bits <= 52:
        raise ValueError("bits must be an integer between 1 and 52")
    if not (np.all(np.isfinite(vmin)) and np.all(np.isfinite(vmax))) or np.any(np.asarray(vmax) <= vmin):
        raise ValueError("finite vmax must exceed vmin")
    return (vmax - vmin) / (2**bits)


def quantize(x, vmin, vmax, bits):
    r"""Ideal uniform mid-rise quantization, saturating outside the ADC range.

    Parameters
    ----------
    x : float or array_like
        Input samples [V or any amplitude unit].
    vmin, vmax : float
        Input range [same unit], ``vmax > vmin``.
    bits : int
        Resolution :math:`N`, between 1 and 52.

    Returns
    -------
    float or ndarray
        One of the :math:`2^N` reconstruction levels at the bin centres
        :math:`v_{min} + (k + 1/2)q` [same unit].

    Notes
    -----
    In range, the error never exceeds :math:`q/2`; for a busy signal it is
    close to uniform on :math:`[-q/2, q/2]`, with power :math:`q^2/12`.
    See :doc:`/guides/signal_processing`.

    References
    ----------
    W. R. Bennett, "Spectra of quantized signals," *Bell Syst. Tech. J.* 27
    (1948) 446-472.

    Examples
    --------
    >>> quantize([-0.9, 0.1], -1.0, 1.0, 3).tolist()
    [-0.875, 0.125]
    """
    q = quantization_step(vmin, vmax, bits)
    x = np.asarray(x, float)
    if not np.all(np.isfinite(x)):
        raise ValueError("samples must be finite")
    codes = np.clip(np.floor((x - vmin) / q), 0, 2**bits - 1)
    return vmin + (codes + 0.5) * q


def _rate_input(x, factor):
    _integer("factor", factor, 1)
    x = np.asarray(x)
    _finite(x=x)
    if x.ndim != 1:
        raise ValueError("x must be a 1D signal")
    return x.astype(np.result_type(x, np.float64))


def _odd_taps(numtaps, factor):
    numtaps = 20 * factor + 1 if numtaps is None else numtaps
    _integer("numtaps", numtaps, 1)
    if numtaps % 2 == 0:
        raise ValueError("numtaps must be odd, so the filter delay is a whole number of samples")
    return numtaps


def downsample(x, factor, phase=0):
    r"""Keep every ``factor``-th sample, :math:`y_m = x_{mM + p}`, with no filtering.

    Parameters
    ----------
    x : array_like
        1D samples [signal units].
    factor : int
        Downsampling factor :math:`M`, at least 1.
    phase : int, optional
        Offset :math:`p` of the first kept sample, from 0 to :math:`M-1`.

    Returns
    -------
    ndarray
        Samples at rate :math:`f_s/M` [signal units].

    Notes
    -----
    The spectrum of the result is the sum of :math:`M` copies of the input
    spectrum stretched by :math:`M` and shifted by multiples of
    :math:`2\pi/M`, so any content above :math:`f_s/(2M)` aliases. Use
    :func:`decimate` to filter first. See :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> downsample([0.0, 1.0, 2.0, 3.0, 4.0], 2).tolist()
    [0.0, 2.0, 4.0]
    """
    x = _rate_input(x, factor)
    _integer("phase", phase, 0)
    if phase >= factor:
        raise ValueError("phase must be smaller than factor")
    return x[phase::factor]


def upsample(x, factor):
    r"""Insert :math:`L - 1` zeros after every sample, :math:`y_{nL} = x_n`.

    Parameters
    ----------
    x : array_like
        1D samples [signal units].
    factor : int
        Upsampling factor :math:`L`, at least 1.

    Returns
    -------
    ndarray
        ``len(x) * L`` samples at rate :math:`L f_s` [signal units].

    Notes
    -----
    :math:`Y(e^{j\omega}) = X(e^{j\omega L})`: the spectrum is compressed by
    :math:`L` and :math:`L - 1` images of it appear up to the new Nyquist
    frequency. :func:`interpolate` removes them. See
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> upsample([1.0, 2.0], 3).tolist()
    [1.0, 0.0, 0.0, 2.0, 0.0, 0.0]
    """
    x = _rate_input(x, factor)
    y = np.zeros(x.size * factor, dtype=x.dtype)
    y[::factor] = x
    return y


def decimate(x, factor, numtaps=None):
    r"""Reduce the sampling rate by an integer factor: anti-alias low-pass FIR filter, then downsample.

    Parameters
    ----------
    x : array_like
        1D samples [signal units].
    factor : int
        Decimation factor :math:`M`, at least 1.
    numtaps : int, optional
        Odd length of the Hamming-window FIR filter, cutoff
        :math:`f_s/(2M)`. Defaults to :math:`20M + 1`.

    Returns
    -------
    ndarray
        ``ceil(len(x)/M)`` samples at rate :math:`f_s/M` [signal units],
        aligned with the input: the filter's :math:`(N-1)/2`-sample delay is
        removed.

    Notes
    -----
    The filter removes what would alias before :func:`downsample` discards
    samples. Its gain is about -6 dB at the new Nyquist frequency and its
    stopband about 53 dB down, so keep wanted content below roughly
    :math:`0.8 f_s/(2M)`. The record is treated as zero outside its ends,
    so the first and last :math:`(N-1)/2` input samples carry a transient.
    For large :math:`M`, decimating in stages is cheaper. See
    :doc:`/guides/signal_processing`.

    References
    ----------
    R. E. Crochiere and L. R. Rabiner, *Multirate Digital Signal
    Processing*, Prentice-Hall (1983), ch. 2.

    Examples
    --------
    >>> np.round(decimate(np.ones(200), 4)[20:22], 6).tolist()
    [1.0, 1.0]
    """
    x = _rate_input(x, factor)
    numtaps = _odd_taps(numtaps, factor)
    h = signal.firwin(numtaps, 1 / factor) if factor > 1 else np.ones(1)
    delay = (numtaps - 1) // 2 if factor > 1 else 0
    return np.convolve(x, h)[delay : delay + x.size : factor]


def interpolate(x, factor, numtaps=None):
    r"""Raise the sampling rate by an integer factor: upsample, then low-pass filter with gain :math:`L`.

    Parameters
    ----------
    x : array_like
        1D samples [signal units].
    factor : int
        Interpolation factor :math:`L`, at least 1.
    numtaps : int, optional
        Odd length of the Hamming-window FIR filter, cutoff
        :math:`f_s/2` of the input. Defaults to :math:`20L + 1`.

    Returns
    -------
    ndarray
        ``len(x) * L`` samples at rate :math:`L f_s` [signal units],
        aligned with the input.

    Notes
    -----
    The filter is the windowed ideal interpolator
    :math:`h_n = w_n \operatorname{sinc}(n/L)`, which is 1 at :math:`n = 0`
    and 0 at every other multiple of :math:`L`, so the original samples pass
    through unchanged and only the new ones are computed. It is Shannon's
    sinc reconstruction, truncated to :math:`N` taps; the images left by
    :func:`upsample` are attenuated by about 53 dB. See
    :doc:`/guides/signal_processing`.

    References
    ----------
    R. W. Schafer and L. R. Rabiner, "A digital signal processing approach
    to interpolation," *Proc. IEEE* 61 (1973) 692-702,
    doi:10.1109/PROC.1973.9150.

    Examples
    --------
    >>> y = interpolate([0.0, 1.0, 0.0, -1.0], 2)
    >>> np.round(y[::2], 12).tolist()
    [0.0, 1.0, 0.0, -1.0]
    """
    x = _rate_input(x, factor)
    numtaps = _odd_taps(numtaps, factor)
    u = upsample(x, factor)
    if factor == 1:
        return u
    h = factor * signal.firwin(numtaps, 1 / factor, scale=False)
    delay = (numtaps - 1) // 2
    return np.convolve(u, h)[delay : delay + u.size]
