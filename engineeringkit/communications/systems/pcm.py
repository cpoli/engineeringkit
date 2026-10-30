r"""Waveform digitization: pulse-code modulation, mu-law companding and delta modulation.

Reeves's pulse-code modulation samples a waveform, rounds each sample to one
of :math:`2^N` levels and sends the level number as :math:`N` binary digits,
so that repeaters can regenerate the pulses without accumulating noise. Each
extra bit buys about 6 dB of signal-to-quantization-noise ratio. Delta
modulation sends one bit per sample instead, saying only whether the signal
has risen or fallen by a fixed step.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive

__all__ = [
    "pcm_encode",
    "pcm_decode",
    "pcm_sqnr_db",
    "mu_law_compress",
    "mu_law_expand",
    "delta_modulate",
    "delta_demodulate",
    "delta_slope_limit",
]


def _resolution(n_bits):
    if isinstance(n_bits, (bool, np.bool_)) or not isinstance(n_bits, (int, np.integer)) or not 1 <= n_bits <= 52:
        raise ValueError("n_bits must be an integer between 1 and 52")
    return int(n_bits)


def _binary(bits, name, block=1):
    b = np.asarray(bits)
    if b.ndim != 1 or not b.size or b.size % block or not np.all((b == 0) | (b == 1)):
        raise ValueError(f"{name} must be a nonempty binary vector whose length is a multiple of {block}")
    return b.astype(np.int64)


def pcm_encode(x, n_bits, full_scale=1.0):
    r"""Pulse-code modulation: quantize samples uniformly and write each level number in binary.

    The quantizer is mid-rise over :math:`[-V, V]` with step
    :math:`q = 2V/2^N`; samples outside the range saturate at the end
    levels.

    Parameters
    ----------
    x : array_like
        Samples [V or any amplitude unit].
    n_bits : int
        Bits per sample :math:`N`, between 1 and 52.
    full_scale : float, optional
        Peak amplitude :math:`V` [same unit], positive.

    Returns
    -------
    ndarray of int
        :math:`N` bits per sample, most significant first, the level number
        counted upward from :math:`-V`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    A. H. Reeves, French Patent 852,183 (filed 1938, published 1939); U.S.
    Patent 2,272,070 (1942).

    Examples
    --------
    >>> pcm_encode([-0.9, 0.1, 0.99], 3).tolist()
    [0, 0, 0, 1, 0, 0, 1, 1, 1]
    """
    N = _resolution(n_bits)
    _positive(full_scale=full_scale)
    samples = np.asarray(x, dtype=float).ravel()
    if not samples.size:
        raise ValueError("x must be nonempty")
    _finite(x=samples)
    q = 2 * full_scale / 2**N
    codes = np.clip(np.floor((samples + full_scale) / q), 0, 2**N - 1).astype(np.int64)
    return ((codes[:, None] >> np.arange(N - 1, -1, -1)) & 1).reshape(-1)


def pcm_decode(bits, n_bits, full_scale=1.0):
    r"""Reconstruct PCM samples at the centres of their quantization intervals.

    Parameters
    ----------
    bits : array_like
        Bit stream from :func:`pcm_encode`, :math:`N` bits per sample.
    n_bits : int
        Bits per sample :math:`N`, between 1 and 52.
    full_scale : float, optional
        Peak amplitude :math:`V`, positive.

    Returns
    -------
    ndarray
        Reconstructed samples :math:`-V + (k + 1/2)q` [same unit as
        ``full_scale``]; in range, each is within :math:`q/2` of the input.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> pcm_decode([0, 0, 0, 1, 0, 0], 3).tolist()
    [-0.875, 0.125]
    """
    N = _resolution(n_bits)
    _positive(full_scale=full_scale)
    words = _binary(bits, "bits", N).reshape(-1, N)
    codes = words @ (1 << np.arange(N - 1, -1, -1))
    q = 2 * full_scale / 2**N
    return -full_scale + (codes + 0.5) * q


def pcm_sqnr_db(n_bits):
    r"""Signal-to-quantization-noise ratio of a full-scale sine through an ``N``-bit PCM quantizer.

    With noise power :math:`q^2/12` and signal power :math:`V^2/2`, the
    ratio is :math:`\tfrac32 4^N`, or :math:`6.02N + 1.76` dB.

    Parameters
    ----------
    n_bits : int
        Bits per sample :math:`N`, between 1 and 52.

    Returns
    -------
    float
        SQNR [dB]. It is the large-:math:`N` value: a simulated full-scale
        sine falls short by about :math:`2.3 \cdot 2^{-N/2}` dB (0.14 dB at
        8 bits), because the sine dwells near its peaks, where the error in
        the top interval is not uniform.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    B. M. Oliver, J. R. Pierce and C. E. Shannon, "The philosophy of PCM,"
    *Proc. IRE* 36 (1948) 1324-1331.

    Examples
    --------
    >>> round(pcm_sqnr_db(8), 2)
    49.93
    """
    N = _resolution(n_bits)
    return float(10 * np.log10(1.5) + 20 * N * np.log10(2))


def mu_law_compress(x, mu=255.0, full_scale=1.0):
    r"""Mu-law compressor :math:`y = V\operatorname{sgn}(x)\ln(1 + \mu|x|/V)/\ln(1 + \mu)`.

    Compressing before a uniform quantizer gives small signals finer steps
    than loud ones, so the SQNR stays nearly constant over a wide range of
    speech levels. North American and Japanese telephony use
    :math:`\mu = 255` with 8-bit PCM.

    Parameters
    ----------
    x : float or array_like
        Samples with :math:`|x| \le V` [V or any amplitude unit].
    mu : float, optional
        Compression parameter :math:`\mu` [-], positive.
    full_scale : float, optional
        Peak amplitude :math:`V` [same unit], positive.

    Returns
    -------
    float or ndarray
        Compressed samples in :math:`[-V, V]`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    B. Smith, "Instantaneous companding of quantized signals," *Bell Syst.
    Tech. J.* 36 (1957) 653-709.

    Examples
    --------
    >>> float(mu_law_compress(1.0)), round(float(mu_law_compress(0.01)), 4)
    (1.0, 0.2285)
    """
    _positive(mu=mu, full_scale=full_scale)
    _interval("x", np.asarray(x, dtype=float) / full_scale, -1, 1)
    x = np.asarray(x, dtype=float)
    return full_scale * np.sign(x) * np.log1p(mu * np.abs(x) / full_scale) / np.log1p(mu)


def mu_law_expand(y, mu=255.0, full_scale=1.0):
    r"""Inverse of :func:`mu_law_compress`, :math:`x = V\operatorname{sgn}(y)[(1 + \mu)^{|y|/V} - 1]/\mu`.

    Parameters
    ----------
    y : float or array_like
        Compressed samples with :math:`|y| \le V`.
    mu : float, optional
        Compression parameter :math:`\mu` [-], positive.
    full_scale : float, optional
        Peak amplitude :math:`V`, positive.

    Returns
    -------
    float or ndarray
        Expanded samples [same unit as ``full_scale``].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> round(float(mu_law_expand(mu_law_compress(0.3))), 12)
    0.3
    """
    _positive(mu=mu, full_scale=full_scale)
    _interval("y", np.asarray(y, dtype=float) / full_scale, -1, 1)
    y = np.asarray(y, dtype=float)
    return full_scale * np.sign(y) * np.expm1(np.abs(y) / full_scale * np.log1p(mu)) / mu


def delta_modulate(x, step, initial=0.0):
    r"""Linear delta modulation: one bit per sample saying whether to step the staircase up or down.

    The encoder tracks its own staircase approximation :math:`\hat x`: bit
    :math:`b_n = 1` if :math:`x_n \ge \hat x_{n-1}`, and then
    :math:`\hat x_n = \hat x_{n-1} \pm \Delta`.

    Parameters
    ----------
    x : array_like
        Samples [V or any amplitude unit], taken well above the Nyquist rate.
    step : float
        Step size :math:`\Delta` [same unit], positive.
    initial : float, optional
        Starting value of the staircase [same unit].

    Returns
    -------
    ndarray of int
        One bit per sample; :func:`delta_demodulate` rebuilds the staircase.

    Notes
    -----
    A small step causes *slope overload* where the signal changes faster
    than :func:`delta_slope_limit`; a large one causes *granular noise* where
    it is flat. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    F. de Jager, "Delta modulation, a method of P.C.M. transmission using
    the 1-unit code," *Philips Res. Rep.* 7 (1952) 442-466.

    Examples
    --------
    >>> delta_modulate([0.0, 0.3, 0.6, 0.6, 0.2], 0.25).tolist()
    [1, 1, 1, 0, 0]
    """
    samples = np.asarray(x, dtype=float)
    if samples.ndim != 1 or not samples.size:
        raise ValueError("x must be a nonempty 1-D array")
    _finite(x=samples, initial=initial)
    _positive(step=step)
    bits = np.empty(samples.size, dtype=np.int64)
    level = float(initial)
    for n, sample in enumerate(samples):
        bits[n] = sample >= level
        level += step if bits[n] else -step
    return bits


def delta_demodulate(bits, step, initial=0.0):
    r"""Staircase reconstruction of a delta-modulated signal, :math:`\hat x_n = \hat x_0 + \Delta\sum(2b_k - 1)`.

    Parameters
    ----------
    bits : array_like
        Bits from :func:`delta_modulate`.
    step : float
        Step size :math:`\Delta`, positive.
    initial : float, optional
        Starting value of the staircase.

    Returns
    -------
    ndarray
        Staircase samples [same unit as ``step``]; a low-pass filter smooths
        them into the recovered waveform.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> delta_demodulate([1, 1, 1, 0, 0], 0.25).tolist()
    [0.25, 0.5, 0.75, 0.5, 0.25]
    """
    b = _binary(bits, "bits")
    _finite(initial=initial)
    _positive(step=step)
    return initial + step * np.cumsum(2 * b - 1)


def delta_slope_limit(step, fs):
    r"""Fastest slope a delta modulator can follow without overload, :math:`\Delta f_s`.

    Parameters
    ----------
    step : float or array_like
        Step size :math:`\Delta` [V or any amplitude unit], positive.
    fs : float or array_like
        Sampling rate [Hz], positive.

    Returns
    -------
    float or ndarray
        Slope limit [amplitude unit per second]. A sine of amplitude
        :math:`A` and frequency :math:`f` avoids overload when
        :math:`2\pi f A \le \Delta f_s`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(delta_slope_limit(0.01, 64e3))
    640.0
    """
    _positive(step=step, fs=fs)
    return np.asarray(step, dtype=float) * fs
