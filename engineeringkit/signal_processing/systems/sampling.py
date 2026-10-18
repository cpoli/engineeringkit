r"""Sampling and quantization: the Nyquist rate, aliasing, sinc reconstruction and the uniform ADC.

Shannon (1949) proved that a signal with no energy above :math:`B` is fixed
by samples taken at :math:`f_s > 2B`, and rebuilt by a sum of sinc pulses.
A tone above :math:`f_s/2` is indistinguishable from its alias in
:math:`[0, f_s/2]`. Bennett (1948) showed that a fine uniform quantizer of
step :math:`q` adds noise of power :math:`q^2/12`, which gives the
:math:`6.02N + 1.76` dB rule for an :math:`N`-bit converter.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive

__all__ = ["sampling_period", "nyquist_frequency", "alias_frequency", "sinc_reconstruct", "quantization_step", "quantize"]


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
    if not np.all(np.isfinite([vmin, vmax])) or vmax <= vmin:
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
