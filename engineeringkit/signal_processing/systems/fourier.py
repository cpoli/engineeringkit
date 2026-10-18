r"""Fourier analysis: the DFT and its inverse, a hand-written radix-2 FFT, and Fourier partial sums.

The DFT :math:`X_k = \sum_n x_n e^{-j2\pi kn/N}` costs :math:`N^2` complex
multiplications when evaluated directly. Cooley and Tukey (1965) split it
into the DFTs of the even and odd samples, recursively, for
:math:`(N/2)\log_2 N`. A truncated Fourier series of a discontinuous
waveform overshoots by about 9% of the jump however many terms are kept,
the phenomenon Gibbs (1899) explained.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive

__all__ = ["dft", "idft", "fft_frequencies", "radix2_fft", "fourier_partial_sum"]


def dft(x):
    r"""Discrete Fourier transform :math:`X_k = \sum_{n=0}^{N-1} x_n e^{-j2\pi kn/N}` (NumPy's FFT).

    Parameters
    ----------
    x : array_like
        Samples, real or complex [signal units].

    Returns
    -------
    ndarray of complex
        Unnormalized coefficients :math:`X_k` [signal units].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> dft([1.0, 1.0, 1.0, 1.0]).real.tolist()
    [4.0, 0.0, 0.0, 0.0]
    """
    return np.fft.fft(x)


def idft(X):
    r"""Inverse DFT :math:`x_n = \frac{1}{N}\sum_k X_k e^{j2\pi kn/N}` (NumPy's inverse FFT).

    Parameters
    ----------
    X : array_like
        DFT coefficients [signal units].

    Returns
    -------
    ndarray of complex
        Samples [signal units].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> idft([4.0, 0.0, 0.0, 0.0]).real.tolist()
    [1.0, 1.0, 1.0, 1.0]
    """
    return np.fft.ifft(X)


def fft_frequencies(n, fs):
    r"""Frequency of each DFT bin, :math:`k f_s/N`, with the upper half negative.

    Parameters
    ----------
    n : int
        Transform length :math:`N`.
    fs : float
        Sampling rate [Hz], positive.

    Returns
    -------
    ndarray
        Bin frequencies [Hz] in NumPy's FFT order.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> fft_frequencies(4, 4.0).tolist()
    [0.0, 1.0, -2.0, -1.0]
    """
    _positive(fs=fs)
    return np.fft.fftfreq(n, 1 / fs)


def _radix2(x):
    n = x.size
    if n == 1:
        return x, 0
    even, m_even = _radix2(x[0::2])
    odd, m_odd = _radix2(x[1::2])
    twiddled = np.exp(-2j * np.pi * np.arange(n // 2) / n) * odd
    return np.concatenate([even + twiddled, even - twiddled]), m_even + m_odd + n // 2


def radix2_fft(x, return_multiplications=False):
    r"""Recursive decimation-in-time radix-2 FFT, written out to show the butterflies.

    Parameters
    ----------
    x : array_like
        Samples, real or complex [signal units]; the length must be a power
        of two.
    return_multiplications : bool, optional
        Also return the number of complex twiddle multiplications.

    Returns
    -------
    ndarray of complex or tuple
        DFT coefficients, equal to :func:`dft` up to rounding, and, if
        requested, the multiplication count :math:`(N/2)\log_2 N`.

    Notes
    -----
    Each stage combines the DFTs :math:`E_k` and :math:`O_k` of the even and
    odd samples as :math:`X_k = E_k + W^k O_k` and
    :math:`X_{k+N/2} = E_k - W^k O_k`, with :math:`W = e^{-j2\pi/N}`. The
    count includes the trivial twiddle :math:`W^0 = 1`. See
    :doc:`/guides/signal_processing`.

    References
    ----------
    J. W. Cooley and J. W. Tukey, "An algorithm for the machine calculation
    of complex Fourier series," *Math. Comp.* 19 (1965) 297-301.

    Examples
    --------
    >>> X, count = radix2_fft([1.0, 0.0, 0.0, 0.0], return_multiplications=True)
    >>> X.real.tolist(), count
    ([1.0, 1.0, 1.0, 1.0], 4)
    """
    x = np.asarray(x, dtype=complex)
    _finite(x=x)
    if x.ndim != 1 or x.size & (x.size - 1):
        raise ValueError("x must be a 1D array whose length is a power of two")
    X, multiplications = _radix2(x)
    return (X, multiplications) if return_multiplications else X


def fourier_partial_sum(t, period, a, b, a0=0.0):
    r"""Truncated Fourier series :math:`a_0 + \sum_{k=1}^{K} [a_k\cos(k\omega t) + b_k\sin(k\omega t)]`.

    Parameters
    ----------
    t : float or array_like
        Times [s].
    period : float
        Period :math:`2\pi/\omega` [s], positive.
    a, b : array_like
        Cosine and sine coefficients :math:`a_1..a_K` and :math:`b_1..b_K`
        [signal units], of equal length.
    a0 : float, optional
        Mean value [signal units].

    Returns
    -------
    float or ndarray
        Partial sum at ``t`` [signal units].

    Notes
    -----
    At a jump discontinuity the partial sums overshoot by
    :math:`\frac{1}{\pi}\operatorname{Si}(\pi) - \frac12 \approx 0.0895`
    of the jump, however large :math:`K` is; the overshoot only moves
    closer to the jump. See :doc:`/guides/signal_processing`.

    References
    ----------
    J. W. Gibbs, "Fourier's series," *Nature* 59 (1899) 606.

    Examples
    --------
    >>> float(fourier_partial_sum(0.0, 1.0, [1.0, 0.5], [0.0, 0.0], a0=0.25))
    1.75
    """
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    _finite(t=t, a0=a0)
    _positive(period=period)
    if a.ndim != 1 or a.shape != b.shape or not a.size:
        raise ValueError("a and b must be nonempty 1D arrays of equal length")
    _finite(a=a, b=b)
    phase = 2 * np.pi / period * np.multiply.outer(np.asarray(t, dtype=float), np.arange(1, a.size + 1))
    return a0 + np.cos(phase) @ a + np.sin(phase) @ b
