r"""Linear time-invariant operations: linear and circular convolution, overlap-add, correlation and the moving average.

A linear time-invariant system is fixed by its impulse response
:math:`h`, and its output is the convolution
:math:`y_n = \sum_k h_k x_{n-k}`. Correlation, convolution with a
time-reversed (conjugated) sequence, measures how alike two signals are at
each lag. Multiplying two DFTs gives the *circular* convolution, in which
the tail wraps around to the start; zero-padding both sequences to at least
:math:`N_x + N_h - 1` samples makes it equal the linear one. Stockham
(1966) used this to filter long records with the FFT, block by block, by
overlap-add.
"""

from __future__ import annotations

import numpy as np
from scipy import signal

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer

__all__ = ["convolution", "circular_convolution", "overlap_add", "correlation", "moving_average"]


def convolution(x, h, mode="full"):
    r"""Linear convolution :math:`y_n = \sum_k h_k x_{n-k}` (NumPy's ``convolve``).

    Parameters
    ----------
    x : array_like
        Input samples [signal units].
    h : array_like
        Impulse response [-].
    mode : {"full", "same", "valid"}, optional
        Output length, as in :func:`numpy.convolve`.

    Returns
    -------
    ndarray
        Output samples [signal units].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> convolution([1, 2], [1, 1]).tolist()
    [1, 3, 2]
    """
    return np.convolve(x, h, mode)


def _sequence(name, value):
    a = np.asarray(value)
    _finite(**{name: a})
    if a.ndim != 1:
        raise ValueError(f"{name} must be a 1D sequence")
    return a.astype(np.result_type(a, np.float64))


def circular_convolution(x, h, n=None):
    r"""Circular (periodic) convolution :math:`y_m = \sum_{k=0}^{N-1} h_k\,x_{(m-k) \bmod N}`, written out.

    Parameters
    ----------
    x : array_like
        1D input samples [signal units].
    h : array_like
        1D impulse response [-].
    n : int, optional
        Period :math:`N`, at least ``max(len(x), len(h))``; both sequences
        are zero-padded to it. Defaults to that maximum.

    Returns
    -------
    ndarray
        The :math:`N` output samples [signal units].

    Notes
    -----
    This is what multiplying two :math:`N`-point DFTs computes:
    :math:`\mathrm{DFT}(y) = \mathrm{DFT}(x)\,\mathrm{DFT}(h)`. The last
    :math:`N_h - 1` samples of the linear convolution wrap around and add
    onto the first ones (time aliasing) unless
    :math:`N \ge N_x + N_h - 1`, in which case the two agree. See
    :func:`overlap_add` and :doc:`/guides/signal_processing`.

    References
    ----------
    A. V. Oppenheim and R. W. Schafer, *Discrete-Time Signal Processing*,
    3rd ed., Pearson (2010), sec. 8.6-8.7.

    Examples
    --------
    >>> circular_convolution([1.0, 2.0, 3.0], [1.0, 1.0]).tolist()
    [4.0, 3.0, 5.0]
    >>> circular_convolution([1.0, 2.0, 3.0], [1.0, 1.0], n=4).tolist()
    [1.0, 3.0, 5.0, 3.0]
    """
    x, h = _sequence("x", x), _sequence("h", h)
    n = max(x.size, h.size) if n is None else n
    _integer("n", n, max(x.size, h.size))
    x, h = np.pad(x, (0, n - x.size)), np.pad(h, (0, n - h.size))
    return x[(np.arange(n)[:, None] - np.arange(n)) % n] @ h


def overlap_add(x, h, block_length=None):
    r"""Linear convolution of a long record with an FIR filter by FFT blocks, overlap-add.

    Parameters
    ----------
    x : array_like
        1D input samples [signal units].
    h : array_like
        1D impulse response of length :math:`M` [-].
    block_length : int, optional
        Input samples per block :math:`L`, at least 1. Defaults to the
        largest :math:`L` for which :math:`L + M - 1` fits in the smallest
        power of two of at least :math:`4M`.

    Returns
    -------
    ndarray
        The :math:`N_x + M - 1` samples of ``convolution(x, h)`` [signal
        units], real when both inputs are real.

    Notes
    -----
    The record is cut into blocks of :math:`L` samples. Each block and
    :math:`h` are zero-padded to :math:`N_{FFT} \ge L + M - 1`, so that the
    circular convolution of the block equals its linear one, multiplied in
    the frequency domain, and the :math:`L + M - 1` output samples are
    added into place, overlapping the next block's by :math:`M - 1`. The
    cost per output sample grows as :math:`\log N_{FFT}` rather than
    :math:`M`. See :func:`circular_convolution` and
    :doc:`/guides/signal_processing`.

    References
    ----------
    T. G. Stockham Jr., "High-speed convolution and correlation," *AFIPS
    Conf. Proc.* 28 (1966 Spring Joint Computer Conf.) 229-233,
    doi:10.1145/1464182.1464209.

    Examples
    --------
    >>> np.round(overlap_add([1.0, 2.0, 3.0, 4.0], [1.0, 1.0], block_length=2), 12).tolist()
    [1.0, 3.0, 5.0, 7.0, 4.0]
    """
    x, h = _sequence("x", x), _sequence("h", h)
    m = h.size
    if block_length is None:
        n_fft = 1 << int(np.ceil(np.log2(4 * m)))
        block_length = n_fft - m + 1
    _integer("block_length", block_length, 1)
    n_fft = 1 << int(np.ceil(np.log2(block_length + m - 1)))
    H = np.fft.fft(h, n_fft)
    y = np.zeros(x.size + n_fft, dtype=complex)
    for start in range(0, x.size, block_length):
        y[start : start + n_fft] += np.fft.ifft(np.fft.fft(x[start : start + block_length], n_fft) * H)
    y = y[: x.size + m - 1]
    return y if np.iscomplexobj(x) or np.iscomplexobj(h) else y.real


def correlation(x, y, mode="full"):
    r"""Cross-correlation :math:`r_k = \sum_n x_{n+k}\,\overline{y_n}` (SciPy's ``correlate``).

    Parameters
    ----------
    x, y : array_like
        Samples [signal units].
    mode : {"full", "same", "valid"}, optional
        Output length, as in :func:`scipy.signal.correlate`; ``"full"``
        returns lags :math:`-(N_y-1)` to :math:`N_x-1`.

    Returns
    -------
    ndarray
        Correlation at each lag [signal units squared].

    Notes
    -----
    ``correlation(x, x)`` is the unnormalized autocorrelation; divide by
    :math:`N` for the biased estimate whose DFT is the periodogram. See
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> correlation([1, 2], [1, 2]).tolist()
    [2, 5, 2]
    """
    return signal.correlate(x, y, mode=mode)


def moving_average(x, n):
    r"""Moving average over ``n`` samples, :math:`y_m = \frac{1}{n}\sum_{k=0}^{n-1} x_{m+k}`.

    Parameters
    ----------
    x : array_like
        1D input samples [signal units].
    n : int
        Window length, between 1 and ``len(x)``.

    Returns
    -------
    ndarray
        The ``len(x) - n + 1`` fully overlapped averages [signal units].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> moving_average([1.0, 2.0, 3.0], 2).tolist()
    [1.5, 2.5]
    """
    _integer("n", n, 1)
    x = np.asarray(x)
    _finite(x=x)
    if x.ndim != 1 or n > len(x):
        raise ValueError("window n must not exceed the 1D signal length")
    return np.convolve(x, np.ones(n) / n, "valid")
