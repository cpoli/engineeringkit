r"""Linear time-invariant operations in the time domain: convolution, correlation and the moving average.

A linear time-invariant system is fixed by its impulse response
:math:`h`, and its output is the convolution
:math:`y_n = \sum_k h_k x_{n-k}`. Correlation, convolution with a
time-reversed (conjugated) sequence, measures how alike two signals are at
each lag.
"""

from __future__ import annotations

import numpy as np
from scipy import signal

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer

__all__ = ["convolution", "correlation", "moving_average"]


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
