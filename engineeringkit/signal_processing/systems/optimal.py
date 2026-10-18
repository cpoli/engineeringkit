r"""Optimal and adaptive filters: the FIR Wiener filter and Widrow and Hoff's LMS algorithm.

Wiener (1949) found the filter that minimizes the mean squared error
between its output and a desired signal, given only their correlations.
For an FIR filter :math:`w` with input :math:`x` and desired output
:math:`d` this is the solution of the normal equations :math:`Rw = p`,
with :math:`R` the Toeplitz autocorrelation matrix of :math:`x` and
:math:`p` its cross-correlation with :math:`d`. Widrow and Hoff's (1960)
least-mean-squares rule :math:`w \leftarrow w + \mu e_n u_n` reaches the
same solution by stochastic gradient descent, without knowing :math:`R`
or :math:`p`, and keeps tracking it if the statistics drift.
"""

from __future__ import annotations

import numpy as np
from scipy import linalg

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive
from engineeringkit.signal_processing.core.base import LMSResult

__all__ = ["wiener_filter", "lms_filter"]


def _pair(x, d, n_taps):
    x, d = np.asarray(x), np.asarray(d)
    _finite(x=x, d=d)
    if x.ndim != 1 or x.shape != d.shape or np.iscomplexobj(x) or np.iscomplexobj(d):
        raise ValueError("x and d must be real 1D signals of equal length")
    _integer("n_taps", n_taps, 1)
    if n_taps > x.size:
        raise ValueError("n_taps must not exceed the signal length")
    return x.astype(float), d.astype(float)


def _tap_inputs(x, n_taps):
    """Rows ``u_n = [x_n, x_{n-1}, ..., x_{n-L+1}]`` with zeros before the record starts."""
    padded = np.concatenate([np.zeros(n_taps - 1), x])
    return padded[np.arange(x.size)[:, None] + np.arange(n_taps - 1, -1, -1)]


def wiener_filter(x, d, n_taps):
    r"""FIR Wiener filter: the weights that minimize :math:`E[(d_n - \sum_k w_k x_{n-k})^2]`.

    Parameters
    ----------
    x : array_like
        Real input samples [signal units].
    d : array_like
        Desired output, same length [signal units].
    n_taps : int
        Filter length :math:`L`, between 1 and ``len(x)``.

    Returns
    -------
    ndarray
        Weights :math:`w_0..w_{L-1}` [-], to be applied as
        ``filter_signal(w, 1, x)``.

    Notes
    -----
    Solves :math:`Rw = p` with the biased sample estimates
    :math:`r_k = \frac1N\sum_n x_n x_{n-k}` and
    :math:`p_k = \frac1N\sum_n d_n x_{n-k}`, by Levinson recursion
    (``scipy.linalg.solve_toeplitz``). The minimum mean squared error is
    :math:`E[d^2] - p^T w`. Signals are assumed stationary. See
    :doc:`/guides/signal_processing`.

    References
    ----------
    N. Wiener, *Extrapolation, Interpolation, and Smoothing of Stationary
    Time Series*, MIT Press and Wiley (1949); N. Levinson, "The Wiener RMS
    (root mean square) error criterion in filter design and prediction,"
    *J. Math. Phys.* 25 (1947) 261-278.

    Examples
    --------
    >>> x = np.array([1.0, -1.0, 2.0, 0.5, -1.5, 1.0])
    >>> np.round(wiener_filter(x, 3 * x, 2), 12).tolist()
    [3.0, 0.0]
    """
    x, d = _pair(x, d, n_taps)
    n = x.size
    r = np.array([x[k:] @ x[: n - k] for k in range(n_taps)]) / n
    p = np.array([d[k:] @ x[: n - k] for k in range(n_taps)]) / n
    if r[0] == 0:
        raise ValueError("x must not be identically zero")
    return linalg.solve_toeplitz(r, p)


def lms_filter(x, d, n_taps, mu, initial_weights=None):
    r"""Widrow-Hoff least-mean-squares adaptive FIR filter.

    Parameters
    ----------
    x : array_like
        Real input samples [signal units].
    d : array_like
        Desired output, same length [signal units].
    n_taps : int
        Filter length :math:`L`, between 1 and ``len(x)``.
    mu : float
        Step size [1/signal units squared], positive.
    initial_weights : array_like, optional
        Starting weights, length ``n_taps`` [-]; zeros by default.

    Returns
    -------
    LMSResult
        Output, error, final weights and the weight history.

    Notes
    -----
    For each sample, with :math:`u_n = [x_n, \ldots, x_{n-L+1}]`:
    :math:`y_n = w^T u_n`, :math:`e_n = d_n - y_n` and
    :math:`w \leftarrow w + \mu e_n u_n`. The mean weights converge to the
    Wiener solution when :math:`0 < \mu < 2/\lambda_{max}(R)`; in practice
    :math:`\mu \ll 2/(L\,E[x^2])` keeps the excess error small. See
    :doc:`/guides/signal_processing`.

    References
    ----------
    B. Widrow and M. E. Hoff, "Adaptive switching circuits," *IRE WESCON
    Convention Record*, part 4 (1960) 96-104.

    Examples
    --------
    >>> rng = np.random.default_rng(0)
    >>> x = rng.normal(size=2000)
    >>> result = lms_filter(x, 0.5 * x, 1, 0.05)
    >>> round(float(result.weights[0]), 6)
    0.5
    """
    x, d = _pair(x, d, n_taps)
    _positive(mu=mu)
    w = np.zeros(n_taps) if initial_weights is None else np.array(initial_weights, dtype=float)
    _finite(initial_weights=w)
    if w.shape != (n_taps,):
        raise ValueError("initial_weights must have length n_taps")
    u = _tap_inputs(x, n_taps)
    y, e, history = np.empty(x.size), np.empty(x.size), np.empty((x.size, n_taps))
    for n in range(x.size):
        history[n] = w
        y[n] = w @ u[n]
        e[n] = d[n] - y[n]
        w = w + mu * e[n] * u[n]
    return LMSResult(y, e, w, history)
