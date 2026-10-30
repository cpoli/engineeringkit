r"""Autoregressive models: the Levinson-Durbin recursion, Yule-Walker estimation and the AR spectrum.

Yule (1927) explained the sunspot cycle not as a hidden sine but as a
disturbed oscillator, :math:`x_n = -a_1 x_{n-1} - a_2 x_{n-2} + e_n`,
driven by random shocks, and fitted :math:`a_1, a_2` from the record's
correlations. Walker (1931) wrote the equations for any order. Their
matrix is Toeplitz, and Levinson (1947) and Durbin (1960) solved it by a
recursion on the order that costs :math:`p^2` operations instead of
:math:`p^3`, producing along the way the reflection coefficients, whose
magnitudes below one certify a stable model. The fitted model gives a
smooth spectrum :math:`\sigma^2/|A(e^{j\omega})|^2` with sharp peaks from
a short record.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit.signal_processing.core.base import ARModel

__all__ = ["levinson_durbin", "yule_walker"]


def levinson_durbin(r, order):
    r"""Levinson-Durbin recursion: AR coefficients of every order up to ``order`` from autocorrelations.

    Parameters
    ----------
    r : array_like
        Real autocorrelation :math:`r_0, r_1, \ldots` [signal units
        squared], with :math:`r_0 > 0` and at least ``order + 1`` values.
    order : int
        Model order :math:`p`, at least 1.

    Returns
    -------
    ARModel
        Coefficients :math:`[1, a_1..a_p]`, noise variance :math:`E_p`,
        reflection coefficients and prediction errors.

    Notes
    -----
    Starting from :math:`E_0 = r_0`, order :math:`m` adds the reflection
    coefficient :math:`k_m = -(r_m + \sum_{i=1}^{m-1} a_i r_{m-i})/E_{m-1}`,
    updates :math:`a_i \leftarrow a_i + k_m a_{m-i}` with :math:`a_m = k_m`,
    and shrinks the error to :math:`E_m = (1 - k_m^2)E_{m-1}`. The result
    solves the Yule-Walker equations :math:`\sum_k a_k r_{|i-k|} = -r_i`,
    :math:`i = 1..p`. A :math:`|k_m| \ge 1` means the sequence is not a
    valid (positive definite) autocorrelation. See
    :doc:`/guides/signal_processing`.

    References
    ----------
    N. Levinson, "The Wiener RMS (root mean square) error criterion in
    filter design and prediction," *J. Math. Phys.* 25 (1947) 261-278;
    J. Durbin, "The fitting of time-series models," *Rev. Int. Stat. Inst.*
    28 (1960) 233-244, doi:10.2307/1401322.

    Examples
    --------
    >>> model = levinson_durbin([1.0, 0.5, 0.25], 2)
    >>> (np.round(model.coefficients, 12) + 0.0).tolist(), round(model.noise_variance, 12)
    ([1.0, -0.5, 0.0], 0.75)
    """
    r = np.asarray(r)
    _finite(r=r)
    _integer("order", order, 1)
    if r.ndim != 1 or np.iscomplexobj(r) or r.size < order + 1 or r[0] <= 0:
        raise ValueError("r must be a real 1D autocorrelation with r[0] > 0 and at least order + 1 lags")
    r = r.astype(float)
    a = np.array([1.0])
    errors, reflections = [r[0]], []
    for m in range(1, order + 1):
        k = -(a @ r[m:0:-1]) / errors[-1]
        if abs(k) >= 1:
            raise ValueError(f"r is not positive definite: |k_{m}| = {abs(k):.3g} >= 1")
        a = np.append(a, 0.0) + k * np.append(a, 0.0)[::-1]
        reflections.append(k)
        errors.append(errors[-1] * (1 - k * k))
    return ARModel(a, float(errors[-1]), np.array(reflections), np.array(errors))


def yule_walker(x, order, demean=True):
    r"""Yule-Walker estimate of an AR model from a record, solved by :func:`levinson_durbin`.

    Parameters
    ----------
    x : array_like
        Real 1D samples [signal units], longer than ``order``.
    order : int
        Model order :math:`p`, at least 1.
    demean : bool, optional
        Subtract the record's mean first.

    Returns
    -------
    ARModel
        The fitted model; ``model.psd(fs)`` is the AR spectrum.

    Notes
    -----
    Uses the biased autocorrelation :math:`r_k = \frac1N\sum_n x_n x_{n+k}`,
    which is always positive definite, so the model is always stable. Its
    spectrum is smooth and can resolve two close peaks from a record too
    short for the periodogram, but a wrong order gives spurious or missing
    peaks. See :doc:`/guides/signal_processing`.

    References
    ----------
    G. U. Yule, "On a method of investigating periodicities in disturbed
    series, with special reference to Wolfer's sunspot numbers," *Phil.
    Trans. R. Soc. A* 226 (1927) 267-298, doi:10.1098/rsta.1927.0007;
    G. Walker, "On periodicity in series of related terms," *Proc. R. Soc.
    A* 131 (1931) 518-532.

    Examples
    --------
    >>> rng = np.random.default_rng(0)
    >>> e = rng.normal(size=20000)
    >>> x = np.zeros_like(e)
    >>> for n in range(1, e.size):
    ...     x[n] = 0.9 * x[n - 1] + e[n]
    >>> round(float(yule_walker(x, 1).coefficients[1]), 2)
    -0.9
    """
    x = np.asarray(x)
    _finite(x=x)
    _integer("order", order, 1)
    if x.ndim != 1 or np.iscomplexobj(x) or x.size <= order:
        raise ValueError("x must be a real 1D record longer than order")
    x = x.astype(float)
    if demean:
        x = x - x.mean()
    if not np.any(x):
        raise ValueError("x must not be constant")
    n = x.size
    r = np.array([x[: n - k] @ x[k:] for k in range(order + 1)]) / n
    return levinson_durbin(r, order)
