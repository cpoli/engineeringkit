r"""Weibull parameters from life data: median ranks, probability plots and maximum likelihood.

Leonard Johnson showed that the :math:`i`-th of :math:`n` failures is best
plotted at its *median rank*, the median of the beta distribution of the
:math:`i`-th order statistic, and Benard's :math:`(i - 0.3)/(n + 0.4)`
approximates it within a fraction of a percent. Johnson's adjusted ranks
carry censored units into the plot. On Weibull probability axes the data
fall on a line whose slope is the shape :math:`\beta`. Cohen wrote down the
maximum-likelihood equations for complete and censored samples, which use
every censored unit's survival time instead of only its effect on the
ranks.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize, stats

from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.reliability.core.base import WeibullFitResult

__all__ = ["median_ranks", "adjusted_ranks", "weibull_rank_regression", "weibull_mle"]

_RANK_METHODS = ("exact", "benard")


def median_ranks(n, order=None, method="exact"):
    r"""Median-rank plotting positions of ordered failures in a sample of ``n``.

    The :math:`i`-th smallest of :math:`n` lives has a population failure
    fraction distributed as :math:`\mathrm{Beta}(i, n - i + 1)`; its median
    is the exact median rank (Johnson). Benard's approximation is
    :math:`(i - 0.3)/(n + 0.4)`.

    Parameters
    ----------
    n : int
        Number of units on test, at least one.
    order : array_like, optional
        Order numbers :math:`i` in [1, n]; fractional values from
        :func:`adjusted_ranks` are allowed. Default: ``1, 2, ..., n``.
    method : {"exact", "benard"}, optional
        Beta-distribution median or Benard's approximation.

    Returns
    -------
    ndarray
        Estimated failure probability at each order number [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    L. G. Johnson, "The median ranks of sample values in their population
    with an application to certain fatigue studies," *Industrial
    Mathematics* 2 (1951) 1-9; A. Benard and E. C. Bos-Levenbach, "Het
    uitzetten van waarnemingen op waarschijnlijkheids-papier," *Statistica
    Neerlandica* 7 (1953) 163-173.

    Examples
    --------
    >>> np.round(median_ranks(5), 4)
    array([0.1294, 0.3138, 0.5   , 0.6862, 0.8706])
    >>> np.round(median_ranks(5, method="benard"), 4)
    array([0.1296, 0.3148, 0.5   , 0.6852, 0.8704])
    """
    _integer("n", n, 1)
    if method not in _RANK_METHODS:
        raise ValueError(f"method must be one of {_RANK_METHODS}")
    i = np.arange(1.0, n + 1) if order is None else _vector("order", order)
    _interval("order", i, 1, n)
    if method == "benard":
        return (i - 0.3) / (n + 0.4)
    return stats.beta.ppf(0.5, i, n - i + 1)


def _sorted_life_data(times, failed):
    """Validated times sorted increasing (failures before censored units at ties) and failure flags."""
    t = _vector("times", times)
    _positive(times=t)
    if failed is None:
        d = np.ones(t.size, dtype=bool)
    else:
        flags = np.asarray(failed)
        if flags.shape != t.shape or not np.all((flags == 0) | (flags == 1)):
            raise ValueError("failed must be a 0/1 vector matching times")
        d = flags.astype(bool)
    order = np.lexsort((~d, t))
    return t[order], d[order]


def adjusted_ranks(times, failed=None):
    r"""Johnson's adjusted order numbers of the failures in a censored sample.

    Units are sorted by time. At the unit in position :math:`j` with reverse
    rank :math:`r = n - j + 1`, a failure gets the order number
    :math:`O_j = (r\,O_{prev} + n + 1)/(r + 1)`, where :math:`O_{prev}` is
    that of the previous failure (zero at first). A censored unit takes no
    order number but makes every later increment larger, since it might
    have failed anywhere after its removal. Without censoring the order
    numbers are ``1, 2, ..., n``.

    Parameters
    ----------
    times : array_like
        Failure or censoring time of each unit [s], positive.
    failed : array_like of {0, 1} or bool, optional
        ``1`` where the unit failed, ``0`` where it was censored. Default:
        every unit failed. At a tie, failures are ranked first.

    Returns
    -------
    tuple of ndarray
        ``(failure_times, order)``: the failure times in increasing order
        [s] and their adjusted order numbers [-], for :func:`median_ranks`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    L. G. Johnson, *The Statistical Treatment of Fatigue Experiments*
    (Elsevier, Amsterdam, 1964), ch. 3.

    Examples
    --------
    >>> failure_times, order = adjusted_ranks([10, 20, 30, 40], [1, 0, 1, 1])
    >>> np.round(order, 3)
    array([1.   , 2.333, 3.667])
    """
    t, d = _sorted_life_data(times, failed)
    n = t.size
    order, previous = [], 0.0
    for j in np.flatnonzero(d):
        reverse = n - j
        previous = (reverse * previous + n + 1) / (reverse + 1)
        order.append(previous)
    return t[d], np.array(order)


def _failure_check(t, d, minimum):
    if np.count_nonzero(d) < minimum or (minimum > 1 and np.ptp(t[d]) == 0):
        raise ValueError(f"need at least {minimum} failures at different times")


def weibull_rank_regression(times, failed=None, regress="y", ranks="exact"):
    r"""Fit a two-parameter Weibull distribution by least squares on a probability plot.

    Each failure time :math:`t_i` is plotted at its median rank
    :math:`F_i` (with Johnson's adjusted ranks for censored data) on the
    axes :math:`x = \ln t`, :math:`y = \ln(-\ln(1 - F))`, where
    :math:`y = \beta x - \beta \ln\eta`. Regressing :math:`y` on :math:`x`
    is the textbook method; regressing :math:`x` on :math:`y` treats the
    times as the uncertain quantity.

    Parameters
    ----------
    times : array_like
        Failure or censoring time of each unit [s], positive.
    failed : array_like of {0, 1} or bool, optional
        ``1`` where the unit failed, ``0`` where it was censored. Default:
        every unit failed.
    regress : {"y", "x"}, optional
        Dependent variable of the least-squares line.
    ranks : {"exact", "benard"}, optional
        Median-rank method, as in :func:`median_ranks`.

    Returns
    -------
    WeibullFitResult
        Scale, shape, the plotted points and ``plot()`` for the probability
        plot.

    Raises
    ------
    ValueError
        With fewer than two failures at different times.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    L. G. Johnson, "The median ranks of sample values in their population
    with an application to certain fatigue studies," *Industrial
    Mathematics* 2 (1951) 1-9; C. E. Ebeling, *An Introduction to
    Reliability and Maintainability Engineering* (McGraw-Hill, 1997), §15.

    Examples
    --------
    >>> fit = weibull_rank_regression([16, 34, 53, 75, 93, 120])
    >>> round(fit.beta, 3), round(fit.eta, 1)
    (1.43, 76.3)
    """
    if regress not in ("y", "x"):
        raise ValueError('regress must be "y" or "x"')
    t, d = _sorted_life_data(times, failed)
    _failure_check(t, d, 2)
    failure_times, order = adjusted_ranks(t, d)
    f = median_ranks(t.size, order, ranks)
    x, y = np.log(failure_times), np.log(-np.log1p(-f))
    if regress == "y":
        beta, intercept = np.polyfit(x, y, 1)
        eta = np.exp(-intercept / beta)
    else:
        slope, log_eta = np.polyfit(y, x, 1)
        beta, eta = 1 / slope, np.exp(log_eta)
    return WeibullFitResult(float(eta), float(beta), f"rank regression on {regress}", failure_times, f, int(t.size))


def weibull_mle(times, failed=None, ranks="exact"):
    r"""Maximum-likelihood Weibull scale and shape from complete or censored life data.

    With :math:`r` failures among :math:`n` units and every unit's time
    :math:`t_i` (failure or censoring), the shape solves Cohen's equation

    .. math::

        \frac{\sum_{i=1}^n t_i^\beta \ln t_i}{\sum_{i=1}^n t_i^\beta}
        - \frac{1}{\beta} - \frac{1}{r}\sum_{\text{failures}} \ln t_i = 0,

    whose left side increases with :math:`\beta`, and then
    :math:`\eta = (\sum_i t_i^\beta / r)^{1/\beta}`.

    Parameters
    ----------
    times : array_like
        Failure or censoring time of each unit [s], positive.
    failed : array_like of {0, 1} or bool, optional
        ``1`` where the unit failed, ``0`` where it was censored. Default:
        every unit failed.
    ranks : {"exact", "benard"}, optional
        Median-rank method for the plotted points of the result; the fit
        does not use them.

    Returns
    -------
    WeibullFitResult
        Scale, shape and the points for ``plot()``.

    Raises
    ------
    ValueError
        With no failures, or when every failure is at the longest time on
        test (the likelihood then grows without bound in :math:`\beta`).

    Notes
    -----
    Times are scaled by the longest one before solving, which leaves the
    equation unchanged and keeps :math:`t^\beta` from overflowing. Inputs are
    raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    A. C. Cohen, "Maximum likelihood estimation in the Weibull distribution
    based on complete and on censored samples," *Technometrics* 7 (1965)
    579-588.

    Examples
    --------
    >>> fit = weibull_mle([16, 34, 53, 75, 93, 120])
    >>> round(fit.beta, 3), round(fit.eta, 1)
    (1.933, 73.5)
    """
    t, d = _sorted_life_data(times, failed)
    _failure_check(t, d, 1)
    if np.all(t[d] == t.max()):
        raise ValueError("every failure is at the longest time on test; the shape estimate is unbounded")
    scale = t.max()
    logs = np.log(t / scale)
    mean_failure_log = logs[d].mean()

    def equation(beta):
        w = np.exp(beta * logs)
        return float(w @ logs / w.sum() - 1 / beta - mean_failure_log)

    low, high = 1.0, 1.0
    while equation(low) > 0:
        low /= 2
    while equation(high) < 0:
        high *= 2
    beta = optimize.brentq(equation, low, high, xtol=1e-14, rtol=1e-12)
    eta = scale * (np.exp(beta * logs).sum() / np.count_nonzero(d)) ** (1 / beta)
    failure_times, order = adjusted_ranks(t, d)
    return WeibullFitResult(float(eta), float(beta), "maximum likelihood", failure_times, median_ranks(t.size, order, ranks), int(t.size))
