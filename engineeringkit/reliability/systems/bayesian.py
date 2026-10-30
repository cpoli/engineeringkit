r"""Bayesian updating of a constant failure rate with a conjugate gamma prior.

Field data on a new part are scarce, but handbooks and similar parts say
roughly what its failure rate should be. A gamma prior on :math:`\lambda`
encodes that knowledge, and because the gamma is conjugate to Poisson
failure counts (Raiffa and Schlaifer, 1961) the update is two additions:
:math:`r` failures in time :math:`T` turn a :math:`\mathrm{Gamma}(a, b)`
prior into a :math:`\mathrm{Gamma}(a + r, b + T)` posterior. The prior acts
like :math:`a` failures seen in :math:`b` seconds of earlier testing.
"""

from __future__ import annotations

from scipy import special

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.reliability.core.base import FailureRatePosterior

__all__ = ["bayesian_failure_rate", "gamma_prior_from_bounds"]


def bayesian_failure_rate(prior_shape, prior_rate, failures, total_time):
    r"""Posterior gamma distribution of a failure rate after a life test.

    Parameters
    ----------
    prior_shape : float
        Prior gamma shape :math:`a` [-], positive.
    prior_rate : float
        Prior gamma rate :math:`b` [s], positive; the prior mean rate is
        :math:`a/b`.
    failures : int
        Failures observed :math:`r`, nonnegative.
    total_time : float
        Total time on test summed over all units [s], nonnegative.

    Returns
    -------
    FailureRatePosterior
        The posterior :math:`\mathrm{Gamma}(a + r, b + T)`, with its mean,
        credible bounds and ``plot()`` of prior and posterior.

    Notes
    -----
    ``prior_shape``, ``prior_rate`` and ``total_time`` are scalar-only,
    since the result describes one posterior distribution; update several
    parts in a loop. The posterior mean :math:`(a + r)/(b + T)` is a weighted average of the
    prior mean and the classical estimate :math:`r/T`, with the data gaining
    weight as :math:`T` grows. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/reliability`.

    References
    ----------
    H. Raiffa and R. Schlaifer, *Applied Statistical Decision Theory*
    (Harvard University, Boston, 1961), ch. 3 and 7; H. F. Martz and R. A.
    Waller, *Bayesian Reliability Analysis* (Wiley, New York, 1982).

    Examples
    --------
    >>> posterior = bayesian_failure_rate(2.0, 1e6, 3, 4e6)
    >>> posterior.mean
    1e-06
    """
    _positive(prior_shape=prior_shape, prior_rate=prior_rate)
    _finite(failures=failures)
    _integer("failures", failures, 0)
    _nonnegative(total_time=total_time)
    return FailureRatePosterior(float(prior_shape + failures), float(prior_rate + total_time), float(prior_shape), float(prior_rate))


def gamma_prior_from_bounds(median, upper, confidence=0.95):
    r"""Gamma prior on a failure rate with a given median and upper percentile.

    An engineer can more easily say "the rate is probably about
    :math:`\lambda_{50}` and almost certainly below :math:`\lambda_{u}`"
    than choose a shape and rate. The ratio :math:`\lambda_u/\lambda_{50}`
    fixes the shape :math:`a` (it falls as :math:`a` grows), found by
    bisection, and then :math:`b = P^{-1}(a, 1/2)/\lambda_{50}`.

    Parameters
    ----------
    median : float
        Prior median failure rate [1/s], positive.
    upper : float
        Rate the prior puts below it with probability ``confidence``
        [1/s], greater than ``median``.
    confidence : float, optional
        Prior probability below ``upper`` [-], in (0.5, 1).

    Returns
    -------
    tuple of float
        ``(prior_shape, prior_rate)``: :math:`a` [-] and :math:`b` [s].

    Raises
    ------
    ValueError
        If the ratio is too wide to fit with a shape above :math:`10^{-3}`.

    Notes
    -----
    ``median``, ``upper`` and ``confidence`` are scalar-only, since the
    shape is a root found for one pair of bounds. Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/reliability`.

    Examples
    --------
    >>> shape, rate = gamma_prior_from_bounds(1e-6, 3e-6, 0.95)
    >>> round(shape, 3)
    1.782
    """
    _positive(median=median, upper=upper)
    _interval("confidence", confidence, 0.5, 1, closed=False)
    if upper <= median:
        raise ValueError("upper must exceed median")
    target = upper / median

    def ratio(shape):
        return special.gammaincinv(shape, confidence) / special.gammaincinv(shape, 0.5)

    low, high = 1e-3, 1.0
    if ratio(low) < target:
        raise ValueError("upper/median is too wide for a gamma prior")
    while ratio(high) > target:
        high *= 2
    for _ in range(200):
        middle = 0.5 * (low + high)
        low, high = (middle, high) if ratio(middle) > target else (low, middle)
    shape = 0.5 * (low + high)
    return float(shape), float(special.gammaincinv(shape, 0.5) / median)
