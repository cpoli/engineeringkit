r"""Power-law fits for the empirical laws of materials testing.

Basquin's S-N curve, the Coffin-Manson relation, Paris's law, Norton's creep
law and the Hall-Petch relation (in :math:`d^{-1/2}`) are all straight lines
on log-log axes. A least-squares line through :math:`(\log x, \log y)` gives
the exponent and coefficient that the test data imply.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import positive as _positive

__all__ = ["power_law_fit"]


def power_law_fit(x, y):
    r"""Least-squares fit of :math:`y = A\,x^p` on log-log axes.

    Parameters
    ----------
    x, y : array_like
        Matching vectors of at least two positive samples, with at least two
        distinct ``x``.

    Returns
    -------
    tuple of float
        ``(A, p)``, the coefficient (in the units of ``y`` over those of
        :math:`x^p`) and the exponent.

    Notes
    -----
    The fit minimises the error in :math:`\log y`, so every decade of data
    carries the same weight. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/materials`.

    Examples
    --------
    >>> A, p = power_law_fit([1.0, 10.0, 100.0], [3.0, 30.0, 300.0])
    >>> round(A, 9), round(p, 9)
    (3.0, 1.0)
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.ndim != 1 or x.shape != y.shape or x.size < 2:
        raise ValueError("x and y must be matching vectors of at least two samples")
    _positive(x=x, y=y)
    if np.ptp(x) == 0:
        raise ValueError("x must contain at least two distinct values")
    p, log_A = np.polyfit(np.log(x), np.log(y), 1)
    return float(np.exp(log_A)), float(p)
