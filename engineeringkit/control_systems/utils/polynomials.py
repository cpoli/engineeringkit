"""Polynomial and transfer-function input handling.

Polynomials are coefficient vectors in descending powers of the Laplace
variable :math:`s`, as used by :func:`numpy.polyval` and
:class:`scipy.signal.TransferFunction`.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy import signal

__all__ = ["as_polynomial", "num_den"]


def as_polynomial(coeffs, name="coeffs", min_degree=0):
    """Validate a descending-power polynomial and strip leading zeros.

    Parameters
    ----------
    coeffs : array_like
        Real coefficients, highest power first.
    name : str, optional
        Argument name used in error messages.
    min_degree : int, optional
        Smallest admissible degree after stripping leading zeros.

    Returns
    -------
    ndarray
        Coefficients with a nonzero leading term.

    Raises
    ------
    ValueError
        If the coefficients are not a finite real vector, are all zero, or the
        degree is below ``min_degree``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/control_systems`.

    Examples
    --------
    >>> as_polynomial([0, 1, 3, 2])
    array([1., 3., 2.])
    """
    c = np.asarray(coeffs, float)
    if c.ndim != 1 or not c.size or not np.all(np.isfinite(c)):
        raise ValueError(f"{name} must be a finite nonempty coefficient vector")
    nonzero = np.flatnonzero(c)
    if not nonzero.size:
        raise ValueError(f"{name} must not be identically zero")
    c = c[nonzero[0] :]
    if len(c) - 1 < min_degree:
        raise ValueError(f"{name} must have degree at least {min_degree}")
    return c


def num_den(sys):
    """Numerator and denominator polynomials of a SISO continuous-time model.

    Parameters
    ----------
    sys : scipy.signal.lti
        A :class:`~scipy.signal.TransferFunction`, or any SciPy LTI model that
        converts to one (state space, zeros-poles-gain).

    Returns
    -------
    tuple of ndarray
        ``(num, den)`` in descending powers of :math:`s`.

    Raises
    ------
    ValueError
        If the model is not single-input single-output.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/control_systems`.

    Examples
    --------
    >>> num_den(signal.TransferFunction([2, 2], [1, 3, 2]))
    (array([2., 2.]), array([1., 3., 2.]))
    """
    if not hasattr(sys, "den"):
        with warnings.catch_warnings():
            # State-space conversion pads the numerator with exact zeros.
            warnings.simplefilter("ignore", signal.BadCoefficients)
            sys = sys.to_tf()
    num, den = np.atleast_1d(np.squeeze(sys.num)), np.atleast_1d(np.squeeze(sys.den))
    if num.ndim != 1 or den.ndim != 1:
        raise ValueError("only single-input single-output models are supported")
    return as_polynomial(num, "num"), as_polynomial(den, "den")
