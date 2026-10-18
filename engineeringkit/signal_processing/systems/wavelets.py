r"""The Haar wavelet transform and its inverse.

Haar (1910) built the first orthonormal basis of square-integrable
functions from a single step: averages and differences of neighbouring
samples, repeated on the averages. Unlike a Fourier basis, each Haar
function is local in time, so a jump in the signal touches only a few
coefficients. The orthonormal scaling :math:`1/\sqrt2` preserves energy,
and the transform is inverted exactly.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer

__all__ = ["haar_transform", "inverse_haar_transform"]


def haar_transform(x, levels=None):
    r"""Orthonormal multilevel Haar wavelet transform.

    Parameters
    ----------
    x : array_like
        Real 1D samples [signal units]; the length must be divisible by
        :math:`2^{levels}`.
    levels : int, optional
        Number of averaging and differencing stages, at least 1. Defaults to
        the largest number the length allows.

    Returns
    -------
    list of ndarray
        ``[a_J, d_J, d_{J-1}, ..., d_1]``: the coarsest approximation and the
        details from coarse to fine [signal units], the order used by
        PyWavelets' ``wavedec``.

    Notes
    -----
    Each stage maps pairs to :math:`a = (x_{2k} + x_{2k+1})/\sqrt2` and
    :math:`d = (x_{2k} - x_{2k+1})/\sqrt2`, then repeats on :math:`a`. The
    sum of squared coefficients equals that of ``x``. See
    :doc:`/guides/signal_processing`.

    References
    ----------
    A. Haar, "Zur Theorie der orthogonalen Funktionensysteme," *Math. Ann.*
    69 (1910) 331-371.

    Examples
    --------
    >>> [c.round(6).tolist() for c in haar_transform([4.0, 2.0, 5.0, 5.0])]
    [[8.0], [-2.0], [1.414214, 0.0]]
    """
    x = np.asarray(x)
    _finite(x=x)
    if x.ndim != 1 or np.iscomplexobj(x) or x.size < 2:
        raise ValueError("x must be a real 1D signal with at least 2 samples")
    max_levels = 0
    while x.size % 2 ** (max_levels + 1) == 0:
        max_levels += 1
    levels = max_levels if levels is None else levels
    _integer("levels", levels, 1)
    if levels > max_levels:
        raise ValueError(f"len(x) = {x.size} allows at most {max_levels} levels")
    approximation, details = x.astype(float), []
    for _ in range(levels):
        even, odd = approximation[0::2], approximation[1::2]
        details.append((even - odd) / np.sqrt(2))
        approximation = (even + odd) / np.sqrt(2)
    return [approximation, *details[::-1]]


def inverse_haar_transform(coefficients):
    r"""Rebuild a signal from :func:`haar_transform` coefficients.

    Parameters
    ----------
    coefficients : sequence of array_like
        ``[a_J, d_J, ..., d_1]``, each detail twice as long as the one before
        it, the first detail as long as ``a_J`` [signal units].

    Returns
    -------
    ndarray
        Samples [signal units].

    Notes
    -----
    Inverts each stage with :math:`x_{2k} = (a + d)/\sqrt2` and
    :math:`x_{2k+1} = (a - d)/\sqrt2`. Zeroing small details before
    inverting is the simplest wavelet compression or denoising. See
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> inverse_haar_transform([[8.0], [-2.0], [np.sqrt(2), 0.0]]).round(12).tolist()
    [4.0, 2.0, 5.0, 5.0]
    """
    if len(coefficients) < 2:
        raise ValueError("coefficients must hold an approximation and at least one detail")
    approximation = np.asarray(coefficients[0], dtype=float)
    _finite(approximation=approximation)
    for detail in coefficients[1:]:
        detail = np.asarray(detail, dtype=float)
        _finite(detail=detail)
        if detail.shape != approximation.shape or detail.ndim != 1:
            raise ValueError("each detail must match the length of the approximation it refines")
        x = np.empty(2 * approximation.size)
        x[0::2] = (approximation + detail) / np.sqrt(2)
        x[1::2] = (approximation - detail) / np.sqrt(2)
        approximation = x
    return approximation
