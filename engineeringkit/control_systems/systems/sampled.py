r"""Sampled-data control: exact zero-order-hold discretisation and Jury's test.

Ragazzini and Zadeh (1952) analysed loops in which a computer samples the
plant output and holds its command constant between samples. For a
zero-order hold of period :math:`T`, the continuous model
:math:`\dot x = Ax + Bu` becomes exactly
:math:`x_{k+1} = A_d x_k + B_d u_k` with :math:`A_d = e^{AT}` and
:math:`B_d = \int_0^T e^{A\tau}\,d\tau\,B`. The sampled loop is stable
when every root of its characteristic polynomial in :math:`z` lies inside
the unit circle, which Jury (1961) tested from the coefficients alone.
"""

from __future__ import annotations

import numpy as np
from scipy import linalg

from engineeringkit._validation import positive as _positive
from engineeringkit.control_systems.utils.matrices import state_matrices
from engineeringkit.control_systems.utils.polynomials import as_polynomial

__all__ = ["discretize_zoh", "jury_stable"]


def discretize_zoh(A, B, dt):
    r"""Exact zero-order-hold discretisation :math:`(A_d, B_d)` of :math:`\dot x = Ax + Bu`.

    Parameters
    ----------
    A : array_like
        Continuous state matrix, shape ``(n, n)`` [1/s].
    B : array_like
        Continuous input matrix, shape ``(n, m)``.
    dt : float
        Sampling period [s], positive.

    Returns
    -------
    tuple of ndarray
        ``(Ad, Bd)`` for :math:`x_{k+1} = A_d x_k + B_d u_k`.

    Notes
    -----
    Both matrices come from one matrix exponential of the augmented matrix
    :math:`\begin{bmatrix}A & B\\ 0 & 0\end{bmatrix} T`, which is valid even
    when :math:`A` is singular. See :doc:`/guides/control_systems`.

    References
    ----------
    J. R. Ragazzini and L. A. Zadeh, "The analysis of sampled-data systems,"
    *Trans. AIEE* 71 (1952) 225-234.

    Examples
    --------
    >>> Ad, Bd = discretize_zoh([[0.0]], [[1.0]], 0.5)
    >>> Ad.tolist(), Bd.tolist()
    ([[1.0]], [[0.5]])
    """
    A, B, _ = state_matrices(A, B=B)
    _positive(dt=dt)
    n, m = B.shape
    M = np.zeros((n + m, n + m))
    M[:n, :n], M[:n, n:] = A, B
    E = linalg.expm(M * dt)
    return E[:n, :n], E[:n, n:]


def jury_stable(coeffs):
    r"""Whether every root of a polynomial in :math:`z` lies strictly inside the unit circle.

    Jury's table reduces :math:`P(z) = a_0 z^n + a_1 z^{n-1} + \dots + a_n`
    one degree at a time. Each row is
    :math:`a_0 P(z) - a_n z^n P(1/z)`, divided by :math:`z`, and the
    polynomial is stable if and only if :math:`|a_n| < |a_0|` holds at every
    step down to degree zero.

    Parameters
    ----------
    coeffs : array_like
        Real coefficients in descending powers of :math:`z`, degree at least
        one.

    Returns
    -------
    bool
        ``True`` if all roots satisfy :math:`|z| < 1`; roots on the unit
        circle count as unstable.

    Notes
    -----
    This is the discrete counterpart of Routh's test: no roots are
    computed. The quick necessary conditions :math:`P(1) > 0` and
    :math:`(-1)^n P(-1) > 0` (for :math:`a_0 > 0`) follow from it. Each row
    is normalised to unit leading coefficient, and a relative tolerance of
    :math:`10^{-12}` decides the boundary cases. See
    :doc:`/guides/control_systems`.

    References
    ----------
    E. I. Jury and J. Blanchard, "A stability test for linear discrete
    systems in table form," *Proc. IRE* 49 (1961) 1947-1948.

    Examples
    --------
    >>> jury_stable([1, -1.0, 0.5]), jury_stable([1, -2.5, 1.0])
    (True, False)
    """
    a = as_polynomial(coeffs, "coeffs", min_degree=1)
    while len(a) > 1:
        a = a / a[0]
        if abs(a[-1]) >= 1 - 1e-12:
            return False
        a = (a - a[-1] * a[::-1])[:-1]
    return True
