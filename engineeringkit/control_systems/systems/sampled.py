r"""Sampled-data control: exact zero-order-hold discretisation.

Ragazzini and Zadeh (1952) analysed loops in which a computer samples the
plant output and holds its command constant between samples. For a
zero-order hold of period :math:`T`, the continuous model
:math:`\dot x = Ax + Bu` becomes exactly
:math:`x_{k+1} = A_d x_k + B_d u_k` with :math:`A_d = e^{AT}` and
:math:`B_d = \int_0^T e^{A\tau}\,d\tau\,B`.
"""

from __future__ import annotations

import numpy as np
from scipy import linalg

from engineeringkit._validation import positive as _positive
from engineeringkit.control_systems.utils.matrices import state_matrices

__all__ = ["discretize_zoh"]


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
