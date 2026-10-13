r"""Optimal control and estimation: LQR, discrete LQR and the Kalman-Bucy filter.

Bellman's dynamic programming (1957) and Kalman's linear-quadratic regulator
(1960) turned controller design into minimising a quadratic cost, solved by
an algebraic Riccati equation. The same equation, transposed, gives the
steady-state Kalman-Bucy estimator (1961).
"""

from __future__ import annotations

import numpy as np
from scipy import linalg

from engineeringkit.control_systems.utils.matrices import state_matrices, weight_matrix

__all__ = ["lqr", "dlqr", "kalman_gain"]


def lqr(A, B, Q, R):
    r"""Continuous-time linear-quadratic regulator.

    Minimises :math:`J = \int_0^\infty (x^\mathsf{T} Q x + u^\mathsf{T} R u)\,dt`
    over :math:`u = -Kx` for :math:`\dot x = Ax + Bu`.

    Parameters
    ----------
    A : array_like
        State matrix, shape ``(n, n)`` [1/s].
    B : array_like
        Input matrix, shape ``(n, m)``.
    Q : array_like
        Symmetric positive semidefinite state weight, shape ``(n, n)``.
    R : array_like
        Symmetric positive definite input weight, shape ``(m, m)``.

    Returns
    -------
    tuple
        ``(K, P, poles)``: gain :math:`K = R^{-1} B^\mathsf{T} P`, the
        Riccati solution :math:`P` (so the optimal cost from :math:`x_0` is
        :math:`x_0^\mathsf{T} P x_0`) and the closed-loop poles [1/s].

    Notes
    -----
    :math:`P` solves :math:`A^\mathsf{T}P + PA - PBR^{-1}B^\mathsf{T}P + Q = 0`.
    A stabilising solution needs :math:`(A, B)` stabilisable and
    :math:`(Q^{1/2}, A)` detectable; SciPy errors are not masked. See
    :doc:`/guides/control_systems`.

    References
    ----------
    R. E. Kalman, "Contributions to the theory of optimal control," *Bol. Soc.
    Mat. Mexicana* 5 (1960) 102-119.

    Examples
    --------
    >>> K, P, poles = lqr([[0]], [[1]], [[4]], [[1]])
    >>> K.round(6).tolist(), poles.real.round(6).tolist()
    ([[2.0]], [-2.0])
    """
    A, B, _ = state_matrices(A, B=B)
    Q = weight_matrix(Q, len(A), "Q")
    R = weight_matrix(R, B.shape[1], "R", strict=True)
    P = linalg.solve_continuous_are(A, B, Q, R)
    K = np.linalg.solve(R, B.T @ P)
    return K, P, np.linalg.eigvals(A - B @ K)


def dlqr(A, B, Q, R):
    r"""Discrete-time infinite-horizon linear-quadratic regulator.

    Minimises :math:`\sum_{k \ge 0} (x_k^\mathsf{T} Q x_k + u_k^\mathsf{T} R u_k)`
    over :math:`u_k = -K x_k` for :math:`x_{k+1} = A x_k + B u_k`; the
    steady state of Bellman's backward dynamic-programming recursion.

    Parameters
    ----------
    A : array_like
        Discrete state matrix, shape ``(n, n)``.
    B : array_like
        Discrete input matrix, shape ``(n, m)``.
    Q : array_like
        Symmetric positive semidefinite state weight, shape ``(n, n)``.
    R : array_like
        Symmetric positive definite input weight, shape ``(m, m)``.

    Returns
    -------
    tuple
        ``(K, P, poles)``: gain :math:`K = (R + B^\mathsf{T}PB)^{-1}B^\mathsf{T}PA`,
        the discrete Riccati solution :math:`P` and the closed-loop
        eigenvalues (inside the unit circle when stabilising).

    Notes
    -----
    See :doc:`/guides/control_systems`.

    References
    ----------
    R. Bellman, *Dynamic Programming* (Princeton University Press, 1957).

    Examples
    --------
    >>> K, P, poles = dlqr([[1.0]], [[1.0]], [[1.0]], [[1.0]])
    >>> round(float(P[0, 0]), 6), round(float(K[0, 0]), 6)
    (1.618034, 0.618034)
    """
    A, B, _ = state_matrices(A, B=B)
    Q = weight_matrix(Q, len(A), "Q")
    R = weight_matrix(R, B.shape[1], "R", strict=True)
    P = linalg.solve_discrete_are(A, B, Q, R)
    K = np.linalg.solve(R + B.T @ P @ B, B.T @ P @ A)
    return K, P, np.linalg.eigvals(A - B @ K)


def kalman_gain(A, C, Q, R):
    r"""Steady-state Kalman-Bucy estimator gain for continuous-time measurements.

    For :math:`\dot x = Ax + w`, :math:`y = Cx + v` with white process and
    measurement noise of spectral densities :math:`Q` and :math:`R`, the
    estimator :math:`\dot{\hat x} = A\hat x + L(y - C\hat x)` with
    :math:`L = PC^\mathsf{T}R^{-1}` minimises the steady error covariance
    :math:`P`.

    Parameters
    ----------
    A : array_like
        State matrix, shape ``(n, n)`` [1/s].
    C : array_like
        Output matrix, shape ``(p, n)``.
    Q : array_like
        Process-noise spectral density, symmetric positive semidefinite,
        shape ``(n, n)``.
    R : array_like
        Measurement-noise spectral density, symmetric positive definite,
        shape ``(p, p)``.

    Returns
    -------
    tuple
        ``(L, P, poles)``: observer gain, error covariance and estimator poles
        [1/s].

    Notes
    -----
    :math:`P` solves :math:`AP + PA^\mathsf{T} - PC^\mathsf{T}R^{-1}CP + Q = 0`.
    :math:`Q` and :math:`R` are spectral densities, not sampled covariances.
    See :doc:`/guides/control_systems`.

    References
    ----------
    R. E. Kalman and R. S. Bucy, "New results in linear filtering and
    prediction theory," *J. Basic Eng.* 83 (1961) 95-108.

    Examples
    --------
    >>> L, P, poles = kalman_gain([[0]], [[1]], [[4]], [[1]])
    >>> L.round(6).tolist(), P.round(6).tolist()
    ([[2.0]], [[2.0]])
    """
    A, _, C = state_matrices(A, C=C)
    Q = weight_matrix(Q, len(A), "Q")
    R = weight_matrix(R, C.shape[0], "R", strict=True)
    P = linalg.solve_continuous_are(A.T, C.T, Q, R)
    L = np.linalg.solve(R, C @ P).T
    return L, P, np.linalg.eigvals(A - L @ C)
