"""Shape, finiteness and definiteness checks for state-space and weighting matrices.

State-space models follow :math:`\\dot x = Ax + Bu`, :math:`y = Cx + Du` with
``A`` of shape ``(n, n)``, ``B`` of shape ``(n, m)`` and ``C`` of shape
``(p, n)``. SISO models still use explicit two-dimensional ``B`` and ``C``.
"""

from __future__ import annotations

import numpy as np

__all__ = ["state_matrices", "weight_matrix"]


def state_matrices(A, B=None, C=None):
    """Validate and convert state-space matrices to float arrays.

    Parameters
    ----------
    A : array_like
        Square state matrix, shape ``(n, n)`` [1/s].
    B : array_like, optional
        Input matrix, shape ``(n, m)``.
    C : array_like, optional
        Output matrix, shape ``(p, n)``.

    Returns
    -------
    tuple of ndarray
        ``(A, B, C)``, with ``None`` for matrices not supplied.

    Raises
    ------
    ValueError
        If a matrix is nonfinite, empty or has an inconsistent shape.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/control_systems`.

    Examples
    --------
    >>> A, B, C = state_matrices([[0, 1], [-2, -3]], B=[[0], [1]])
    >>> A.shape, B.shape, C
    ((2, 2), (2, 1), None)
    """
    A = np.asarray(A, float)
    if A.ndim != 2 or not A.shape[0] or A.shape[0] != A.shape[1] or not np.all(np.isfinite(A)):
        raise ValueError("A must be a finite nonempty square matrix")
    n = A.shape[0]
    if B is not None:
        B = np.asarray(B, float)
        if B.ndim != 2 or B.shape[0] != n or not B.shape[1] or not np.all(np.isfinite(B)):
            raise ValueError("B must have shape (states, inputs)")
    if C is not None:
        C = np.asarray(C, float)
        if C.ndim != 2 or C.shape[1] != n or not C.shape[0] or not np.all(np.isfinite(C)):
            raise ValueError("C must have shape (outputs, states)")
    return A, B, C


def weight_matrix(matrix, size, name, strict=False):
    """Validate a symmetric positive (semi)definite weighting or covariance matrix.

    Parameters
    ----------
    matrix : array_like
        Candidate matrix, shape ``(size, size)``.
    size : int
        Required dimension.
    name : str
        Argument name used in error messages.
    strict : bool, optional
        Require positive *definite* (``True``) rather than semidefinite.

    Returns
    -------
    ndarray
        The matrix as a float array.

    Raises
    ------
    ValueError
        If the matrix has the wrong shape, is not symmetric, or is indefinite.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/control_systems`.

    Examples
    --------
    >>> weight_matrix([[2.0]], 1, "R", strict=True)
    array([[2.]])
    """
    matrix = np.asarray(matrix, float)
    if matrix.shape != (size, size) or not np.all(np.isfinite(matrix)):
        raise ValueError(f"{name} must be finite with shape {(size, size)}")
    if not np.allclose(matrix, matrix.T, rtol=1e-10, atol=1e-12):
        raise ValueError(f"{name} must be symmetric")
    eig = np.linalg.eigvalsh(matrix)
    if np.any(eig <= 0) if strict else np.any(eig < -1e-12 * max(1.0, np.max(abs(eig)))):
        raise ValueError(f"{name} must be positive {'definite' if strict else 'semidefinite'}")
    return matrix
