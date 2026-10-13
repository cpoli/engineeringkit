r"""State-space models: controllability, observability, pole placement and observers.

Kalman (1960) showed that the internal state of :math:`\dot x = Ax + Bu`,
:math:`y = Cx + Du` can be fully steered only if the controllability matrix
has full rank, and fully reconstructed only if the observability matrix
does. With both, state feedback places the closed-loop poles anywhere
(Wonham, 1967) and a Luenberger observer (1964) supplies the unmeasured
states.
"""

from __future__ import annotations

import numpy as np
from scipy import signal

from engineeringkit._validation import finite as _finite
from engineeringkit.control_systems.utils.matrices import state_matrices

__all__ = [
    "state_space",
    "controllability_matrix",
    "observability_matrix",
    "is_controllable",
    "is_observable",
    "state_feedback_poles",
    "observer_gain",
    "feedforward_gain",
]


def state_space(A, B, C, D):
    """Continuous-time state-space model :math:`\\dot x = Ax + Bu,\\ y = Cx + Du`.

    Parameters
    ----------
    A : array_like
        State matrix, shape ``(n, n)`` [1/s].
    B : array_like
        Input matrix, shape ``(n, m)``.
    C : array_like
        Output matrix, shape ``(p, n)``.
    D : array_like
        Feedthrough matrix, shape ``(p, m)``.

    Returns
    -------
    scipy.signal.StateSpace
        Continuous-time model.

    Notes
    -----
    Matrix units follow the chosen state, input and output units; see
    :doc:`/guides/control_systems`.

    Examples
    --------
    >>> sys = state_space([[-2]], [[2]], [[1]], [[0]])
    >>> float((-sys.C @ np.linalg.solve(sys.A, sys.B) + sys.D).item())
    1.0
    """
    A, B, C = state_matrices(A, B, C)
    D = np.asarray(D, float)
    if D.shape != (C.shape[0], B.shape[1]) or not np.all(np.isfinite(D)):
        raise ValueError("D must have shape (outputs, inputs)")
    return signal.StateSpace(A, B, C, D)


def controllability_matrix(A, B):
    r"""Kalman controllability matrix :math:`[B,\ AB,\ \dots,\ A^{n-1}B]`.

    Parameters
    ----------
    A : array_like
        State matrix, shape ``(n, n)``.
    B : array_like
        Input matrix, shape ``(n, m)``.

    Returns
    -------
    ndarray
        Shape ``(n, n m)``.

    Notes
    -----
    See :doc:`/guides/control_systems`.

    References
    ----------
    R. E. Kalman, "On the general theory of control systems," *Proc. 1st IFAC
    Congress* (Moscow, 1960) 481-492.

    Examples
    --------
    >>> controllability_matrix([[0, 1], [0, 0]], [[0], [1]]).tolist()
    [[0.0, 1.0], [1.0, 0.0]]
    """
    A, B, _ = state_matrices(A, B=B)
    return np.hstack([np.linalg.matrix_power(A, k) @ B for k in range(A.shape[0])])


def observability_matrix(A, C):
    r"""Kalman observability matrix :math:`[C;\ CA;\ \dots;\ CA^{n-1}]`.

    Parameters
    ----------
    A : array_like
        State matrix, shape ``(n, n)``.
    C : array_like
        Output matrix, shape ``(p, n)``.

    Returns
    -------
    ndarray
        Shape ``(n p, n)``.

    Notes
    -----
    See :doc:`/guides/control_systems`.

    Examples
    --------
    >>> observability_matrix([[0, 1], [0, 0]], [[1, 0]]).tolist()
    [[1.0, 0.0], [0.0, 1.0]]
    """
    A, _, C = state_matrices(A, C=C)
    return np.vstack([C @ np.linalg.matrix_power(A, k) for k in range(A.shape[0])])


def is_controllable(A, B):
    """Whether the controllability matrix has full state rank.

    Parameters
    ----------
    A : array_like
        State matrix, shape ``(n, n)``.
    B : array_like
        Input matrix, shape ``(n, m)``.

    Returns
    -------
    bool
        ``True`` if every state can be steered by the inputs.

    Notes
    -----
    Rank is computed numerically with SVD tolerances. See
    :doc:`/guides/control_systems`.

    Examples
    --------
    >>> is_controllable([[-1, 0], [0, -1]], [[1], [1]])
    False
    """
    return bool(np.linalg.matrix_rank(controllability_matrix(A, B)) == np.asarray(A).shape[0])


def is_observable(A, C):
    """Whether the observability matrix has full state rank.

    Parameters
    ----------
    A : array_like
        State matrix, shape ``(n, n)``.
    C : array_like
        Output matrix, shape ``(p, n)``.

    Returns
    -------
    bool
        ``True`` if every state can be reconstructed from the outputs.

    Notes
    -----
    See :doc:`/guides/control_systems`.

    Examples
    --------
    >>> is_observable([[0, 1], [0, 0]], [[1, 0]])
    True
    """
    return bool(np.linalg.matrix_rank(observability_matrix(A, C)) == np.asarray(A).shape[0])


def state_feedback_poles(A, B, desired_poles):
    """State-feedback gain :math:`K` placing the eigenvalues of :math:`A - BK`.

    Parameters
    ----------
    A : array_like
        State matrix, shape ``(n, n)``.
    B : array_like
        Input matrix, shape ``(n, m)``; ``(A, B)`` must be controllable.
    desired_poles : array_like
        ``n`` closed-loop poles [1/s]; complex poles must come in conjugate
        pairs.

    Returns
    -------
    ndarray
        Gain :math:`K`, shape ``(m, n)``, for the control law :math:`u = -Kx`.

    Notes
    -----
    Wraps :func:`scipy.signal.place_poles`. See :doc:`/guides/control_systems`.

    References
    ----------
    W. M. Wonham, "On pole assignment in multi-input controllable linear
    systems," *IEEE Trans. Autom. Control* 12 (1967) 660-665.

    Examples
    --------
    >>> state_feedback_poles([[0, 1], [0, 0]], [[0], [1]], [-1, -2]).round(6).tolist()
    [[2.0, 3.0]]
    """
    A, B, _ = state_matrices(A, B=B)
    return signal.place_poles(A, B, desired_poles).gain_matrix


def observer_gain(A, C, desired_poles):
    r"""Luenberger observer gain :math:`L` placing the eigenvalues of :math:`A - LC`.

    The observer :math:`\dot{\hat x} = A\hat x + Bu + L(y - C\hat x)` makes
    the estimation error obey :math:`\dot e = (A - LC)e`.

    Parameters
    ----------
    A : array_like
        State matrix, shape ``(n, n)``.
    C : array_like
        Output matrix, shape ``(p, n)``; ``(A, C)`` must be observable.
    desired_poles : array_like
        ``n`` error-dynamics poles [1/s], conjugate pairs for complex values.

    Returns
    -------
    ndarray
        Gain :math:`L`, shape ``(n, p)``.

    Notes
    -----
    Computed by duality, placing the poles of :math:`A^\mathsf{T} -
    C^\mathsf{T} L^\mathsf{T}`. See :doc:`/guides/control_systems`.

    References
    ----------
    D. G. Luenberger, "Observing the state of a linear system," *IEEE Trans.
    Mil. Electron.* 8 (1964) 74-80.

    Examples
    --------
    >>> observer_gain([[0, 1], [0, 0]], [[1, 0]], [-1, -2]).round(6).tolist()
    [[3.0], [2.0]]
    """
    A, _, C = state_matrices(A, C=C)
    return signal.place_poles(A.T, C.T, desired_poles).gain_matrix.T


def feedforward_gain(A, B, C, K):
    r"""Reference precompensator :math:`N` giving unit DC gain under SISO state feedback.

    With :math:`u = -Kx + Nr`, the closed loop
    :math:`\dot x = (A - BK)x + BNr` tracks a constant reference when
    :math:`N = -1/\left(C (A - BK)^{-1} B\right)`.

    Parameters
    ----------
    A : array_like
        State matrix, shape ``(n, n)``.
    B : array_like
        Input matrix, shape ``(n, 1)``.
    C : array_like
        Output matrix, shape ``(1, n)``.
    K : array_like
        State-feedback gain, shape ``(1, n)``.

    Returns
    -------
    float
        Precompensator gain [input/reference units].

    Raises
    ------
    ValueError
        If the system is not SISO, or the closed loop has a zero at the
        origin so no finite :math:`N` exists.

    Notes
    -----
    Exact only for a perfect model; integral action is needed for robust
    tracking. See :doc:`/guides/control_systems`.

    Examples
    --------
    >>> feedforward_gain([[-1]], [[1]], [[1]], [[2]])
    3.0
    """
    A, B, C = state_matrices(A, B, C)
    K = np.asarray(K, float)
    _finite(K=K)
    if B.shape[1] != 1 or C.shape[0] != 1 or K.shape != (1, A.shape[0]):
        raise ValueError("feedforward_gain requires a SISO system with K of shape (1, states)")
    dc = (C @ np.linalg.solve(A - B @ K, B)).item()
    if dc == 0:
        raise ValueError("the closed loop has zero DC gain; no finite precompensator exists")
    return float(-1 / dc)
