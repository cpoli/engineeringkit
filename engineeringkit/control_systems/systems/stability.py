r"""Algebraic and energy-based stability tests.

Maxwell (1868) reduced the stability of a governed engine to the signs of
the roots of its characteristic polynomial. Routh (1877) and Hurwitz (1895)
found coefficient tests that avoid computing the roots; Lyapunov (1892)
replaced roots with an energy-like function that must decrease along every
trajectory.
"""

from __future__ import annotations

import numpy as np
from scipy import linalg

from engineeringkit.control_systems.utils.matrices import state_matrices, weight_matrix
from engineeringkit.control_systems.utils.polynomials import as_polynomial

__all__ = ["routh_first_column", "routh_stable", "hurwitz_determinants", "lyapunov_equation", "lyapunov_stable"]


def routh_first_column(coeffs):
    """First column of the Routh array of a real polynomial.

    The number of sign changes down the column equals the number of roots in
    the open right half-plane.

    Parameters
    ----------
    coeffs : array_like
        Real coefficients in descending powers of :math:`s`, with a nonzero
        leading coefficient and degree at least one.

    Returns
    -------
    ndarray
        First-column entries for rows :math:`s^n` down to :math:`s^0`,
        normalised so the leading coefficient is one.

    Notes
    -----
    A row of zeros is replaced by the derivative of the auxiliary polynomial
    formed from the row above. A zero pivot is replaced by a small positive
    :math:`\\varepsilon`, so this is a numerical teaching method rather than a
    symbolic proof for roots on or near the imaginary axis. See
    :doc:`/guides/control_systems`.

    References
    ----------
    E. J. Routh, *A Treatise on the Stability of a Given State of Motion*
    (Macmillan, 1877).

    Examples
    --------
    >>> routh_first_column([1, 3, 3, 1]).round(4).tolist()
    [1.0, 3.0, 2.6667, 1.0]
    """
    c = np.asarray(coeffs, float)
    if c.ndim != 1 or len(c) < 2 or not np.all(np.isfinite(c)) or c[0] == 0:
        raise ValueError("coeffs must be finite with nonzero leading coefficient")
    c = c / c[0]
    n = len(c) - 1
    m = (len(c) + 1) // 2
    tab = np.zeros((n + 1, m))
    tab[0, : len(c[::2])] = c[::2]
    tab[1, : len(c[1::2])] = c[1::2]
    eps = 1e-12
    for i in range(2, n + 1):
        if np.all(tab[i - 1] == 0):
            powers = n - (i - 2) - 2 * np.arange(m)
            tab[i - 1] = tab[i - 2] * np.maximum(powers, 0)
        if tab[i - 1, 0] == 0:
            tab[i - 1, 0] = eps
        for j in range(m - 1):
            tab[i, j] = (tab[i - 1, 0] * tab[i - 2, j + 1] - tab[i - 2, 0] * tab[i - 1, j + 1]) / tab[i - 1, 0]
    return tab[:, 0]


def routh_stable(coeffs):
    """Whether every root of a polynomial lies strictly in the left half-plane.

    Parameters
    ----------
    coeffs : array_like
        Real coefficients in descending powers of :math:`s`.

    Returns
    -------
    bool
        ``True`` if the Routh column has no sign change and no root lies on the
        imaginary axis; marginal stability counts as unstable.

    Notes
    -----
    Routh signs alone cannot distinguish the auxiliary-row marginal case, so
    root locations confirm it with a relative tolerance. See
    :doc:`/guides/control_systems`.

    Examples
    --------
    >>> routh_stable([1, 3, 3, 1]), routh_stable([1, 0, 1])
    (True, False)
    """
    column = routh_first_column(coeffs)
    roots = np.roots(coeffs)
    tolerance = 1e-10 * max(1.0, float(np.max(np.abs(roots))))
    return bool(np.all(column > 0) and np.all(roots.real < -tolerance))


def hurwitz_determinants(coeffs):
    r"""Leading principal minors :math:`\Delta_1, \dots, \Delta_n` of the Hurwitz matrix.

    For :math:`a_0 s^n + a_1 s^{n-1} + \dots + a_n` with :math:`a_0 > 0`, all
    roots lie in the open left half-plane if and only if every
    :math:`\Delta_k > 0`. The Hurwitz matrix has entries
    :math:`H_{ij} = a_{2j-i}` (1-based, :math:`a_k = 0` outside
    :math:`0 \le k \le n`).

    Parameters
    ----------
    coeffs : array_like
        Real coefficients in descending powers of :math:`s`, degree at least
        one. They are normalised by the leading coefficient.

    Returns
    -------
    ndarray
        The :math:`n` determinants.

    Notes
    -----
    On the stability boundary :math:`\Delta_{n-1} = 0` (a pair of imaginary
    roots) or :math:`\Delta_n = 0` (a root at the origin). See
    :doc:`/guides/control_systems`.

    References
    ----------
    A. Hurwitz, "Ueber die Bedingungen, unter welchen eine Gleichung nur
    Wurzeln mit negativen reellen Theilen besitzt," *Math. Ann.* 46 (1895)
    273-284.

    Examples
    --------
    >>> hurwitz_determinants([1, 3, 3, 1]).round(6).tolist()
    [3.0, 8.0, 8.0]
    """
    a = as_polynomial(coeffs, "coeffs", min_degree=1)
    a = a / a[0]
    n = len(a) - 1
    idx = 2 * np.arange(1, n + 1)[None, :] - np.arange(1, n + 1)[:, None]
    H = np.where((idx >= 0) & (idx <= n), a[np.clip(idx, 0, n)], 0.0)
    return np.array([np.linalg.det(H[:k, :k]) for k in range(1, n + 1)])


def lyapunov_equation(A, Q=None):
    r"""Solve the continuous Lyapunov equation :math:`A^\mathsf{T} P + P A = -Q`.

    If :math:`A` is stable, :math:`V(x) = x^\mathsf{T} P x` is a Lyapunov
    function: it is positive and decreases at the rate
    :math:`\dot V = -x^\mathsf{T} Q x` along every trajectory of
    :math:`\dot x = Ax`.

    Parameters
    ----------
    A : array_like
        State matrix, shape ``(n, n)`` [1/s].
    Q : array_like, optional
        Symmetric positive semidefinite decay-rate matrix, shape ``(n, n)``;
        the identity by default.

    Returns
    -------
    ndarray
        Symmetric solution :math:`P`.

    Raises
    ------
    ValueError
        If two eigenvalues of :math:`A` sum to zero, so the equation has no
        unique solution.

    Notes
    -----
    See :doc:`/guides/control_systems`.

    References
    ----------
    A. M. Lyapunov, *The General Problem of the Stability of Motion* (Kharkov,
    1892; English translation, Taylor & Francis, 1992).

    Examples
    --------
    >>> lyapunov_equation([[-1.0]]).tolist()
    [[0.5]]
    """
    A, _, _ = state_matrices(A)
    Q = np.eye(len(A)) if Q is None else weight_matrix(Q, len(A), "Q")
    eig = np.linalg.eigvals(A)
    gap = np.min(np.abs(eig[:, None] + eig[None, :]))
    if gap <= 1e-10 * max(1.0, float(np.max(np.abs(eig)))):
        raise ValueError("eigenvalues of A sum to zero: the Lyapunov equation is singular")
    P = linalg.solve_continuous_lyapunov(A.T, -Q)
    return 0.5 * (P + P.T)


def lyapunov_stable(A):
    r"""Asymptotic stability of :math:`\dot x = Ax` by Lyapunov's direct method.

    Parameters
    ----------
    A : array_like
        State matrix, shape ``(n, n)`` [1/s].

    Returns
    -------
    bool
        ``True`` if :math:`A^\mathsf{T} P + PA = -I` has a positive definite
        solution :math:`P`.

    Notes
    -----
    A singular equation (eigenvalues summing to zero, e.g. an undamped
    oscillator) means the origin is not asymptotically stable. See
    :doc:`/guides/control_systems`.

    Examples
    --------
    >>> lyapunov_stable([[0, 1], [-2, -3]]), lyapunov_stable([[0, 1], [-1, 0]])
    (True, False)
    """
    try:
        P = lyapunov_equation(A)
    except ValueError as exc:
        if "singular" in str(exc):
            return False
        raise
    return bool(np.all(np.linalg.eigvalsh(P) > 0))
