r"""Structural reliability: stress-strength interference and the Hasofer-Lind index.

Freudenthal argued that a safety factor hides the scatter of loads and
strengths, and that the probability of the load exceeding the strength is
the quantity to design for. Hasofer and Lind defined a reliability index
that does not depend on how the limit state is written, the basis of the
first-order reliability method (FORM).
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.reliability.core.base import HasoferLindResult
from engineeringkit.reliability.utils.normal import failure_probability_from_index

__all__ = ["stress_strength_reliability", "hasofer_lind_index"]


def stress_strength_reliability(strength_mean, strength_std, stress_mean, stress_std):
    r"""Probability that a normal strength exceeds an independent normal stress.

    :math:`R = \Phi\left((\mu_S - \mu_s)/\sqrt{\sigma_S^2 + \sigma_s^2}\right)`.

    Parameters
    ----------
    strength_mean, strength_std : float or array_like
        Mean and standard deviation of the strength [Pa] (or any load unit);
        the standard deviation is nonnegative.
    stress_mean, stress_std : float or array_like
        Mean and standard deviation of the applied stress, in the same unit.

    Returns
    -------
    float or ndarray
        Reliability [-].

    Notes
    -----
    The two standard deviations may not both be zero. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/reliability`.

    References
    ----------
    A. M. Freudenthal, "The safety of structures," *Trans. ASCE* 112 (1947)
    125-180.

    Examples
    --------
    >>> round(float(stress_strength_reliability(500e6, 30e6, 400e6, 40e6)), 4)
    0.9772
    """
    _finite(strength_mean=strength_mean, stress_mean=stress_mean)
    _nonnegative(strength_std=strength_std, stress_std=stress_std)
    spread = np.hypot(strength_std, stress_std)
    _positive(combined_std=spread)
    beta = (np.asarray(strength_mean) - stress_mean) / spread
    return 1 - failure_probability_from_index(beta)


def hasofer_lind_index(limit_state, mean, std, tol=1e-8, max_iter=100):
    r"""Hasofer-Lind reliability index for independent normal variables (HL-RF iteration).

    Each variable is standardised, :math:`u_i = (x_i - \mu_i)/\sigma_i`, and
    :math:`\beta` is the distance from the origin to the nearest point of
    the failure surface :math:`g = 0`. The Rackwitz-Fiessler update
    :math:`u_{k+1} = [(\nabla g \cdot u_k - g)/|\nabla g|^2]\,\nabla g`
    finds it, with central-difference gradients.

    Parameters
    ----------
    limit_state : callable
        :math:`g(x)` for a 1-D array ``x``; failure is :math:`g < 0`.
    mean : array_like
        Means of the variables (their own units).
    std : array_like
        Standard deviations, positive, same length as ``mean``.
    tol : float, optional
        Convergence tolerance on the design point in standard space [-].
    max_iter : int, optional
        Maximum number of iterations.

    Returns
    -------
    HasoferLindResult
        Index, design point, sensitivity vector and iteration count. The
        index is negative if the mean point already lies in the failure
        domain.

    Raises
    ------
    ValueError
        If the gradient vanishes or the iteration does not converge.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    A. M. Hasofer and N. C. Lind, "Exact and invariant second-moment code
    format," *J. Eng. Mech. Div. ASCE* 100 (1974) 111-121; R. Rackwitz and
    B. Fiessler, "Structural reliability under combined random load
    sequences," *Comput. Struct.* 9 (1978) 489-494.

    Examples
    --------
    >>> result = hasofer_lind_index(lambda x: x[0] - x[1], [500.0, 400.0], [30.0, 40.0])
    >>> round(result.beta, 4)
    2.0
    """
    if not callable(limit_state):
        raise ValueError("limit_state must be callable")
    mu = _vector("mean", mean)
    sigma = _vector("std", std)
    if mu.shape != sigma.shape:
        raise ValueError("mean and std must have the same length")
    _positive(std=sigma, tol=tol)
    _integer("max_iter", max_iter, 1)

    def g(u):
        return float(limit_state(mu + sigma * u))

    def gradient(u, h=1e-5):
        return np.array([(g(u + h * e) - g(u - h * e)) / (2 * h) for e in np.eye(u.size)])

    u = np.zeros(mu.size)
    for iteration in range(1, max_iter + 1):
        grad = gradient(u)
        norm2 = float(grad @ grad)
        if norm2 == 0:
            raise ValueError("limit-state gradient vanished; the index is undefined")
        u_new = (grad @ u - g(u)) / norm2 * grad
        converged = np.linalg.norm(u_new - u) <= tol * (1 + np.linalg.norm(u_new))
        u = u_new
        if converged:
            break
    else:
        raise ValueError(f"HL-RF iteration did not converge in {max_iter} iterations")
    grad = gradient(u)
    alpha = -grad / np.linalg.norm(grad)
    beta = float(np.copysign(np.linalg.norm(u), g(np.zeros(mu.size))))
    return HasoferLindResult(beta, mu + sigma * u, alpha, iteration)
