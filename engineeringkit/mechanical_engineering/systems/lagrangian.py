r"""Lagrange's equations: equations of motion from a single scalar function.

Lagrange's *Mécanique analytique* wrote the motion of any system with
generalised coordinates :math:`q` and Lagrangian :math:`L = T - V` as

.. math:: \frac{d}{dt}\frac{\partial L}{\partial \dot q} - \frac{\partial L}{\partial q} = Q

with :math:`Q` the generalised non-conservative forces. Expanding the time
derivative gives :math:`H\ddot q = \partial L/\partial q - J\dot q + Q`, with
:math:`H = \partial^2 L/\partial\dot q^2` (the mass matrix) and
:math:`J = \partial^2 L/\partial\dot q\,\partial q`. Here the derivatives are
taken by central differences, so a model is just a Python function for
:math:`L`.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.integrators import integrate_ode
from engineeringkit.mechanical_engineering.core.base import MDOFResult

__all__ = ["lagrange_accelerations", "simulate_lagrangian"]


def _accelerations(lagrangian, q, qdot, Q, h):
    n = q.size
    eye = np.eye(n)
    hq, hv = h * (1 + np.abs(q)), h * (1 + np.abs(qdot))

    def L(a, b):
        return float(lagrangian(a, b))

    grad = np.array([(L(q + hq[i] * eye[i], qdot) - L(q - hq[i] * eye[i], qdot)) / (2 * hq[i]) for i in range(n)])
    H = np.empty((n, n))
    J = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            di, dj = hv[i] * eye[i], hv[j] * eye[j]
            H[i, j] = (L(q, qdot + di + dj) - L(q, qdot + di - dj) - L(q, qdot - di + dj) + L(q, qdot - di - dj)) / (4 * hv[i] * hv[j])
            dq = hq[j] * eye[j]
            J[i, j] = (L(q + dq, qdot + di) - L(q + dq, qdot - di) - L(q - dq, qdot + di) + L(q - dq, qdot - di)) / (4 * hv[i] * hq[j])
    return np.linalg.solve(H, grad - J @ qdot + Q)


def _forces(forces, t, q, qdot):
    if forces is None:
        return np.zeros(q.size)
    Q = np.asarray(forces(t, q, qdot), dtype=float).reshape(-1)
    if Q.size != q.size or not np.all(np.isfinite(Q)):
        raise ValueError("forces must return one finite generalised force per coordinate")
    return Q


def lagrange_accelerations(lagrangian, q, qdot, t=0.0, forces=None, step=1e-4):
    r"""Generalised accelerations :math:`\ddot q` from Lagrange's equations at one state.

    Parameters
    ----------
    lagrangian : callable
        ``lagrangian(q, qdot) -> float``, the Lagrangian :math:`T - V` [J]
        of the coordinate and rate vectors.
    q : array_like
        Generalised coordinates [m or rad].
    qdot : array_like
        Their rates [m/s or rad/s], same length as ``q``.
    t : float, optional
        Time [s], passed to ``forces``.
    forces : callable, optional
        ``forces(t, q, qdot) -> Q``, the generalised non-conservative forces
        [N or N m], one per coordinate; none if omitted.
    step : float, optional
        Relative central-difference step [-], positive.

    Returns
    -------
    ndarray
        :math:`\ddot q` [m/s^2 or rad/s^2].

    Raises
    ------
    numpy.linalg.LinAlgError
        If :math:`\partial^2 L/\partial\dot q^2` is singular, so the
        coordinates are not independent.

    Notes
    -----
    The derivatives are second-order central differences, exact to
    rounding for a Lagrangian quadratic in ``qdot`` and accurate to about
    ``step**2`` in ``q``. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/mechanical_engineering`.

    References
    ----------
    J. L. Lagrange, *Mécanique analytique* (Desaint, Paris, 1788), part II,
    section IV.

    Examples
    --------
    A pendulum of length 1 m, :math:`L = \tfrac12\dot\theta^2 + g\cos\theta`:

    >>> L = lambda q, v: 0.5 * v[0] ** 2 + 9.81 * np.cos(q[0])
    >>> round(float(lagrange_accelerations(L, [np.pi / 2], [0.0])[0]), 6)
    -9.81
    """
    qq, vv = _vector("q", q), _vector("qdot", qdot)
    if qq.size != vv.size:
        raise ValueError("q and qdot must have the same length")
    _finite(t=t)
    _positive(step=step)
    return _accelerations(lagrangian, qq, vv, _forces(forces, t, qq, vv), step)


def simulate_lagrangian(lagrangian, q0, qdot0, t=None, forces=None, step=1e-4) -> MDOFResult:
    r"""Integrate Lagrange's equations of motion from an initial state.

    Parameters
    ----------
    lagrangian : callable
        ``lagrangian(q, qdot) -> float`` [J]; see
        :func:`lagrange_accelerations`.
    q0, qdot0 : array_like
        Initial coordinates [m or rad] and rates [m/s or rad/s].
    t : array_like, optional
        Strictly increasing output times [s]; 1001 samples on [0, 10] s if
        omitted.
    forces : callable, optional
        ``forces(t, q, qdot) -> Q`` [N or N m]; none if omitted.
    step : float, optional
        Relative central-difference step [-], positive.

    Returns
    -------
    MDOFResult
        Coordinates and rates at each time, shape ``(len(t), n)``.

    Notes
    -----
    Uses :func:`engineeringkit.integrators.integrate_ode` (adaptive
    Runge-Kutta). Inputs are raw SI numbers with no dimensional metadata;
    see :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> L = lambda q, v: 0.5 * v[0] ** 2 - 0.5 * 4.0 * q[0] ** 2
    >>> result = simulate_lagrangian(L, [0.1], [0.0], t=np.linspace(0, np.pi, 3))
    >>> np.round(result.displacement[:, 0], 6).tolist()
    [0.1, -0.1, 0.1]
    """
    qq, vv = _vector("q0", q0), _vector("qdot0", qdot0)
    if qq.size != vv.size:
        raise ValueError("q0 and qdot0 must have the same length")
    _positive(step=step)
    n = qq.size
    times = np.linspace(0, 10, 1001) if t is None else t

    def rhs(tt, y):
        q, v = y[:n], y[n:]
        return np.concatenate((v, _accelerations(lagrangian, q, v, _forces(forces, tt, q, v), step)))

    result = integrate_ode(rhs, np.concatenate((qq, vv)), times)
    return MDOFResult(t=result.t, displacement=result.history[:, :n], velocity=result.history[:, n:])
