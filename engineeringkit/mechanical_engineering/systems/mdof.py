r"""Multi-degree-of-freedom vibration: normal modes, modal superposition and frequency estimates.

Daniel Bernoulli held that a vibrating string sounds a superposition of simple
modes, each at its own frequency. For a system of masses and springs,
:math:`M\ddot x + Kx = 0`, the modes are the eigenvectors of
:math:`K\phi = \omega^2 M\phi`; they are orthogonal through :math:`M` and
:math:`K`, so any motion splits into independent single-degree-of-freedom
oscillators. Rayleigh's quotient bounds the fundamental frequency from above
for any trial shape, and Dunkerley's sum bounds it from below.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.mechanical_engineering.core.base import MDOFResult, ModalResult

__all__ = [
    "normal_modes",
    "modal_response",
    "harmonic_response",
    "rayleigh_frequency",
    "dunkerley_frequency",
]


def _square(name, value, size=None):
    a = np.asarray(value)
    if a.ndim != 2 or a.shape[0] != a.shape[1] or not a.size or (size is not None and a.shape[0] != size):
        raise ValueError(f"{name} must be a square matrix" + ("" if size is None else f" of size {size}"))
    _finite(**{name: a})
    a = a.astype(float)
    if not np.allclose(a, a.T, rtol=1e-12, atol=1e-12 * np.max(np.abs(a))):
        raise ValueError(f"{name} must be symmetric")
    return a


def _system(M, K):
    mass = _square("M", M)
    stiffness = _square("K", K, mass.shape[0])
    if np.linalg.eigvalsh(mass)[0] <= 0:
        raise ValueError("M must be positive definite")
    if np.linalg.eigvalsh(stiffness)[0] < -1e-10 * np.max(np.abs(stiffness)):
        raise ValueError("K must be positive semidefinite")
    return mass, stiffness


def normal_modes(M, K) -> ModalResult:
    r"""Natural frequencies and mode shapes from the eigenproblem :math:`K\phi = \omega^2 M\phi`.

    Parameters
    ----------
    M : array_like
        Mass matrix [kg or kg m^2], symmetric positive definite, ``(n, n)``.
    K : array_like
        Stiffness matrix [N/m or N m/rad], symmetric positive semidefinite,
        ``(n, n)``.

    Returns
    -------
    ModalResult
        All ``n`` frequencies [rad/s] in increasing order, zero for rigid-body
        modes, and mode shapes normalised so :math:`\Phi^T M \Phi = I` and
        :math:`\Phi^T K \Phi = \operatorname{diag}(\omega^2)`.

    Notes
    -----
    Damping is ignored; the frequencies are undamped natural frequencies.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    J. L. Lagrange, *Mécanique analytique* (Desaint, Paris, 1788), part II,
    section VI; D. Bernoulli, "Réflexions et éclaircissemens sur les
    nouvelles vibrations des cordes," *Hist. Acad. R. Sci. Belles Lettres
    Berlin* 9 (1753, published 1755) 147-172.

    Examples
    --------
    >>> modes = normal_modes(np.eye(2), [[2.0, -1.0], [-1.0, 2.0]])
    >>> np.round(modes.frequencies ** 2, 12).tolist()
    [1.0, 3.0]
    """
    mass, stiffness = _system(M, K)
    w2, phi = eigh(stiffness, mass)
    shapes = phi.T
    shapes = shapes * np.where(shapes[np.arange(len(shapes)), np.argmax(np.abs(shapes), axis=1)] < 0, -1.0, 1.0)[:, None]
    return ModalResult(frequencies=np.sqrt(np.clip(w2, 0, None)), mode_shapes=shapes)


def modal_response(M, K, x0, v0=None, t=None, modal_damping=0.0) -> MDOFResult:
    r"""Free vibration by modal superposition, :math:`x(t) = \sum_r \phi_r q_r(t)`.

    Each modal coordinate :math:`q_r = \phi_r^T M x` starts from the
    projection of the initial state and decays as a damped single-degree-of-
    freedom oscillator, :math:`\ddot q_r + 2\zeta_r\omega_r\dot q_r +
    \omega_r^2 q_r = 0`, in closed form.

    Parameters
    ----------
    M, K : array_like
        Mass and stiffness matrices, as for :func:`normal_modes`.
    x0 : array_like
        Initial displacements [m or rad], length ``n``.
    v0 : array_like, optional
        Initial velocities [m/s or rad/s]; zero if omitted.
    t : array_like, optional
        Output times [s]; 1001 samples on [0, 10] s if omitted.
    modal_damping : float or array_like, optional
        Damping ratio of each mode [-] in [0, 1), one value for all modes or
        one per mode. Rigid-body modes are left undamped.

    Returns
    -------
    MDOFResult
        Displacements and velocities at each time, shape ``(len(t), n)``.

    Notes
    -----
    Modal damping keeps the modes uncoupled; a general damping matrix does
    not. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> result = modal_response(np.eye(2), [[2.0, -1.0], [-1.0, 2.0]], [1.0, 1.0], t=[0.0, np.pi])
    >>> np.round(result.displacement[-1], 12).tolist()
    [-1.0, -1.0]
    """
    mass, stiffness = _system(M, K)
    n = mass.shape[0]
    x = _vector("x0", x0)
    v = np.zeros(n) if v0 is None else _vector("v0", v0)
    if x.size != n or v.size != n:
        raise ValueError("x0 and v0 must have one entry per degree of freedom")
    times = np.linspace(0, 10, 1001) if t is None else _vector("t", t)
    zeta = (
        np.broadcast_to(np.asarray(modal_damping, dtype=float), (n,))
        if np.ndim(modal_damping) == 0
        else _vector("modal_damping", modal_damping)
    )
    if zeta.size != n:
        raise ValueError("modal_damping must be a scalar or have one entry per mode")
    _interval("modal_damping", zeta, 0, 1, closed="left")
    modes = normal_modes(mass, stiffness)
    phi, w = modes.mode_shapes, modes.frequencies
    q0, p0 = phi @ mass @ x, phi @ mass @ v
    tt = times[:, None]
    rigid = w < 1e-9 * max(1.0, float(np.max(w)))
    z = np.where(rigid, 0.0, zeta)
    wd = w * np.sqrt(1 - z * z)
    safe = np.where(rigid, 1.0, wd)
    decay = np.exp(-z * w * tt)
    b = (p0 + z * w * q0) / safe
    q = np.where(rigid, q0 + p0 * tt, decay * (q0 * np.cos(wd * tt) + b * np.sin(wd * tt)))
    qdot = np.where(rigid, p0 + 0 * tt, decay * ((b * wd - z * w * q0) * np.cos(wd * tt) - (q0 * wd + z * w * b) * np.sin(wd * tt)))
    return MDOFResult(t=times, displacement=q @ phi, velocity=qdot @ phi)


def harmonic_response(M, K, F, omega, C=None):
    r"""Steady-state complex amplitudes :math:`X = (K - \omega^2 M + i\omega C)^{-1} F`.

    Parameters
    ----------
    M, K : array_like
        Mass and stiffness matrices, as for :func:`normal_modes`.
    F : array_like
        Force amplitudes on each coordinate [N or N m], length ``n``.
    omega : float or array_like
        Forcing frequencies [rad/s], nonnegative.
    C : array_like, optional
        Viscous damping matrix [N s/m], symmetric; undamped if omitted.

    Returns
    -------
    ndarray
        Complex displacement amplitudes [m or rad], shape ``(n,)`` for a
        scalar ``omega`` or ``(len(omega), n)``; the response is
        :math:`\operatorname{Re}(X e^{i\omega t})`.

    Raises
    ------
    ValueError
        If the dynamic stiffness is singular, as at an undamped resonance.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> X = harmonic_response(np.eye(2), [[2.0, -1.0], [-1.0, 2.0]], [1.0, 0.0], 0.0)
    >>> np.round(X.real, 12).tolist()
    [0.666666666667, 0.333333333333]
    """
    mass, stiffness = _system(M, K)
    n = mass.shape[0]
    damping = np.zeros_like(mass) if C is None else _square("C", C, n)
    force = _vector("F", F)
    if force.size != n:
        raise ValueError("F must have one entry per degree of freedom")
    _nonnegative(omega=omega)
    w = np.asarray(omega, dtype=float)
    dynamic = stiffness - (w.reshape(-1, 1, 1) ** 2) * mass + 1j * w.reshape(-1, 1, 1) * damping
    if np.any(np.linalg.cond(dynamic) > 1e13):
        raise ValueError("the dynamic stiffness is singular at some omega (undamped resonance)")
    X = np.linalg.solve(dynamic, np.broadcast_to(force.astype(complex), (w.size, n))[..., None])[..., 0]
    return X[0] if w.ndim == 0 else X


def rayleigh_frequency(M, K, shape):
    r"""Rayleigh's quotient :math:`\omega^2 \approx x^T K x / x^T M x` for a trial shape.

    Parameters
    ----------
    M, K : array_like
        Mass and stiffness matrices, as for :func:`normal_modes`.
    shape : array_like
        Assumed mode shape [-], length ``n``, not zero.

    Returns
    -------
    float
        Estimated frequency [rad/s]. It is never below the lowest natural
        frequency, and its error is second order in the error of the shape.

    Notes
    -----
    The static deflection under the weights is a good trial shape for the
    fundamental mode of a shaft or beam. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/mechanical_engineering`.

    References
    ----------
    J. W. Strutt, Baron Rayleigh, "Some general theorems relating to
    vibrations," *Proc. Lond. Math. Soc.* 4 (1873) 357-368; *The Theory of
    Sound*, vol. 1 (Macmillan, London, 1877), sections 88-89.

    Examples
    --------
    >>> round(rayleigh_frequency(np.eye(2), [[2.0, -1.0], [-1.0, 2.0]], [1.0, 1.0]), 12)
    1.0
    """
    mass, stiffness = _system(M, K)
    x = _vector("shape", shape)
    if x.size != mass.shape[0] or not np.any(x):
        raise ValueError("shape must be a nonzero vector with one entry per degree of freedom")
    return float(np.sqrt((x @ stiffness @ x) / (x @ mass @ x)))


def dunkerley_frequency(masses, flexibility):
    r"""Dunkerley's lower bound :math:`1/\omega^2 \approx \sum_i m_i a_{ii}` on the fundamental frequency.

    Each term :math:`m_i a_{ii} = 1/\omega_i^2` is the inverse square of the
    frequency the shaft would have with mass :math:`i` alone on it.

    Parameters
    ----------
    masses : array_like
        Lumped masses [kg], positive.
    flexibility : array_like
        Influence coefficients :math:`a_{ij}` [m/N]: deflection at station
        :math:`i` per unit force at :math:`j`, symmetric ``(n, n)``.

    Returns
    -------
    float
        Estimated fundamental frequency [rad/s], never above the true one.

    Notes
    -----
    The sum is the trace of the dynamic matrix :math:`A M`, the sum of every
    :math:`1/\omega_r^2`, so the estimate is close when the higher modes are
    well separated from the first. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/mechanical_engineering`.

    References
    ----------
    S. Dunkerley, "On the whirling and vibration of shafts," *Phil. Trans.
    R. Soc. Lond. A* 185 (1894) 279-360.

    Examples
    --------
    >>> round(dunkerley_frequency([1.0, 1.0], [[2 / 3, 1 / 3], [1 / 3, 2 / 3]]), 12)
    0.866025403784
    """
    m = _vector("masses", masses)
    _positive(masses=m)
    a = _square("flexibility", flexibility, m.size)
    if np.any(np.diag(a) <= 0):
        raise ValueError("the diagonal influence coefficients must be positive")
    return float(1 / np.sqrt(np.sum(m * np.diag(a))))
