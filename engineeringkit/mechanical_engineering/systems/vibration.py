r"""Vibration: the damped oscillator, the tuned absorber and torsional chains.

Rayleigh's *Theory of Sound* set out the free and forced vibration of linear
systems with natural frequency :math:`\omega_n = \sqrt{k/m}` and damping
ratio :math:`\zeta = c/2\sqrt{km}`. Frahm hung a small tuned oscillator on a
vibrating body to cancel its motion, and Den Hartog chose the absorber's
tuning and damping to flatten the response. Holzer found the torsional
frequencies of a shaft line by marching a trial twist from disk to disk.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.integrators import integrate_ode
from engineeringkit.mechanical_engineering.core.base import HolzerResult, SDOFResult, TunedMassDamper
from engineeringkit.mechanical_engineering.utils.roots import bracketed_roots

__all__ = [
    "natural_frequency",
    "natural_frequency_hz",
    "damping_ratio",
    "damped_frequency",
    "transmissibility",
    "simulate_sdof",
    "tuned_mass_damper",
    "holzer_residual",
    "holzer_frequencies",
]


def natural_frequency(k, m):
    r"""Undamped natural frequency :math:`\omega_n = \sqrt{k/m}`.

    Parameters
    ----------
    k : float or array_like
        Stiffness [N/m], positive.
    m : float or array_like
        Mass [kg], positive.

    Returns
    -------
    float or ndarray
        Natural frequency [rad/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    J. W. Strutt, Baron Rayleigh, *The Theory of Sound*, vol. 1 (Macmillan,
    London, 1877).

    Examples
    --------
    >>> float(natural_frequency(100.0, 4.0))
    5.0
    """
    _finite(k=k, m=m)
    _positive(k=k, m=m)
    return np.sqrt(k / m)


def natural_frequency_hz(k, m):
    r"""Undamped natural frequency in hertz, :math:`f_n = \sqrt{k/m}/2\pi`.

    Parameters
    ----------
    k : float or array_like
        Stiffness [N/m], positive.
    m : float or array_like
        Mass [kg], positive.

    Returns
    -------
    float or ndarray
        Natural frequency [Hz].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(float(natural_frequency_hz(4.0, 1.0)), 6)
    0.31831
    """
    return natural_frequency(k, m) / (2 * np.pi)


def damping_ratio(c, k, m):
    r"""Damping ratio :math:`\zeta = c / (2\sqrt{km})`, the fraction of critical damping.

    Parameters
    ----------
    c : float or array_like
        Viscous damping coefficient [N s/m], nonnegative.
    k : float or array_like
        Stiffness [N/m], positive.
    m : float or array_like
        Mass [kg], positive.

    Returns
    -------
    float or ndarray
        Damping ratio [-]; one is critical damping.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> float(damping_ratio(4.0, 4.0, 1.0))
    1.0
    """
    _finite(c=c, k=k, m=m)
    _nonnegative(c=c)
    _positive(k=k, m=m)
    return c / (2 * np.sqrt(k * m))


def damped_frequency(k, m, c):
    r"""Damped natural frequency :math:`\omega_d = \omega_n\sqrt{1 - \zeta^2}`.

    Parameters
    ----------
    k : float or array_like
        Stiffness [N/m], positive.
    m : float or array_like
        Mass [kg], positive.
    c : float or array_like
        Viscous damping coefficient [N s/m], nonnegative.

    Returns
    -------
    float or ndarray
        Frequency of the decaying oscillation [rad/s]; zero at or above
        critical damping, where the motion does not oscillate.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(float(damped_frequency(4.0, 1.0, 0.4)), 6)
    1.989975
    """
    wn = natural_frequency(k, m)
    z = damping_ratio(c, k, m)
    return wn * np.sqrt(np.maximum(0, 1 - z * z))


def transmissibility(r, zeta):
    r"""Force (or base-motion) transmissibility of a damped linear isolator.

    .. math:: T = \sqrt{\frac{1 + (2\zeta r)^2}{(1 - r^2)^2 + (2\zeta r)^2}}

    Parameters
    ----------
    r : float or array_like
        Excitation frequency over natural frequency [-], nonnegative.
    zeta : float or array_like
        Damping ratio [-], nonnegative.

    Returns
    -------
    float or ndarray
        Transmitted over applied amplitude [-]; below one only for
        :math:`r > \sqrt 2`, and infinite at undamped resonance.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(float(transmissibility(np.sqrt(2), 0.1)), 12)
    1.0
    """
    _nonnegative(r=r, zeta=zeta)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.sqrt(np.divide(1 + (2 * zeta * r) ** 2, (1 - r * r) ** 2 + (2 * zeta * r) ** 2))


def simulate_sdof(m, c, k, x0=0, v0=0, t=None, force=None):
    r"""Integrate :math:`m\ddot x + c\dot x + kx = F(t)` for a passive linear oscillator.

    Parameters
    ----------
    m : float
        Mass [kg], positive.
    c : float
        Viscous damping [N s/m], nonnegative.
    k : float
        Stiffness [N/m], nonnegative.
    x0 : float, optional
        Initial displacement [m].
    v0 : float, optional
        Initial velocity [m/s].
    t : array_like, optional
        Strictly increasing output times [s]; defaults to 1001 samples on
        [0, 10] s.
    force : callable, optional
        ``force(time [s]) -> force [N]``; omit for free vibration.

    Returns
    -------
    SDOFResult
        Time [s], displacement [m] and velocity [m/s].

    Notes
    -----
    Constant coefficients and small linear displacements are assumed.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> result = simulate_sdof(1.0, 0.0, 4.0, x0=0.1, t=np.linspace(0, np.pi, 3))
    >>> np.round(result.displacement, 6).tolist()
    [0.1, -0.1, 0.1]
    """
    if not np.all(np.isfinite([m, c, k, x0, v0])) or m <= 0 or c < 0 or k < 0:
        raise ValueError("mass must be positive; damping/stiffness nonnegative; inputs finite")
    t = np.linspace(0, 10, 1001) if t is None else np.asarray(t)

    def rhs(tt, y):
        f = 0 if force is None else force(tt)
        return [y[1], (f - c * y[1] - k * y[0]) / m]

    result = integrate_ode(rhs, [x0, v0], t)
    return SDOFResult(result.t, result.history[:, 0], result.history[:, 1])


def tuned_mass_damper(m, k, mass_ratio, frequency_ratio=None, damping_ratio=None) -> TunedMassDamper:
    r"""Design a dynamic vibration absorber for a mass-spring system.

    With ``frequency_ratio`` and ``damping_ratio`` omitted, Den Hartog's
    optimum tuning is used:

    .. math:: f = \frac{1}{1 + \mu}, \qquad
              \zeta = \sqrt{\frac{3\mu}{8(1 + \mu)^3}}

    which makes the two fixed points of the response curve equal in height,
    :math:`\sqrt{1 + 2/\mu}` times the static deflection, and puts the
    peaks near them.

    Parameters
    ----------
    m : float
        Primary mass :math:`m_1` [kg], positive.
    k : float
        Primary stiffness :math:`k_1` [N/m], positive.
    mass_ratio : float
        :math:`\mu = m_2/m_1` [-], positive.
    frequency_ratio : float, optional
        Absorber frequency :math:`\sqrt{k_2/m_2}` over :math:`\omega_1` [-],
        positive. Frahm's undamped absorber uses 1, tuned to the forcing.
    damping_ratio : float, optional
        :math:`c_2/(2 m_2 \omega_1)` [-], nonnegative.

    Returns
    -------
    TunedMassDamper
        Absorber mass, stiffness and damping, with ``amplification(r)`` for
        the primary mass's frequency response.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    H. Frahm, "Device for damping vibrations of bodies," U.S. Patent 989,958
    (filed 1909, issued 1911); J. Ormondroyd and J. P. Den Hartog, "The
    theory of the dynamic vibration absorber," *Trans. ASME* 50 (1928)
    APM-50-7, 9-22.

    Examples
    --------
    >>> tmd = tuned_mass_damper(100.0, 1e5, 0.05)
    >>> round(tmd.frequency_ratio, 6), round(tmd.damping_ratio, 6)
    (0.952381, 0.127267)
    """
    _positive(m=m, k=k, mass_ratio=mass_ratio)
    f = 1 / (1 + mass_ratio) if frequency_ratio is None else frequency_ratio
    z = np.sqrt(3 * mass_ratio / (8 * (1 + mass_ratio) ** 3)) if damping_ratio is None else damping_ratio
    _positive(frequency_ratio=f)
    _nonnegative(damping_ratio=z)
    w1 = np.sqrt(k / m)
    m2 = mass_ratio * m
    return TunedMassDamper(
        primary_mass=float(m),
        primary_stiffness=float(k),
        absorber_mass=float(m2),
        absorber_stiffness=float(m2 * (f * w1) ** 2),
        absorber_damping=float(2 * m2 * w1 * z),
        mass_ratio=float(mass_ratio),
        frequency_ratio=float(f),
        damping_ratio=float(z),
    )


def _chain(inertias, stiffnesses):
    J = _vector("inertias", inertias, min_length=2)
    K = _vector("stiffnesses", stiffnesses)
    if K.size != J.size - 1:
        raise ValueError("one shaft stiffness is required between each pair of disks")
    _positive(inertias=J, stiffnesses=K)
    return J, K


def _holzer_march(omega, J, K):
    w2 = np.asarray(omega, dtype=float)[..., None] ** 2
    twist = np.ones(w2.shape[:-1] + (J.size,))
    torque = np.zeros(w2.shape[:-1])
    for i in range(J.size):
        if i:
            twist[..., i] = twist[..., i - 1] - torque / K[i - 1]
        torque = torque + J[i] * w2[..., 0] * twist[..., i]
    return twist, torque


def holzer_residual(omega, inertias, stiffnesses):
    r"""Holzer's residual torque beyond the last disk of a free-free shaft.

    Start with unit twist at disk 1 and no torque; at each disk add its
    inertia torque :math:`J_i\omega^2\theta_i`, then twist the next shaft
    by the accumulated torque over its stiffness. ``omega`` is a natural
    frequency when the torque left over at the free end is zero.

    Parameters
    ----------
    omega : float or array_like
        Trial frequency [rad/s].
    inertias : array_like
        Disk polar moments of inertia [kg m^2], positive, at least two.
    stiffnesses : array_like
        Torsional stiffness of each shaft between disks [N m/rad], positive,
        one fewer than ``inertias``.

    Returns
    -------
    float or ndarray
        Residual torque per unit twist of disk 1 [N m/rad].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(abs(holzer_residual(np.sqrt(2.0), [1.0, 1.0], [1.0])), 12)
    0.0
    """
    _finite(omega=omega)
    J, K = _chain(inertias, stiffnesses)
    torque = _holzer_march(omega, J, K)[1]
    return float(torque) if np.ndim(omega) == 0 else torque


def holzer_frequencies(inertias, stiffnesses, omega_max, n_grid=4000) -> HolzerResult:
    r"""Torsional natural frequencies of a free-free shaft by Holzer's method.

    Parameters
    ----------
    inertias : array_like
        Disk polar moments of inertia [kg m^2], positive, at least two.
    stiffnesses : array_like
        Shaft torsional stiffnesses [N m/rad], positive, one fewer than
        ``inertias``.
    omega_max : float
        Upper end of the frequency search [rad/s], positive.
    n_grid : int, optional
        Trial frequencies scanned for sign changes of the residual torque,
        at least 10; closely spaced modes need a finer grid.

    Returns
    -------
    HolzerResult
        Nonzero natural frequencies below ``omega_max`` [rad/s] and each
        mode's disk twists, scaled so the largest is one. The zero-frequency
        rigid-body rotation is excluded.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    H. Holzer, *Die Berechnung der Drehschwingungen* (Springer, Berlin,
    1921).

    Examples
    --------
    >>> result = holzer_frequencies([1.0, 1.0], [1.0], 5.0)
    >>> np.round(result.frequencies ** 2, 9).tolist()
    [2.0]
    """
    J, K = _chain(inertias, stiffnesses)
    _positive(omega_max=omega_max)
    _integer("n_grid", n_grid, 10)
    grid = np.linspace(0, omega_max, n_grid + 1)[1:]
    frequencies = bracketed_roots(lambda w: _holzer_march(w, J, K)[1], grid)
    twist = _holzer_march(frequencies, J, K)[0].reshape(frequencies.size, J.size)
    shapes = twist / np.max(np.abs(twist), axis=1, keepdims=True) if frequencies.size else twist
    return HolzerResult(frequencies=frequencies, mode_shapes=shapes)
