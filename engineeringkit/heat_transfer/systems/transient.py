r"""Transient conduction: lumped bodies, the semi-infinite solid and the plane wall.

Biot (1804) saw that a heated bar's temperature profile depends on the
competition between internal conduction and surface loss, the ratio now
called the Biot number. Kelvin (1862) used the error-function solution for a
semi-infinite solid to estimate the age of the Earth, and Heisler (1947)
charted the series solutions for slabs, long cylinders and spheres.
"""

from __future__ import annotations

import numpy as np
from scipy import special

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.heat_transfer.utils.eigenvalues import cylinder_eigenvalues, plane_wall_eigenvalues, sphere_eigenvalues

__all__ = [
    "biot_number",
    "fourier_number",
    "lumped_temperature",
    "semi_infinite_temperature",
    "plane_wall_transient",
    "cylinder_transient",
    "sphere_transient",
]


def biot_number(h, Lc, k):
    r"""Biot number :math:`Bi = hL_c/k`, internal conduction resistance over surface convection resistance.

    Parameters
    ----------
    h : float or array_like
        Surface convection coefficient [W/(m^2 K)], nonnegative.
    Lc : float or array_like
        Characteristic length [m], positive: :math:`V/A` for a lumped body,
        the half-thickness for a plane wall.
    k : float or array_like
        Conductivity of the *solid* [W/(m K)], positive.

    Returns
    -------
    float or ndarray
        Biot number [-]; below about 0.1 the body is nearly isothermal.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    References
    ----------
    J.-B. Biot, "Mémoire sur la propagation de la chaleur," *Bibliothèque
    britannique* 27 (1804) 310-329.

    Examples
    --------
    >>> biot_number(10.0, 0.01, 1.0)
    0.1
    """
    _finite(h=h, Lc=Lc, k=k)
    _nonnegative(h=h)
    _positive(Lc=Lc, k=k)
    return h * Lc / k


def fourier_number(alpha, t, Lc):
    r"""Fourier number :math:`Fo = \alpha t/L_c^2`, dimensionless time of transient conduction.

    Parameters
    ----------
    alpha : float or array_like
        Thermal diffusivity [m^2/s], positive.
    t : float or array_like
        Time [s], nonnegative.
    Lc : float or array_like
        Characteristic length [m], positive.

    Returns
    -------
    float or ndarray
        Fourier number [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    Examples
    --------
    >>> fourier_number(1e-5, 100.0, 0.01)
    10.0
    """
    _finite(alpha=alpha, t=t, Lc=Lc)
    _positive(alpha=alpha, Lc=Lc)
    _nonnegative(t=t)
    return alpha * t / Lc**2


def lumped_temperature(Ti, Tinf, h, A, rho, cp, V, t):
    r"""Temperature of a lumped (isothermal) body cooling by convection.

    .. math::

        \frac{T - T_\infty}{T_i - T_\infty} = \exp\left(-\frac{hA}{\rho c_p V}\,t\right)

    Parameters
    ----------
    Ti : float
        Initial temperature [K].
    Tinf : float
        Fluid temperature [K].
    h : float
        Convection coefficient [W/(m^2 K)], nonnegative.
    A : float
        Surface area [m^2], positive.
    rho : float
        Density [kg/m^3], positive.
    cp : float
        Specific heat [J/(kg K)], positive.
    V : float
        Volume [m^3], positive.
    t : float or array_like
        Time [s], nonnegative.

    Returns
    -------
    float or ndarray
        Body temperature [K].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Valid when :func:`biot_number` with
    :math:`L_c = V/A` is below about 0.1; that is not checked here.

    Examples
    --------
    >>> round(float(lumped_temperature(400.0, 300.0, 10.0, 1.0, 1000.0, 1000.0, 0.01, 1000.0)), 2)
    336.79
    """
    _finite(Ti=Ti, Tinf=Tinf, h=h, A=A, rho=rho, cp=cp, V=V, t=t)
    _nonnegative(Ti=Ti, Tinf=Tinf, h=h, t=t)
    _positive(A=A, rho=rho, cp=cp, V=V)
    return Tinf + (Ti - Tinf) * np.exp(-h * A * np.asarray(t) / (rho * cp * V))


def semi_infinite_temperature(Ti, Ts, x, alpha, t):
    r"""Semi-infinite solid after a step change of surface temperature.

    .. math::

        \frac{T(x,t) - T_s}{T_i - T_s} = \operatorname{erf}\frac{x}{2\sqrt{\alpha t}}

    Parameters
    ----------
    Ti : float
        Initial uniform temperature [K].
    Ts : float
        Surface temperature imposed at :math:`t = 0` [K].
    x : float or array_like
        Depth below the surface [m], nonnegative.
    alpha : float
        Thermal diffusivity [m^2/s], positive.
    t : float or array_like
        Time since the step [s], positive.

    Returns
    -------
    float or ndarray
        Temperature [K]. The surface gradient is
        :math:`(T_i - T_s)/\sqrt{\pi\alpha t}`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`.

    References
    ----------
    W. Thomson, "On the secular cooling of the earth," *Trans. R. Soc.
    Edinburgh* 23 (1864) 157-169 (read 1862).

    Examples
    --------
    >>> round(float(semi_infinite_temperature(300.0, 400.0, 0.01, 1e-5, 2.5)), 2)
    315.73
    """
    _finite(Ti=Ti, Ts=Ts, x=x, alpha=alpha, t=t)
    _nonnegative(Ti=Ti, Ts=Ts, x=x)
    _positive(alpha=alpha, t=t)
    return Ts + (Ti - Ts) * special.erf(np.asarray(x) / (2 * np.sqrt(alpha * np.asarray(t))))


def plane_wall_transient(Fo, Bi, x_star=0.0, terms=100):
    r"""Dimensionless temperature of a plane wall cooled by convection on both faces (Heisler's charts).

    .. math::

        \theta^* = \frac{T - T_\infty}{T_i - T_\infty}
                 = \sum_{n=1}^{N} C_n e^{-\zeta_n^2 Fo}\cos(\zeta_n x^*),
        \qquad C_n = \frac{4\sin\zeta_n}{2\zeta_n + \sin 2\zeta_n},

    with :math:`\zeta_n\tan\zeta_n = Bi`. ``terms=1`` gives the one-term
    approximation behind Heisler's charts, within 2% at the midplane
    for :math:`Fo > 0.2` (about 3% off-centre).

    Parameters
    ----------
    Fo : float or array_like
        Fourier number :math:`\alpha t/L^2` on the half-thickness :math:`L` [-], nonnegative.
    Bi : float
        Biot number :math:`hL/k` on the half-thickness [-], positive.
    x_star : float or array_like, optional
        Position :math:`x/L` from the midplane [-], in [0, 1]; 0 is the centre.
    terms : int, optional
        Number of series terms :math:`N`, at least one.

    Returns
    -------
    float or ndarray
        :math:`\theta^*` [-], broadcast over ``Fo`` and ``x_star``.

    Notes
    -----
    ``Bi`` is scalar-only, since the eigenvalues are roots found for one Biot
    number; ``Fo`` and ``x_star`` broadcast.

    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Near :math:`Fo = 0` the truncated series
    oscillates about the initial value; use more terms there.

    References
    ----------
    M. P. Heisler, "Temperature charts for induction and constant-temperature
    heating," *Trans. ASME* 69 (1947) 227-236.

    Examples
    --------
    >>> round(float(plane_wall_transient(0.5, 1.0)), 4)
    0.7725
    """
    _nonnegative(Fo=Fo)
    _interval("x_star", x_star, 0, 1)
    zeta = plane_wall_eigenvalues(Bi, terms)
    C = 4 * np.sin(zeta) / (2 * zeta + np.sin(2 * zeta))
    Fo, x_star = np.broadcast_arrays(np.asarray(Fo, dtype=float), np.asarray(x_star, dtype=float))
    modes = C * np.exp(-(zeta**2) * Fo[..., None]) * np.cos(zeta * x_star[..., None])
    return modes.sum(axis=-1)[()]


def cylinder_transient(Fo, Bi, r_star=0.0, terms=100):
    r"""Dimensionless temperature of a long cylinder cooled by convection (Heisler's charts).

    .. math::

        \theta^* = \sum_{n=1}^{N} C_n e^{-\zeta_n^2 Fo} J_0(\zeta_n r^*),
        \qquad C_n = \frac{2}{\zeta_n}\frac{J_1(\zeta_n)}{J_0^2(\zeta_n) + J_1^2(\zeta_n)},

    with :math:`\zeta_n J_1(\zeta_n)/J_0(\zeta_n) = Bi`.

    Parameters
    ----------
    Fo : float or array_like
        Fourier number :math:`\alpha t/r_o^2` on the outer radius [-], nonnegative.
    Bi : float
        Biot number :math:`hr_o/k` on the outer radius [-], positive.
    r_star : float or array_like, optional
        Radial position :math:`r/r_o` [-], in [0, 1]; 0 is the axis.
    terms : int, optional
        Number of series terms :math:`N`, at least one; ``terms=1`` is
        Heisler's one-term chart, good for :math:`Fo > 0.2`.

    Returns
    -------
    float or ndarray
        :math:`\theta^* = (T - T_\infty)/(T_i - T_\infty)` [-], broadcast over ``Fo`` and ``r_star``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. The cylinder is infinitely long (no end
    losses); near :math:`Fo = 0` use more terms.

    References
    ----------
    M. P. Heisler, "Temperature charts for induction and constant-temperature
    heating," *Trans. ASME* 69 (1947) 227-236.

    Examples
    --------
    >>> round(float(cylinder_transient(0.5, 1.0)), 4)
    0.5486
    """
    _nonnegative(Fo=Fo)
    _interval("r_star", r_star, 0, 1)
    zeta = cylinder_eigenvalues(Bi, terms)
    j0, j1 = special.j0(zeta), special.j1(zeta)
    C = 2 / zeta * j1 / (j0**2 + j1**2)
    Fo, r_star = np.broadcast_arrays(np.asarray(Fo, dtype=float), np.asarray(r_star, dtype=float))
    modes = C * np.exp(-(zeta**2) * Fo[..., None]) * special.j0(zeta * r_star[..., None])
    return modes.sum(axis=-1)[()]


def sphere_transient(Fo, Bi, r_star=0.0, terms=100):
    r"""Dimensionless temperature of a sphere cooled by convection (Heisler's charts).

    .. math::

        \theta^* = \sum_{n=1}^{N} C_n e^{-\zeta_n^2 Fo}\,\frac{\sin(\zeta_n r^*)}{\zeta_n r^*},
        \qquad C_n = \frac{4(\sin\zeta_n - \zeta_n\cos\zeta_n)}{2\zeta_n - \sin 2\zeta_n},

    with :math:`1 - \zeta_n\cot\zeta_n = Bi`; at the centre
    :math:`\sin(\zeta r^*)/(\zeta r^*) = 1`.

    Parameters
    ----------
    Fo : float or array_like
        Fourier number :math:`\alpha t/r_o^2` on the outer radius [-], nonnegative.
    Bi : float
        Biot number :math:`hr_o/k` on the outer radius [-], positive.
    r_star : float or array_like, optional
        Radial position :math:`r/r_o` [-], in [0, 1]; 0 is the centre.
    terms : int, optional
        Number of series terms :math:`N`, at least one; ``terms=1`` is
        Heisler's one-term chart, good for :math:`Fo > 0.2`.

    Returns
    -------
    float or ndarray
        :math:`\theta^* = (T - T_\infty)/(T_i - T_\infty)` [-], broadcast over ``Fo`` and ``r_star``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/heat_transfer`. Near :math:`Fo = 0` use more terms.

    References
    ----------
    M. P. Heisler, "Temperature charts for induction and constant-temperature
    heating," *Trans. ASME* 69 (1947) 227-236.

    Examples
    --------
    >>> round(float(sphere_transient(0.5, 1.0)), 4)
    0.3708
    """
    _nonnegative(Fo=Fo)
    _interval("r_star", r_star, 0, 1)
    zeta = sphere_eigenvalues(Bi, terms)
    C = 4 * (np.sin(zeta) - zeta * np.cos(zeta)) / (2 * zeta - np.sin(2 * zeta))
    Fo, r_star = np.broadcast_arrays(np.asarray(Fo, dtype=float), np.asarray(r_star, dtype=float))
    # np.sinc(x) = sin(pi x)/(pi x) is 1 at x = 0, so the centre needs no special case.
    modes = C * np.exp(-(zeta**2) * Fo[..., None]) * np.sinc(zeta * r_star[..., None] / np.pi)
    return modes.sum(axis=-1)[()]
