r"""Stress and strain at a point: Mohr's circle, three-dimensional principal stresses and strain-gauge rosettes.

Rotating a plane-stress element by :math:`\theta` moves its stress point
around a circle of centre :math:`(\sigma_x + \sigma_y)/2` and radius
:math:`\sqrt{((\sigma_x - \sigma_y)/2)^2 + \tau_{xy}^2}` through :math:`2\theta`.
In three dimensions the principal stresses are the eigenvalues of the stress
tensor, and the absolute maximum shear stress is :math:`(\sigma_1 - \sigma_3)/2`,
which for plane stress may lie out of the plane. Strains transform the same
way with :math:`\gamma_{xy}/2` in place of :math:`\tau_{xy}`, so three gauges at
known angles determine the plane strain state.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector

__all__ = ["mohr_circle", "stress_transformation", "principal_stresses_3d", "strain_rosette", "plane_stress_from_strains"]


def mohr_circle(sx, sy, txy):
    """Centre, radius and principal stresses of Mohr's circle for plane stress.

    Parameters
    ----------
    sx, sy : float or array_like
        Normal stresses on the x and y faces [Pa], tension positive.
    txy : float or array_like
        Shear stress [Pa].

    Returns
    -------
    tuple of float or ndarray
        ``(center, radius, sigma_1, sigma_2)``, all [Pa]: the circle's
        centre and radius, and the maximum and minimum principal stresses.
        The radius is also the maximum in-plane shear stress.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.

    References
    ----------
    O. Mohr, "Über die Darstellung des Spannungszustandes und des
    Deformationszustandes eines Körperelementes," *Zivilingenieur* 28
    (1882) 113-156.
    """
    _finite(sx=sx, sy=sy, txy=txy)
    c = (sx + sy) / 2
    r = np.hypot((sx - sy) / 2, txy)
    return c, r, c + r, c - r


def stress_transformation(sx, sy, txy, theta):
    r"""Stresses on an element rotated counterclockwise by ``theta``.

    .. math::

        \sigma_{x'} = \frac{\sigma_x + \sigma_y}{2} + \frac{\sigma_x - \sigma_y}{2}\cos 2\theta + \tau_{xy}\sin 2\theta,\\
        \tau_{x'y'} = -\frac{\sigma_x - \sigma_y}{2}\sin 2\theta + \tau_{xy}\cos 2\theta,

    and :math:`\sigma_{y'} = \sigma_x + \sigma_y - \sigma_{x'}`.

    Parameters
    ----------
    sx, sy : float or array_like
        Normal stresses [Pa], tension positive.
    txy : float or array_like
        Shear stress [Pa].
    theta : float or array_like
        Rotation of the element [rad], counterclockwise positive.

    Returns
    -------
    tuple of float or ndarray
        ``(sigma_x', sigma_y', tau_x'y')`` [Pa].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(sx=sx, sy=sy, txy=txy, theta=theta)
    mean, half = (np.asarray(sx) + sy) / 2, (np.asarray(sx) - sy) / 2
    c2, s2 = np.cos(2 * np.asarray(theta)), np.sin(2 * np.asarray(theta))
    sxp = mean + half * c2 + txy * s2
    return sxp, 2 * mean - sxp, -half * s2 + txy * c2


def principal_stresses_3d(sx, sy, sz, txy=0.0, tyz=0.0, tzx=0.0):
    r"""Principal stresses of a three-dimensional stress state.

    They are the eigenvalues of the symmetric stress tensor

    .. math::

        \boldsymbol\sigma = \begin{pmatrix}\sigma_x & \tau_{xy} & \tau_{zx}\\
        \tau_{xy} & \sigma_y & \tau_{yz}\\ \tau_{zx} & \tau_{yz} & \sigma_z\end{pmatrix},

    the roots of :math:`\sigma^3 - I_1\sigma^2 + I_2\sigma - I_3 = 0`. Mohr
    drew the three circles through each pair of them; the largest has
    radius :math:`\tau_\mathrm{abs} = (\sigma_1 - \sigma_3)/2`.

    Parameters
    ----------
    sx, sy, sz : float or array_like
        Normal stresses [Pa], tension positive.
    txy, tyz, tzx : float or array_like, optional
        Shear stresses [Pa].

    Returns
    -------
    tuple of float or ndarray
        ``(sigma_1, sigma_2, sigma_3)`` [Pa], in decreasing order.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. ``numpy.linalg.eigvalsh`` computes the eigenvalues.

    References
    ----------
    O. Mohr, "Über die Darstellung des Spannungszustandes und des
    Deformationszustandes eines Körperelementes," *Zivilingenieur* 28
    (1882) 113-156.
    """
    _finite(sx=sx, sy=sy, sz=sz, txy=txy, tyz=tyz, tzx=tzx)
    sx, sy, sz, txy, tyz, tzx = np.broadcast_arrays(*(np.asarray(v, float) for v in (sx, sy, sz, txy, tyz, tzx)))
    tensor = np.stack([np.stack([sx, txy, tzx], -1), np.stack([txy, sy, tyz], -1), np.stack([tzx, tyz, sz], -1)], -2)
    values = np.linalg.eigvalsh(tensor)[..., ::-1]
    if values.ndim == 1:
        return tuple(float(v) for v in values)
    return values[..., 0], values[..., 1], values[..., 2]


def strain_rosette(strains, angles=(0.0, np.pi / 4, np.pi / 2)):
    r"""Plane strain components from three strain-gauge readings.

    A gauge at angle :math:`\theta` from the :math:`x` axis reads

    .. math::

        \varepsilon_\theta = \varepsilon_x\cos^2\theta + \varepsilon_y\sin^2\theta
                              + \gamma_{xy}\sin\theta\cos\theta,

    and three such equations are solved for :math:`(\varepsilon_x,
    \varepsilon_y, \gamma_{xy})`.

    Parameters
    ----------
    strains : array_like, shape (3,) or (3, k)
        Readings of the three gauges [-], optionally ``k`` readings each.
    angles : array_like of 3 float, optional
        Gauge directions [rad], counterclockwise from :math:`x`, distinct
        modulo :math:`\pi`. The default is the rectangular (0°, 45°, 90°)
        rosette; ``(0, pi/3, 2*pi/3)`` is the delta rosette.

    Returns
    -------
    tuple of float or ndarray
        ``(eps_x, eps_y, gamma_xy)`` [-], with :math:`\gamma_{xy}` the
        engineering shear strain.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. The principal strains are
    ``mohr_circle(eps_x, eps_y, gamma_xy / 2)``.

    References
    ----------
    E. E. Simmons, US Patent 2,292,549 (11 August 1942), the bonded
    resistance strain gauge, found independently by A. C. Ruge at MIT in 1938.
    """
    theta = _vector("angles", angles)
    if theta.size != 3:
        raise ValueError("three gauge angles are required")
    readings = np.asarray(strains, float)
    if readings.ndim not in (1, 2) or readings.shape[0] != 3:
        raise ValueError("strains must have shape (3,) or (3, k)")
    _finite(strains=readings)
    c, s = np.cos(theta), np.sin(theta)
    A = np.column_stack([c**2, s**2, s * c])
    if abs(np.linalg.det(A)) < 1e-12:
        raise ValueError("gauge angles must be distinct modulo pi")
    ex, ey, gxy = np.linalg.solve(A, readings)
    if readings.ndim == 1:
        return float(ex), float(ey), float(gxy)
    return ex, ey, gxy


def plane_stress_from_strains(ex, ey, gxy, E, nu):
    r"""Plane stress from measured strains by Hooke's law for an isotropic material.

    .. math::

        \sigma_x = \frac{E}{1 - \nu^2}(\varepsilon_x + \nu\varepsilon_y),\quad
        \sigma_y = \frac{E}{1 - \nu^2}(\varepsilon_y + \nu\varepsilon_x),\quad
        \tau_{xy} = \frac{E}{2(1 + \nu)}\gamma_{xy}.

    Parameters
    ----------
    ex, ey : float or array_like
        Normal strains [-].
    gxy : float or array_like
        Engineering shear strain [-].
    E : float or array_like
        Young's modulus [Pa], positive.
    nu : float or array_like
        Poisson's ratio [-], in (-1, 0.5).

    Returns
    -------
    tuple of float or ndarray
        ``(sigma_x, sigma_y, tau_xy)`` [Pa].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. Valid at a free surface, where the gauges are.
    """
    _finite(ex=ex, ey=ey, gxy=gxy)
    _positive(E=E)
    _interval("nu", nu, -1, 0.5, closed=False)
    ex, ey, nu = np.asarray(ex, float), np.asarray(ey, float), np.asarray(nu, float)
    factor = np.asarray(E) / (1 - nu**2)
    return factor * (ex + nu * ey), factor * (ey + nu * ex), np.asarray(E) / (2 * (1 + nu)) * gxy
