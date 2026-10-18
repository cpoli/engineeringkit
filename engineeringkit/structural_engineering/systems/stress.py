r"""Plane stress transformation and Mohr's circle.

Rotating a plane-stress element by :math:`\theta` moves its stress point
around a circle of centre :math:`(\sigma_x + \sigma_y)/2` and radius
:math:`\sqrt{((\sigma_x - \sigma_y)/2)^2 + \tau_{xy}^2}` through :math:`2\theta`.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite

__all__ = ["mohr_circle", "stress_transformation"]


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
