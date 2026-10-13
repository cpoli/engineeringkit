r"""Reynolds-number similarity and flat-plate boundary layers.

Reynolds (1883) showed that the ratio of inertial to viscous forces governs
the transition from laminar to turbulent flow; Prandtl (1904) confined the
viscous effects to a thin boundary layer, which Blasius (1908) solved exactly
for a flat plate.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["reynolds_number", "skin_friction_laminar", "skin_friction_turbulent", "boundary_layer_laminar"]


def reynolds_number(rho, V, L, mu):
    r"""Reynolds number :math:`Re = \rho V L / \mu`.

    Parameters
    ----------
    rho : float or array_like
        Density [kg/m^3], positive.
    V : float or array_like
        Speed [m/s], nonnegative.
    L : float or array_like
        Reference length (chord, plate distance) [m], positive.
    mu : float or array_like
        Dynamic viscosity [Pa s], positive.

    Returns
    -------
    float or ndarray
        Reynolds number [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(reynolds_number(1.0, 10.0, 1.0, 1e-5)))
    1000000
    """
    _positive(rho=rho, L=L, mu=mu)
    _nonnegative(V=V)
    return rho * np.asarray(V) * L / mu


def skin_friction_laminar(Re):
    r"""Blasius average skin-friction coefficient of a laminar flat plate, :math:`C_f = 1.328/\sqrt{Re_L}`.

    Parameters
    ----------
    Re : float or array_like
        Plate-length Reynolds number [-], positive; laminar flow
        (typically :math:`Re_L < 5\times10^5`) is assumed, not checked.

    Returns
    -------
    float or ndarray
        Average skin-friction coefficient (both sides not included) [-].

    Notes
    -----
    Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    References
    ----------
    H. Blasius, "Grenzschichten in Flüssigkeiten mit kleiner Reibung,"
    *Z. Math. Phys.* 56 (1908) 1-37.

    Examples
    --------
    >>> round(float(skin_friction_laminar(1e6)), 6)
    0.001328
    """
    _positive(Re=Re)
    return 1.328 / np.sqrt(Re)


def skin_friction_turbulent(Re):
    r"""Prandtl's one-fifth-power turbulent flat-plate skin friction, :math:`C_f = 0.074/Re_L^{1/5}`.

    Parameters
    ----------
    Re : float or array_like
        Plate-length Reynolds number [-], positive; fully turbulent flow
        from the leading edge, roughly :math:`5\times10^5 < Re_L < 10^7`.

    Returns
    -------
    float or ndarray
        Average skin-friction coefficient [-].

    Notes
    -----
    Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    References
    ----------
    H. Schlichting and K. Gersten, *Boundary-Layer Theory*, 9th ed.
    (Springer, 2017), ch. 18.

    Examples
    --------
    >>> round(float(skin_friction_turbulent(1e5)), 6)
    0.0074
    """
    _positive(Re=Re)
    return 0.074 / np.asarray(Re) ** 0.2


def boundary_layer_laminar(x, Re_x):
    r"""Laminar boundary-layer thickness :math:`\delta \approx 5x/\sqrt{Re_x}`.

    Parameters
    ----------
    x : float or array_like
        Distance from the leading edge [m], nonnegative.
    Re_x : float or array_like
        Local Reynolds number :math:`\rho V x/\mu` [-], positive.

    Returns
    -------
    float or ndarray
        99%-velocity thickness [m] (Blasius gives 4.91; 5 is the textbook rounding).

    Notes
    -----
    Inputs are raw SI numbers; see :doc:`/guides/aerodynamics`.

    Examples
    --------
    >>> round(float(boundary_layer_laminar(1.0, 1e6)), 6)
    0.005
    """
    _nonnegative(x=x)
    _positive(Re_x=Re_x)
    return 5 * np.asarray(x) / np.sqrt(Re_x)
