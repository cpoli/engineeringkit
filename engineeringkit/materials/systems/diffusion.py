r"""Diffusion: Fick's laws, the Arrhenius diffusivity and carburizing.

Fick (1855) borrowed Fourier's law of heat conduction for dissolved salt:
the flux runs down the concentration gradient, :math:`J = -D\,\partial C/\partial x`,
and conservation turns it into :math:`\partial C/\partial t = D\,\partial^2 C/\partial x^2`.
Into a semi-infinite solid whose surface is held at :math:`C_s` the second law
gives the error-function profile used to time carburizing and doping.
"""

from __future__ import annotations

import numpy as np
from scipy.special import erf, erfinv

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.constants import UNIVERSAL_GAS_CONSTANT

__all__ = ["diffusion_coefficient", "fick_flux", "semi_infinite_concentration", "diffusion_time"]


def diffusion_coefficient(D0, Qd, T):
    r"""Arrhenius diffusion coefficient :math:`D = D_0\,e^{-Q_d/RT}`.

    Parameters
    ----------
    D0 : float or array_like
        Pre-exponential factor [m^2/s], positive.
    Qd : float or array_like
        Activation energy for diffusion [J/mol], nonnegative.
    T : float or array_like
        Absolute temperature [K], positive.

    Returns
    -------
    float or ndarray
        Diffusion coefficient [m^2/s].

    Notes
    -----
    Uses the molar gas constant 8.314 J/(mol K). Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/materials`.

    Examples
    --------
    >>> f"{float(diffusion_coefficient(1.2e-4, 130e3, 823.0)):.2e}"
    '6.74e-13'
    """
    _positive(D0=D0, T=T)
    _nonnegative(Qd=Qd)
    return D0 * np.exp(-np.asarray(Qd) / (UNIVERSAL_GAS_CONSTANT * np.asarray(T)))


def fick_flux(D, concentration_gradient):
    r"""Fick's first law, the steady diffusion flux :math:`J = -D\,dC/dx`.

    Parameters
    ----------
    D : float or array_like
        Diffusion coefficient [m^2/s], positive.
    concentration_gradient : float or array_like
        :math:`dC/dx` [kg/m^4] (or [mol/m^4]).

    Returns
    -------
    float or ndarray
        Flux [kg/(m^2 s)] (or [mol/(m^2 s)]), positive down the gradient.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    A. Fick, "Ueber Diffusion," *Ann. Phys.* 170 (1855) 59-86.

    Examples
    --------
    >>> f"{float(fick_flux(3e-11, (0.8 - 1.2) / 5e-3)):.1e}"
    '2.4e-09'
    """
    _positive(D=D)
    _finite(concentration_gradient=concentration_gradient)
    return -np.asarray(D) * concentration_gradient


def semi_infinite_concentration(x, t, D, C0, Cs):
    r"""Fick's second law in a semi-infinite solid with a fixed surface concentration.

    .. math:: \frac{C(x, t) - C_0}{C_s - C_0} = 1 - \operatorname{erf}\frac{x}{2\sqrt{Dt}}

    Parameters
    ----------
    x : float or array_like
        Depth below the surface [m], nonnegative.
    t : float or array_like
        Time [s], positive.
    D : float or array_like
        Diffusion coefficient [m^2/s], positive, independent of concentration.
    C0 : float or array_like
        Initial uniform concentration [any unit].
    Cs : float or array_like
        Surface concentration held from ``t = 0`` [same unit].

    Returns
    -------
    float or ndarray
        Concentration at depth ``x`` [same unit as ``C0``].

    Notes
    -----
    The solid must be thicker than a few :math:`\sqrt{Dt}`. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/materials`.

    Examples
    --------
    >>> round(float(semi_infinite_concentration(0.0, 1.0, 1e-11, 0.2, 1.0)), 6)
    1.0
    """
    _nonnegative(x=x)
    _positive(t=t, D=D)
    _finite(C0=C0, Cs=Cs)
    z = np.asarray(x) / (2 * np.sqrt(np.asarray(D) * t))
    return Cs - (np.asarray(Cs) - C0) * erf(z)


def diffusion_time(x, concentration, D, C0, Cs):
    r"""Time for a semi-infinite solid to reach ``concentration`` at depth ``x``.

    Inverts :func:`semi_infinite_concentration`:
    :math:`t = x^2 / (4Dz^2)` with :math:`\operatorname{erf} z = (C_s - C)/(C_s - C_0)`.

    Parameters
    ----------
    x : float or array_like
        Depth [m], positive.
    concentration : float or array_like
        Target concentration [any unit], strictly between ``C0`` and ``Cs``.
    D : float or array_like
        Diffusion coefficient [m^2/s], positive.
    C0, Cs : float or array_like
        Initial and surface concentrations [same unit], different.

    Returns
    -------
    float or ndarray
        Time [s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    Examples
    --------
    >>> round(float(diffusion_time(5e-4, 0.80, 1.6e-11, 0.25, 1.20)) / 3600, 2)
    7.05
    """
    _positive(x=x, D=D)
    _finite(concentration=concentration, C0=C0, Cs=Cs)
    if np.any(np.asarray(Cs) == C0):
        raise ValueError("C0 and Cs must differ")
    ratio = (np.asarray(Cs) - concentration) / (np.asarray(Cs) - C0)
    if np.any((ratio <= 0) | (ratio >= 1)):
        raise ValueError("concentration must lie strictly between C0 and Cs")
    z = erfinv(ratio)
    return np.asarray(x) ** 2 / (4 * np.asarray(D) * z * z)
