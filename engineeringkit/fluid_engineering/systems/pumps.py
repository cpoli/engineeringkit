r"""Pumps and cavitation: hydraulic power, specific speed and Thoma's parameter.

A pump adds head :math:`H` to a flow :math:`Q`, so the power given to the
fluid is :math:`\rho g Q H`. Where the local pressure falls to the vapour
pressure the liquid boils and the collapsing bubbles pit the blades. Thoma
(1924) found that turbines of one design all start to cavitate at the same
ratio :math:`\sigma = \mathrm{NPSH}/H` of the suction head margin to the
machine head.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.fluid_engineering.utils.constants import STANDARD_GRAVITY as G

__all__ = [
    "pump_hydraulic_power",
    "pump_input_power",
    "pump_specific_speed",
    "cavitation_number",
    "npsh_available",
    "thoma_sigma",
]


def pump_hydraulic_power(rho, Q, head, g=G):
    r"""Power delivered to the fluid, :math:`P = \rho g Q H`.

    Parameters
    ----------
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    Q : float or array_like
        Volume flow [m^3/s], nonnegative.
    head : float or array_like
        Head added by the pump [m], nonnegative.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Hydraulic power [W].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> pump_hydraulic_power(1000.0, 0.01, 10.0, g=10.0)
    1000.0
    """
    _finite(rho=rho, Q=Q, head=head, g=g)
    _positive(rho=rho, g=g)
    _nonnegative(Q=Q, head=head)
    return rho * g * Q * head


def pump_input_power(rho, Q, head, efficiency, g=G):
    r"""Shaft power a pump needs, :math:`P = \rho g Q H/\eta`.

    Parameters
    ----------
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    Q : float or array_like
        Volume flow [m^3/s], nonnegative.
    head : float or array_like
        Head added by the pump [m], nonnegative.
    efficiency : float or array_like
        Pump efficiency [-], in (0, 1].
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Shaft power [W].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> pump_input_power(1000.0, 0.01, 10.0, 0.5, g=10.0)
    2000.0
    """
    _finite(rho=rho, Q=Q, head=head, efficiency=efficiency, g=g)
    _interval("efficiency", efficiency, 0, 1)
    if np.any(np.asarray(efficiency) == 0):
        raise ValueError("efficiency must be positive")
    return pump_hydraulic_power(rho, Q, head, g) / efficiency


def pump_specific_speed(omega, Q, head, g=G):
    r"""Dimensionless specific speed :math:`\Omega_s = \omega\sqrt Q/(gH)^{3/4}`.

    Parameters
    ----------
    omega : float or array_like
        Shaft speed [rad/s], nonnegative.
    Q : float or array_like
        Volume flow at best efficiency [m^3/s], nonnegative.
    head : float or array_like
        Head at best efficiency [m], positive.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Specific speed [-]: below about 1 for radial pumps, above about 4
        for axial ones.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Industry values in rpm, gpm and feet
    differ by a constant factor.

    Examples
    --------
    >>> float(pump_specific_speed(10.0, 4.0, 1.0, g=1.0))
    20.0
    """
    _finite(omega=omega, Q=Q, head=head, g=g)
    _positive(head=head, g=g)
    _nonnegative(Q=Q, omega=omega)
    return omega * np.sqrt(Q) / (g * head) ** 0.75


def cavitation_number(p, p_vapor, rho, v):
    r"""Cavitation number :math:`\mathrm{Ca} = (p - p_v)/(\tfrac12\rho v^2)`.

    Parameters
    ----------
    p : float or array_like
        Reference absolute pressure [Pa].
    p_vapor : float or array_like
        Vapour pressure of the liquid [Pa].
    rho : float or array_like
        Liquid density [kg/m^3], positive.
    v : float or array_like
        Reference speed [m/s], nonzero.

    Returns
    -------
    float or ndarray
        Cavitation number [-]; negative when the pressure is below the
        vapour pressure.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> cavitation_number(3000.0, 1000.0, 1000.0, 2.0)
    1.0
    """
    _finite(p=p, p_vapor=p_vapor, rho=rho, v=v)
    _positive(rho=rho)
    if np.any(np.asarray(v) == 0):
        raise ValueError("velocity must be nonzero")
    return (p - p_vapor) / (0.5 * rho * v * v)


def npsh_available(p_inlet, v_inlet, p_vapor, rho, g=G):
    r"""Net positive suction head available at a pump inlet.

    .. math::

        \mathrm{NPSH}_A = \frac{p - p_v}{\rho g} + \frac{v^2}{2g}

    Parameters
    ----------
    p_inlet : float or array_like
        Absolute static pressure at the pump inlet [Pa], positive.
    v_inlet : float or array_like
        Mean velocity at the inlet [m/s].
    p_vapor : float or array_like
        Vapour pressure of the liquid [Pa], nonnegative.
    rho : float or array_like
        Liquid density [kg/m^3], positive.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        NPSH available [m]; the pump cavitates when it falls below the
        NPSH it requires.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. For a pump drawing from an open sump,
    :math:`p/(\rho g) + v^2/(2g) = p_{atm}/(\rho g) - z_s - h_{L}`, the
    atmosphere less the suction lift and the suction-line losses.

    Examples
    --------
    >>> npsh_available(103000.0, 0.0, 3000.0, 1000.0, g=10.0)
    10.0
    """
    _finite(p_inlet=p_inlet, v_inlet=v_inlet, p_vapor=p_vapor, rho=rho, g=g)
    _positive(p_inlet=p_inlet, rho=rho, g=g)
    _nonnegative(p_vapor=p_vapor)
    return (p_inlet - p_vapor) / (rho * g) + v_inlet * v_inlet / (2 * g)


def thoma_sigma(npsh, head):
    r"""Thoma's cavitation parameter :math:`\sigma = \mathrm{NPSH}/H`.

    Parameters
    ----------
    npsh : float or array_like
        Net positive suction head [m].
    head : float or array_like
        Head of the machine [m], positive.

    Returns
    -------
    float or ndarray
        Thoma's :math:`\sigma` [-]. A machine runs free of cavitation while
        its plant :math:`\sigma` stays above the critical value found by
        test for its design.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    References
    ----------
    D. Thoma, "Die Kavitation bei Wasserturbinen," *Wasserkraft-Jahrbuch*
    (Munich, 1924).

    Examples
    --------
    >>> thoma_sigma(5.0, 50.0)
    0.1
    """
    _finite(npsh=npsh, head=head)
    _positive(head=head)
    return npsh / head
