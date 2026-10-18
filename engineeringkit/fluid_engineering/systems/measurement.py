r"""Flow measurement: Pitot's tube and Venturi's contraction.

Both instruments turn a speed into a pressure difference through Bernoulli's
equation. Pitot's bent tube (1732) faces the flow and reads the stagnation
rise :math:`\tfrac12\rho v^2`. Venturi's contraction (1797), made into a
water meter by Herschel (1887), speeds the flow up in a throat, and the
pressure drop there, with continuity, gives the volume flow.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["pitot_velocity", "venturi_flow"]


def pitot_velocity(delta_p, rho):
    r"""Flow speed from a Pitot-static reading, :math:`v = \sqrt{2\Delta p/\rho}`.

    Parameters
    ----------
    delta_p : float or array_like
        Stagnation minus static pressure [Pa], nonnegative.
    rho : float or array_like
        Fluid density [kg/m^3], positive.

    Returns
    -------
    float or ndarray
        Flow speed [m/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Incompressible flow; for gases the
    error is below 1% up to about Mach 0.2.

    References
    ----------
    H. Pitot, "Description d'une machine pour mesurer la vitesse des eaux
    courantes et le sillage des vaisseaux," *Hist. Acad. R. Sci.* (1732)
    363-376.

    Examples
    --------
    >>> float(pitot_velocity(2000.0, 1000.0))
    2.0
    """
    _finite(delta_p=delta_p, rho=rho)
    _nonnegative(delta_p=delta_p)
    _positive(rho=rho)
    return np.sqrt(2 * delta_p / rho)


def venturi_flow(delta_p, rho, A1, A2, discharge_coefficient=1.0):
    r"""Volume flow through a Venturi meter.

    .. math::

        Q = C_d\,A_2\sqrt{\frac{2\Delta p}{\rho\,(1 - (A_2/A_1)^2)}}

    Parameters
    ----------
    delta_p : float or array_like
        Inlet minus throat pressure [Pa], nonnegative.
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    A1 : float or array_like
        Inlet area [m^2], positive.
    A2 : float or array_like
        Throat area [m^2], positive and smaller than ``A1``.
    discharge_coefficient : float or array_like, optional
        :math:`C_d` in (0, 1]; the default 1 is the ideal meter, and 0.98 is
        typical of a machined Venturi tube.

    Returns
    -------
    float or ndarray
        Volume flow [m^3/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Horizontal meter and incompressible
    flow.

    References
    ----------
    G. B. Venturi, *Recherches expérimentales sur le principe de la
    communication latérale du mouvement dans les fluides* (Paris, 1797);
    C. Herschel, "The Venturi water meter," *Trans. ASCE* 17 (1887) 228-258.

    Examples
    --------
    >>> round(float(venturi_flow(15000.0, 1000.0, 0.02, 0.01)), 4)
    0.0632
    """
    _finite(delta_p=delta_p, rho=rho, A1=A1, A2=A2)
    _nonnegative(delta_p=delta_p)
    _positive(rho=rho, A1=A1, A2=A2)
    _interval("discharge_coefficient", discharge_coefficient, 0, 1, closed="right")
    if np.any(np.asarray(A2) >= np.asarray(A1)):
        raise ValueError("throat A2 must be smaller than A1")
    return discharge_coefficient * A1 * A2 * np.sqrt(2 * delta_p / (rho * (A1 * A1 - A2 * A2)))
