r"""Hydraulic transients: Joukowsky's water hammer.

Closing a valve stops the water next to it, and a pressure wave carries the
news upstream at the speed of sound in the pipe. Joukowsky (1898) measured
the jump :math:`\Delta p = \rho a \Delta v` in the Moscow water mains and
used Korteweg's formula (1878) for the wave speed in an elastic pipe. A
closure is "fast" when it takes less than the round trip :math:`2L/a`.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["water_hammer_pressure", "pipe_wave_speed", "critical_closure_time"]


def water_hammer_pressure(rho, wave_speed, delta_v):
    r"""Joukowsky pressure rise for a sudden velocity change, :math:`\Delta p = \rho a\,\Delta v`.

    Parameters
    ----------
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    wave_speed : float or array_like
        Pressure-wave speed in the pipe [m/s], positive (see :func:`pipe_wave_speed`).
    delta_v : float or array_like
        Velocity reduction [m/s]; negative for a sudden opening.

    Returns
    -------
    float or ndarray
        Pressure rise [Pa].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Holds for closures shorter than
    :math:`2L/a` and neglects friction.

    References
    ----------
    N. Joukowsky, "Über den hydraulischen Stoss in Wasserleitungsröhren,"
    *Mém. Acad. Imp. Sci. St.-Pétersbourg*, 8th ser., 9(5) (1900) 1-71
    (read 1898); English translation by O. Simin, *Proc. AWWA* 24 (1904)
    341-424.

    Examples
    --------
    >>> water_hammer_pressure(1000.0, 1200.0, 1.0)
    1200000.0
    """
    _finite(rho=rho, wave_speed=wave_speed, delta_v=delta_v)
    _positive(rho=rho, wave_speed=wave_speed)
    return rho * wave_speed * delta_v


def pipe_wave_speed(bulk_modulus, rho, D, wall_thickness, E, c1=1.0):
    r"""Korteweg's pressure-wave speed in a liquid-filled elastic pipe.

    .. math::

        a = \sqrt{\frac{K/\rho}{1 + c_1\,K D/(E e)}}

    Parameters
    ----------
    bulk_modulus : float or array_like
        Bulk modulus of the liquid :math:`K` [Pa], positive (2.2 GPa for water).
    rho : float or array_like
        Liquid density [kg/m^3], positive.
    D : float or array_like
        Pipe inner diameter [m], positive.
    wall_thickness : float or array_like
        Wall thickness :math:`e` [m], positive.
    E : float or array_like
        Young's modulus of the wall [Pa], positive.
    c1 : float or array_like, optional
        Pipe-support factor [-], nonnegative: 1 for expansion joints
        throughout, :math:`1 - \nu^2` for a pipe anchored against axial
        movement, 0 for a rigid pipe.

    Returns
    -------
    float or ndarray
        Wave speed [m/s]; at most the speed of sound :math:`\sqrt{K/\rho}` in the liquid.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Thin-walled pipe, :math:`D/e \gtrsim 25`.

    References
    ----------
    D. J. Korteweg, "Ueber die Fortpflanzungsgeschwindigkeit des Schalles in
    elastischen Röhren," *Ann. Phys.* 241 (1878) 525-542.

    Examples
    --------
    >>> round(float(pipe_wave_speed(2.2e9, 1000.0, 0.5, 0.01, 200e9)))
    1191
    """
    _finite(bulk_modulus=bulk_modulus, rho=rho, D=D, wall_thickness=wall_thickness, E=E, c1=c1)
    _positive(bulk_modulus=bulk_modulus, rho=rho, D=D, wall_thickness=wall_thickness, E=E)
    _nonnegative(c1=c1)
    stiffness = np.asarray(E, dtype=float) * np.asarray(wall_thickness, dtype=float)
    return np.sqrt((bulk_modulus / rho) / (1 + c1 * bulk_modulus * np.asarray(D, dtype=float) / stiffness))


def critical_closure_time(length, wave_speed):
    r"""Round-trip time :math:`2L/a` of a pressure wave along a pipe.

    Parameters
    ----------
    length : float or array_like
        Pipe length from valve to reservoir [m], positive.
    wave_speed : float or array_like
        Pressure-wave speed [m/s], positive.

    Returns
    -------
    float or ndarray
        Time [s]; a valve closed faster than this sees the full Joukowsky rise.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> critical_closure_time(600.0, 1200.0)
    1.0
    """
    _finite(length=length, wave_speed=wave_speed)
    _positive(length=length, wave_speed=wave_speed)
    return 2 * length / wave_speed
