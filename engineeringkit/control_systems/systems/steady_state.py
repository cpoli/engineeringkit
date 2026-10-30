r"""System type, static error constants and steady-state tracking errors.

Write a unity-feedback loop as :math:`L(s) = K_0 N_0(s)/(s^N D_0(s))` with
:math:`N_0(0) = D_0(0) = 1`. The number :math:`N` of integrators is the
*system type*. The final-value theorem gives the error of a stable loop
to the test inputs :math:`1/s`, :math:`1/s^2` and :math:`1/s^3` through
the position, velocity and acceleration constants
:math:`K_p = \lim L`, :math:`K_v = \lim sL` and :math:`K_a = \lim s^2L`
as :math:`s \to 0`.
"""

from __future__ import annotations

import numpy as np

from engineeringkit.control_systems.systems.transfer import feedback, poles
from engineeringkit.control_systems.utils.polynomials import num_den

__all__ = ["system_type", "error_constants"]


def _trailing_zeros(p):
    return len(p) - len(np.trim_zeros(p, "b"))


def system_type(loop):
    """Number of free integrators (poles at :math:`s = 0`) in an open loop.

    Parameters
    ----------
    loop : scipy.signal.lti
        SISO open-loop model :math:`L(s)`.

    Returns
    -------
    int
        Poles at the origin minus zeros at the origin; negative when the loop
        differentiates.

    Notes
    -----
    Only exact zero coefficients count, so a pole at :math:`10^{-9}` is not
    an integrator. See :doc:`/guides/control_systems`.

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> system_type(transfer_function([1, 2], [1, 3, 0, 0]))
    2
    """
    num, den = num_den(loop)
    return _trailing_zeros(den) - _trailing_zeros(num)


def error_constants(loop):
    r"""Static error constants and steady-state errors of a stable unity-feedback loop.

    Parameters
    ----------
    loop : scipy.signal.lti
        SISO open-loop model :math:`L(s)`; the closed loop :math:`L/(1+L)`
        must be stable.

    Returns
    -------
    dict
        ``type`` (int); ``Kp`` [-], ``Kv`` [1/s] and ``Ka`` [1/s²]; and the
        steady-state errors ``step_error`` :math:`1/(1+K_p)`,
        ``ramp_error`` :math:`1/K_v` [s] and ``parabola_error``
        :math:`1/K_a` [s²], per unit step, unit-slope ramp and
        :math:`t^2/2` input. Constants are ``inf`` above the system type and
        ``0`` below it; errors are then ``0`` and ``inf``.

    Raises
    ------
    ValueError
        If the closed loop is not asymptotically stable, so that no steady
        state exists.

    Notes
    -----
    A type-:math:`N` loop follows a polynomial input of degree below
    :math:`N` with zero error, one of degree :math:`N` with a finite error,
    and loses higher ones. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/control_systems`.

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> e = error_constants(transfer_function([10], [1, 5, 0]))
    >>> e["type"], e["Kv"], e["ramp_error"]
    (1, 2.0, 0.5)
    """
    if np.any(poles(feedback(loop)).real >= 0):
        raise ValueError("the closed loop is not stable, so it has no steady-state error")
    num, den = num_den(loop)
    n = system_type(loop)
    k0 = float(np.trim_zeros(num, "b")[-1] / np.trim_zeros(den, "b")[-1])
    constants = [k0 if n == k else (np.inf if n > k else 0.0) for k in range(3)]
    kp, kv, ka = constants
    return {
        "type": n,
        "Kp": kp,
        "Kv": kv,
        "Ka": ka,
        "step_error": 1 / (1 + kp),
        "ramp_error": np.inf if kv == 0 else 1 / kv,
        "parabola_error": np.inf if ka == 0 else 1 / ka,
    }
