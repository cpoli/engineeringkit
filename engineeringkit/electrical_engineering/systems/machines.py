r"""Electric machines: the DC machine and the induction motor equivalent circuits.

The Hopkinsons reduced the dynamo to a magnetic circuit and an armature
whose EMF is proportional to flux and speed, :math:`E = k\Phi\omega`.
Steinmetz treated the induction motor as a transformer whose secondary
resistance :math:`R_2/s` depends on the slip, which gives its torque-speed
curve and its breakdown torque.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["dc_motor_torque", "dc_motor_speed", "induction_motor_torque", "induction_motor_breakdown"]


def dc_motor_torque(v_terminal, armature_resistance, k_phi, omega):
    r"""Torque :math:`T = k\Phi(V - k\Phi\omega)/R_a` of a separately excited DC machine.

    The armature is its back EMF :math:`E = k\Phi\omega` behind the
    armature resistance, so the current is :math:`(V - E)/R_a` and the
    torque :math:`k\Phi I_a`. The torque-speed curve is a straight line from
    the stall torque :math:`k\Phi V/R_a` to the no-load speed
    :math:`V/k\Phi`.

    Parameters
    ----------
    v_terminal : float or array_like
        Armature terminal voltage [V].
    armature_resistance : float or array_like
        Armature resistance :math:`R_a` [ohm], positive.
    k_phi : float or array_like
        Machine constant :math:`k\Phi` at the set field current [V s/rad]
        (equally N m/A), positive.
    omega : float or array_like
        Shaft speed [rad/s].

    Returns
    -------
    float or ndarray
        Electromagnetic torque [N m]; negative when the machine runs above
        its no-load speed and generates.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`. Armature reaction, brush drop and
    rotational losses are neglected.

    References
    ----------
    J. Hopkinson and E. Hopkinson, "Dynamo-electric machinery," *Phil.
    Trans. R. Soc.* 177 (1886) 331.

    Examples
    --------
    >>> float(dc_motor_torque(240.0, 0.5, 2.0, 100.0))
    160.0
    """
    _finite(v_terminal=v_terminal, omega=omega)
    _positive(armature_resistance=armature_resistance, k_phi=k_phi)
    return k_phi * (v_terminal - k_phi * np.asarray(omega)) / armature_resistance


def dc_motor_speed(v_terminal, armature_resistance, k_phi, torque):
    r"""Speed :math:`\omega = V/k\Phi - R_a T/(k\Phi)^2` of a separately excited DC motor; the inverse of :func:`dc_motor_torque`.

    Parameters
    ----------
    v_terminal : float or array_like
        Armature terminal voltage [V].
    armature_resistance : float or array_like
        Armature resistance [ohm], positive.
    k_phi : float or array_like
        Machine constant [V s/rad], positive.
    torque : float or array_like
        Load torque [N m].

    Returns
    -------
    float or ndarray
        Shaft speed [rad/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> float(dc_motor_speed(240.0, 0.5, 2.0, 160.0))
    100.0
    """
    _finite(v_terminal=v_terminal, torque=torque)
    _positive(armature_resistance=armature_resistance, k_phi=k_phi)
    return v_terminal / k_phi - armature_resistance * np.asarray(torque) / k_phi**2


def _thevenin(v_phase, r1, x1, xm):
    """Stator and magnetizing branch seen from the rotor: (V_th, R_th, X_th)."""
    v_phase, r1, x1, xm = (np.asarray(x, dtype=float) for x in (v_phase, r1, x1, xm))
    z = 1j * xm * (r1 + 1j * x1) / (r1 + 1j * (x1 + xm))
    return np.abs(v_phase * 1j * xm / (r1 + 1j * (x1 + xm))), z.real, z.imag


def _check_motor(v_phase, synchronous_speed, r1, x1, r2, x2, xm, phases):
    _finite(v_phase=v_phase)
    _nonnegative(v_phase=v_phase, r1=r1)
    _positive(synchronous_speed=synchronous_speed, x1=x1, r2=r2, x2=x2, xm=xm)
    _integer("phases", phases, minimum=1)


def induction_motor_torque(v_phase, slip, synchronous_speed, r1, x1, r2, x2, xm, phases=3):
    r"""Torque of an induction motor from Steinmetz's per-phase equivalent circuit.

    Replacing the stator impedance and magnetizing reactance by their
    Thévenin equivalent :math:`V_{th}, R_{th} + jX_{th}` seen from the rotor,

    .. math::

        T = \frac{m V_{th}^2\, R_2/s}{\omega_s\left[(R_{th} + R_2/s)^2 + (X_{th} + X_2)^2\right]},

    computed in the equivalent form multiplied through by :math:`s^2` so
    that :math:`s = 0` gives zero torque.

    Parameters
    ----------
    v_phase : float or array_like
        RMS phase voltage at the stator [V], nonnegative.
    slip : float or array_like
        Slip :math:`s = (\omega_s - \omega)/\omega_s` [-]: 0 at synchronous
        speed, 1 at standstill, negative when generating and above 1 when
        plugging.
    synchronous_speed : float or array_like
        Synchronous shaft speed :math:`\omega_s` [rad/s], positive.
    r1, x1 : float or array_like
        Stator resistance [ohm] (nonnegative) and leakage reactance [ohm]
        (positive).
    r2, x2 : float or array_like
        Rotor resistance and leakage reactance referred to the stator
        [ohm], positive.
    xm : float or array_like
        Magnetizing reactance [ohm], positive.
    phases : int, optional
        Number of phases :math:`m`, at least one.

    Returns
    -------
    float or ndarray
        Electromagnetic torque [N m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`. Core loss is neglected (no
    shunt resistance) and the parameters are constant: deep-bar and
    double-cage rotors, whose :math:`R_2` rises with slip, are not modelled.

    References
    ----------
    C. P. Steinmetz, "The alternating current induction motor," *Trans.
    AIEE* 14 (1897) 183-217.

    Examples
    --------
    >>> round(float(induction_motor_torque(230.0, 0.03, 2 * np.pi * 25, 0.3, 0.6, 0.25, 0.6, 20.0)), 1)
    104.9
    """
    _check_motor(v_phase, synchronous_speed, r1, x1, r2, x2, xm, phases)
    _finite(slip=slip)
    s = np.asarray(slip, dtype=float)
    v_th, r_th, x_th = _thevenin(v_phase, r1, x1, xm)
    return phases * v_th**2 * r2 * s / (synchronous_speed * ((r_th * s + r2) ** 2 + (s * (x_th + x2)) ** 2))


def induction_motor_breakdown(v_phase, synchronous_speed, r1, x1, r2, x2, xm, phases=3):
    r"""Slip and torque at the breakdown (pull-out) point of an induction motor.

    Torque is greatest when the rotor resistance :math:`R_2/s` matches the
    impedance in front of it, at

    .. math::

        s_{max} = \frac{R_2}{\sqrt{R_{th}^2 + (X_{th} + X_2)^2}},\qquad
        T_{max} = \frac{m V_{th}^2}{2\omega_s\left[R_{th} + \sqrt{R_{th}^2 + (X_{th} + X_2)^2}\right]},

    so the breakdown torque does not depend on :math:`R_2`, only where it
    occurs.

    Parameters
    ----------
    v_phase : float or array_like
        RMS phase voltage [V], nonnegative.
    synchronous_speed : float or array_like
        Synchronous shaft speed [rad/s], positive.
    r1, x1, r2, x2, xm : float or array_like
        Equivalent-circuit parameters [ohm]; see :func:`induction_motor_torque`.
    phases : int, optional
        Number of phases, at least one.

    Returns
    -------
    tuple of float or ndarray
        ``(slip_max, torque_max)``: slip [-] and torque [N m] at breakdown.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> s, T = induction_motor_breakdown(230.0, 2 * np.pi * 25, 0.3, 0.6, 0.25, 0.6, 20.0)
    >>> round(float(s), 4), round(float(T), 1)
    (0.2049, 316.8)
    """
    _check_motor(v_phase, synchronous_speed, r1, x1, r2, x2, xm, phases)
    v_th, r_th, x_th = _thevenin(v_phase, r1, x1, xm)
    z = np.hypot(r_th, x_th + x2)
    return r2 / z, phases * v_th**2 / (2 * synchronous_speed * (r_th + z))
