r"""Proportional-integral-derivative control and Ziegler-Nichols tuning.

Minorsky (1922) analysed automatic ship steering and showed that a
controller combining the heading error, its integral and its rate of change,
:math:`u = K_p e + K_i\int e\,dt + K_d\dot e`, both holds course and removes
the steady offset caused by a constant disturbance. Ziegler and Nichols
(1942) gave rules for choosing the gains from a single experiment.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive
from engineeringkit.control_systems.core.base import ResponseResult
from engineeringkit.control_systems.utils.polynomials import num_den

__all__ = ["PID", "ziegler_nichols"]


@dataclass(frozen=True)
class PID:
    r"""Ideal parallel PID controller :math:`C(s) = K_p + K_i/s + K_d s`.

    Derivative action is ideal (no roll-off filter), and there is no actuator
    saturation or anti-windup. With ``ki == 0`` no integrator is included.

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> r = PID(1, 1).simulate(transfer_function([1], [1, 1]), np.linspace(0, 5, 51))
    >>> round(r.summary()["final"], 4)
    0.9933
    """

    kp: float
    """float: Proportional gain [input/error]."""

    ki: float = 0.0
    """float: Integral gain [input/(error s)]."""

    kd: float = 0.0
    """float: Derivative gain [input s/error]."""

    def __post_init__(self):
        _finite(kp=self.kp, ki=self.ki, kd=self.kd)

    def _polynomials(self):
        if self.ki == 0:
            return [self.kd, self.kp], [1]
        return [self.kd, self.kp, self.ki], [1, 0]

    def transfer_function(self):
        """Controller transfer function :math:`C(s)`; all-zero gains have none."""
        cnum, cden = self._polynomials()
        cnum = np.trim_zeros(np.asarray(cnum, float), "f")
        if not cnum.size:
            raise ValueError("a controller with all gains zero has no transfer function")
        return signal.TransferFunction(cnum, cden)

    def closed_loop(self, plant):
        """Unity-feedback closed loop :math:`CG/(1 + CG)` from reference to output."""
        cnum, cden = self._polynomials()
        pnum, pden = num_den(plant)
        ol_num = np.polymul(cnum, pnum)
        ol_den = np.polymul(cden, pden)
        return signal.TransferFunction(ol_num, np.polyadd(ol_den, ol_num))

    def simulate(self, plant, T=None):
        """Unit-step reference response of the closed loop, as a :class:`ResponseResult`."""
        t, y = signal.step(self.closed_loop(plant), T=T)
        return ResponseResult(t, y)


_ZN_RULES = {"P": (0.5, None, None), "PI": (0.45, 1 / 1.2, None), "PID": (0.6, 0.5, 0.125)}


def ziegler_nichols(ku, tu, controller="PID"):
    r"""Ziegler-Nichols ultimate-sensitivity tuning.

    Raise a proportional gain until the loop oscillates steadily: that gain
    is the ultimate gain :math:`K_u` and the oscillation period is
    :math:`T_u`. The classic rules then set

    ======== =============== ================ ================
    type     :math:`K_p`     :math:`T_i`      :math:`T_d`
    ======== =============== ================ ================
    P        :math:`0.5K_u`
    PI       :math:`0.45K_u` :math:`T_u/1.2`
    PID      :math:`0.6K_u`  :math:`T_u/2`    :math:`T_u/8`
    ======== =============== ================ ================

    with :math:`K_i = K_p/T_i` and :math:`K_d = K_p T_d`.

    Parameters
    ----------
    ku : float
        Ultimate proportional gain [input/error], positive.
    tu : float
        Ultimate oscillation period [s], positive.
    controller : {"P", "PI", "PID"}, optional
        Controller structure.

    Returns
    -------
    PID
        Tuned controller.

    Notes
    -----
    The rules aim for roughly quarter-amplitude decay and are only a starting
    point; they give large overshoot on many plants. See
    :doc:`/guides/control_systems`.

    References
    ----------
    J. G. Ziegler and N. B. Nichols, "Optimum settings for automatic
    controllers," *Trans. ASME* 64 (1942) 759-768.

    Examples
    --------
    >>> ziegler_nichols(8.0, 4.0)
    PID(kp=4.8, ki=2.4, kd=2.4)
    """
    _positive(ku=ku, tu=tu)
    if controller not in _ZN_RULES:
        raise ValueError("controller must be 'P', 'PI' or 'PID'")
    kp_ratio, ti_ratio, td_ratio = _ZN_RULES[controller]
    kp = kp_ratio * ku
    ki = 0.0 if ti_ratio is None else kp / (ti_ratio * tu)
    kd = 0.0 if td_ratio is None else kp * td_ratio * tu
    return PID(float(kp), float(ki), float(kd))
