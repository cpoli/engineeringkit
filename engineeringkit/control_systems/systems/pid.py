r"""Proportional-integral-derivative control and Ziegler-Nichols tuning.

Minorsky (1922) analysed automatic ship steering and showed that a
controller combining the heading error, its integral and its rate of change,
:math:`u = K_p e + K_i\int e\,dt + K_d\dot e`, both holds course and removes
the steady offset caused by a constant disturbance. Ziegler and Nichols
(1942) gave rules for choosing the gains from a single experiment. A real
actuator saturates, and while it does the integrator keeps charging on an
error the loop cannot act on; anti-windup stops it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive
from engineeringkit.control_systems.core.base import ResponseResult
from engineeringkit.control_systems.utils.polynomials import num_den
from engineeringkit.integrators import integrate_ode

__all__ = ["PID", "ziegler_nichols"]


@dataclass(frozen=True)
class PID:
    r"""Ideal parallel PID controller :math:`C(s) = K_p + K_i/s + K_d s`.

    Derivative action is ideal (no roll-off filter). The linear methods have
    no actuator limit; :meth:`simulate_saturated` adds one, with optional
    anti-windup. With ``ki == 0`` no integrator is included.

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
        """Unity-feedback closed loop :math:`CG/(1 + CG)` from reference to output.

        Parameters
        ----------
        plant : scipy.signal.TransferFunction or tuple
            Plant :math:`G(s)`, as a transfer function or ``(num, den)``
            polynomial coefficients in descending powers of :math:`s`.

        Returns
        -------
        scipy.signal.TransferFunction
            Closed-loop transfer function.
        """
        cnum, cden = self._polynomials()
        pnum, pden = num_den(plant)
        ol_num = np.polymul(cnum, pnum)
        ol_den = np.polymul(cden, pden)
        return signal.TransferFunction(ol_num, np.polyadd(ol_den, ol_num))

    def simulate(self, plant, T=None):
        """Unit-step reference response of the closed loop.

        Parameters
        ----------
        plant : scipy.signal.TransferFunction or tuple
            Plant :math:`G(s)`, as in :meth:`closed_loop`.
        T : array_like, optional
            Sample times [s]; chosen by :func:`scipy.signal.step` if omitted.

        Returns
        -------
        ResponseResult
            Time [s] and output for a unit step in the reference.
        """
        t, y = signal.step(self.closed_loop(plant), T=T)
        return ResponseResult(t, y)

    def simulate_saturated(self, plant, T, u_min, u_max, anti_windup=None, tracking_time=None, reference=1.0):
        r"""Step response with a saturating actuator, without or with anti-windup.

        The controller computes :math:`v = K_p e + z - K_d\dot y`, the
        actuator delivers :math:`u = \min(\max(v, u_{min}), u_{max})`, and
        the integrator state :math:`z` obeys

        - ``None``: :math:`\dot z = K_i e`, which winds up while
          :math:`u \ne v`;
        - ``"clamping"`` (conditional integration): :math:`\dot z = 0`
          while the actuator is saturated and :math:`K_i e` would drive it
          further, else :math:`K_i e`;
        - ``"back_calculation"``: :math:`\dot z = K_i e + (u - v)/T_t`,
          bleeding the integrator off at the tracking time constant
          :math:`T_t`.

        Parameters
        ----------
        plant : scipy.signal.TransferFunction or tuple
            Strictly proper plant :math:`G(s)`; with ``kd != 0`` its
            relative degree must be at least 2.
        T : array_like
            At least two increasing sample times [s], fine enough for a
            fixed-step RK4 integration.
        u_min, u_max : float
            Actuator limits [input units], ``u_min < u_max``.
        anti_windup : {None, "clamping", "back_calculation"}, optional
            Integrator protection.
        tracking_time : float, optional
            :math:`T_t` [s] for back-calculation; defaults to
            :math:`\sqrt{T_iT_d}` with :math:`T_i = K_p/K_i`,
            :math:`T_d = K_d/K_p`, or :math:`T_i` without derivative action.
        reference : float, optional
            Height of the reference step [output units].

        Returns
        -------
        ResponseResult
            Time [s], output and actuator output ``u``.

        Raises
        ------
        ValueError
            For an improper or biproper plant, a derivative on a plant of
            relative degree one, bad limits or an unknown anti-windup mode.

        Notes
        -----
        Derivative action acts on the measurement, so the reference step
        gives no derivative kick; within the limits the response then
        differs from :meth:`simulate` by that kick only. Without an integral
        gain the integrator is absent and ``anti_windup`` has no effect.
        Inputs are raw SI numbers with no dimensional metadata; see
        :doc:`/guides/control_systems`.

        References
        ----------
        K. J. Åström and T. Hägglund, *PID Controllers: Theory, Design, and
        Tuning*, 2nd ed. (ISA, 1995), sec. 3.5.

        Examples
        --------
        >>> from engineeringkit.control_systems import transfer_function
        >>> G = transfer_function([1], [1, 0])
        >>> r = PID(2, 1).simulate_saturated(G, np.linspace(0, 20, 2001), -0.5, 0.5, "back_calculation")
        >>> round(r.summary()["final"], 3), r.summary()["input_max"]
        (1.0, 0.5)
        """
        if anti_windup not in (None, "clamping", "back_calculation"):
            raise ValueError("anti_windup must be None, 'clamping' or 'back_calculation'")
        _finite(u_min=u_min, u_max=u_max, reference=reference)
        if not u_min < u_max:
            raise ValueError("u_min must be below u_max")
        num, den = num_den(plant)
        if len(num) >= len(den):
            raise ValueError("the plant must be strictly proper")
        if self.kd != 0 and len(den) - len(num) < 2:
            raise ValueError("derivative action needs a plant of relative degree at least 2")
        if anti_windup == "back_calculation" and self.ki != 0:
            if tracking_time is None:
                if self.kp <= 0:
                    raise ValueError("give tracking_time when kp is not positive")
                ti = self.kp / self.ki
                tracking_time = np.sqrt(ti * self.kd / self.kp) if self.kd > 0 else ti
            _positive(tracking_time=tracking_time)
        A, B, C, _ = signal.tf2ss(num, den)
        CA, n = (C @ A).ravel(), A.shape[0]

        def control(state):
            x, z = state[:n], state[n]
            y = float(C.ravel() @ x)
            e = reference - y
            v = self.kp * e + z - self.kd * float(CA @ x)
            return x, e, v, min(max(v, u_min), u_max), y

        def rhs(_, state):
            x, e, v, u, _ = control(state)
            if self.ki == 0:
                zdot = 0.0
            elif anti_windup == "back_calculation":
                zdot = self.ki * e + (u - v) / tracking_time
            elif anti_windup == "clamping" and ((v > u_max and self.ki * e > 0) or (v < u_min and self.ki * e < 0)):
                zdot = 0.0
            else:
                zdot = self.ki * e
            return np.append(A @ x + B.ravel() * u, zdot)

        sol = integrate_ode(rhs, np.zeros(n + 1), T, method="rk4")
        signals = np.array([control(state)[3:] for state in sol.history])
        return ResponseResult(sol.t, signals[:, 1], u=signals[:, 0])


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
    ku, tu = np.broadcast_arrays(np.asarray(ku, dtype=float), np.asarray(tu, dtype=float))
    kp = kp_ratio * ku
    ki = np.zeros_like(kp) if ti_ratio is None else kp / (ti_ratio * tu)
    kd = np.zeros_like(kp) if td_ratio is None else kp * td_ratio * tu
    return PID(kp[()], ki[()], kd[()])
