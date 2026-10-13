r"""Transient specifications of the standard underdamped second-order system.

Hazen's servomechanism theory (1934) characterised a position servo by its
damping ratio :math:`\zeta` and natural frequency :math:`\omega_n`, with
closed-loop transfer function
:math:`\omega_n^2/(s^2 + 2\zeta\omega_n s + \omega_n^2)`.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive

__all__ = ["second_order_metrics"]


def second_order_metrics(zeta, omega_n):
    r"""Peak time, percent overshoot and 2 % settling time of an underdamped second-order step response.

    Parameters
    ----------
    zeta : float or array_like
        Damping ratio [-], strictly between 0 and 1.
    omega_n : float or array_like
        Undamped natural frequency [rad/s], positive.

    Returns
    -------
    dict
        ``peak_time`` :math:`\pi/\omega_d` [s], ``percent_overshoot``
        :math:`100\,e^{-\pi\zeta/\sqrt{1-\zeta^2}}` [%] and
        ``settling_time_2pct`` :math:`4/(\zeta\omega_n)` [s], with
        :math:`\omega_d = \omega_n\sqrt{1-\zeta^2}`.

    Notes
    -----
    Peak time and overshoot are exact; the settling time is the usual
    envelope estimate, not the exact 2 % crossing. See
    :doc:`/guides/control_systems`.

    References
    ----------
    H. L. Hazen, "Theory of servo-mechanisms," *J. Franklin Inst.* 218 (1934)
    279-331.

    Examples
    --------
    >>> m = second_order_metrics(0.5, 2.0)
    >>> round(float(m["percent_overshoot"]), 2), float(m["settling_time_2pct"])
    (16.3, 4.0)
    """
    _interval("zeta", zeta, 0, 1, closed=False)
    _positive(omega_n=omega_n)
    zeta, omega_n = np.asarray(zeta, float), np.asarray(omega_n, float)
    wd = omega_n * np.sqrt(1 - zeta * zeta)
    return {
        "peak_time": np.pi / wd,
        "percent_overshoot": 100 * np.exp(-np.pi * zeta / np.sqrt(1 - zeta * zeta)),
        "settling_time_2pct": 4 / (zeta * omega_n),
    }
