r"""Sensitivity and complementary sensitivity of a feedback loop.

Bode (1945) defined the sensitivity of a closed loop to a change in the
plant as :math:`S = 1/(1 + L)`: a fractional change :math:`dG/G` reaches
the closed-loop gain reduced to :math:`S\,dG/G`. :math:`S` is also the map
from an output disturbance to the output, and the complementary
sensitivity :math:`T = L/(1 + L)` the map from the reference and from
sensor noise. Since :math:`S + T = 1`, they cannot both be small at one
frequency. :math:`1/|S|` is the distance of :math:`L(j\omega)` from
:math:`-1`, so its peak bounds both stability margins.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import positive as _positive
from engineeringkit.control_systems.systems.frequency import nyquist
from engineeringkit.control_systems.systems.transfer import transfer_function
from engineeringkit.control_systems.utils.polynomials import num_den

__all__ = ["sensitivity", "complementary_sensitivity", "sensitivity_peak"]


def sensitivity(loop):
    r"""Sensitivity function :math:`S = 1/(1 + L)` of a unity-feedback loop.

    Parameters
    ----------
    loop : scipy.signal.lti
        SISO open-loop model :math:`L(s) = N(s)/D(s)`.

    Returns
    -------
    scipy.signal.TransferFunction
        :math:`D/(D + N)`, from an output disturbance to the output and from
        the reference to the error.

    Notes
    -----
    Pole-zero cancellations are not removed. See
    :doc:`/guides/control_systems`.

    References
    ----------
    H. W. Bode, *Network Analysis and Feedback Amplifier Design* (Van
    Nostrand, 1945).

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> S = sensitivity(transfer_function([1], [1, 0]))
    >>> S.num.tolist(), S.den.tolist()
    ([1.0, 0.0], [1.0, 1.0])
    """
    num, den = num_den(loop)
    return transfer_function(den, np.polyadd(den, num))


def complementary_sensitivity(loop):
    r"""Complementary sensitivity :math:`T = L/(1 + L)` of a unity-feedback loop.

    Parameters
    ----------
    loop : scipy.signal.lti
        SISO open-loop model :math:`L(s) = N(s)/D(s)`.

    Returns
    -------
    scipy.signal.TransferFunction
        :math:`N/(D + N)`, from the reference to the output and, with a sign
        change, from sensor noise to the output.

    Notes
    -----
    :math:`S + T = 1` identically. See :doc:`/guides/control_systems`.

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> T = complementary_sensitivity(transfer_function([1], [1, 0]))
    >>> T.num.tolist(), T.den.tolist()
    ([1.0], [1.0, 1.0])
    """
    num, den = num_den(loop)
    return transfer_function(num, np.polyadd(den, num))


def sensitivity_peak(loop, w=None):
    r"""Peaks of :math:`|S|` and :math:`|T|` and the stability margins they guarantee.

    Parameters
    ----------
    loop : scipy.signal.lti
        SISO open-loop model :math:`L(s)` whose closed loop is stable.
    w : array_like, optional
        Positive angular frequencies [rad/s]; defaults to 20000 points from
        :math:`10^{-4}` to :math:`10^4`.

    Returns
    -------
    dict
        ``Ms`` and ``Mt`` [-], the peaks of :math:`|S(j\omega)|` and
        :math:`|T(j\omega)|`; ``Ms_frequency`` and ``Mt_frequency`` [rad/s],
        where they occur; ``gain_margin_bound`` :math:`M_s/(M_s - 1)` [-]
        and ``phase_margin_bound_deg`` :math:`2\arcsin(1/(2M_s))` [deg],
        lower bounds on the margins.

    Notes
    -----
    :math:`L(j\omega)` stays outside the circle of radius :math:`1/M_s`
    about :math:`-1`, which gives both bounds; :math:`M_s` between 1.3 and 2
    is the usual design range. Peaks are read on the grid, and they mean
    nothing for an unstable closed loop. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/control_systems`.

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> p = sensitivity_peak(transfer_function([1], [1, 0]))
    >>> round(p["Ms"], 6), round(p["Mt"], 6)
    (1.0, 1.0)
    """
    w = np.logspace(-4, 4, 20000) if w is None else np.sort(np.asarray(w, float))
    _positive(w=w)
    h = nyquist(loop, w)
    s, t = np.abs(1 / (1 + h)), np.abs(h / (1 + h))
    i, j = int(np.argmax(s)), int(np.argmax(t))
    ms = float(s[i])
    return {
        "Ms": ms,
        "Ms_frequency": float(w[i]),
        "Mt": float(t[j]),
        "Mt_frequency": float(w[j]),
        "gain_margin_bound": np.inf if ms <= 1 else ms / (ms - 1),
        "phase_margin_bound_deg": float(np.degrees(2 * np.arcsin(min(1.0, 1 / (2 * ms))))),
    }
