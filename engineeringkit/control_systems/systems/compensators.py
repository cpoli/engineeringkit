r"""Lead, lag and lead-lag compensators designed from the loop's frequency response.

A first-order compensator :math:`C(s) = (Ts + 1)/(\alpha Ts + 1)` has unit
gain at :math:`s = 0`, so it leaves the error constants alone. With
:math:`\alpha < 1` (lead) it adds up to
:math:`\phi_m = \arcsin\frac{1-\alpha}{1+\alpha}` of phase at
:math:`\omega_m = 1/(T\sqrt\alpha)`, where its gain is :math:`1/\sqrt\alpha`;
with :math:`\alpha = \beta > 1` (lag) it cuts the high-frequency gain by
:math:`\beta` and moves the crossover down to where the loop has phase to
spare. These are the textbook procedures of Ogata (*Modern Control
Engineering*, ch. 7) and Franklin, Powell and Emami-Naeini (ch. 6).
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive
from engineeringkit.control_systems.core.base import CompensatorResult
from engineeringkit.control_systems.systems.frequency import _crossings, nyquist, stability_margins
from engineeringkit.control_systems.systems.transfer import transfer_function
from engineeringkit.control_systems.utils.polynomials import num_den

__all__ = ["lead_compensator", "lag_compensator", "lead_lag_compensator", "design_lead", "design_lag", "design_lead_lag"]


def _series(*systems):
    num, den = np.array([1.0]), np.array([1.0])
    for sys in systems:
        n, d = num_den(sys)
        num, den = np.polymul(num, n), np.polymul(den, d)
    return transfer_function(num, den)


def _lead_alpha(phase_lead_deg):
    sin_phi = np.sin(np.radians(phase_lead_deg))
    return (1 - sin_phi) / (1 + sin_phi)


def lead_compensator(phase_lead_deg, frequency):
    r"""Phase-lead compensator :math:`(Ts + 1)/(\alpha Ts + 1)` with its peak lead at a given frequency.

    Parameters
    ----------
    phase_lead_deg : float
        Maximum phase lead :math:`\phi_m` [deg], strictly between 0 and 90.
    frequency : float
        Frequency of the maximum lead, :math:`\omega_m` [rad/s], positive.

    Returns
    -------
    scipy.signal.TransferFunction
        :math:`C(s)` with :math:`\alpha = (1 - \sin\phi_m)/(1 + \sin\phi_m)`
        and :math:`T = 1/(\omega_m\sqrt\alpha)`: unit gain at DC, gain
        :math:`1/\sqrt\alpha` at :math:`\omega_m` and :math:`1/\alpha` at
        high frequency.

    Notes
    -----
    The zero :math:`1/T` and pole :math:`1/(\alpha T)` sit symmetrically
    about :math:`\omega_m` on a log scale. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/control_systems`.

    Examples
    --------
    >>> C = lead_compensator(30.0, 2.0)
    >>> [round(float(c), 4) for c in C.num / C.den[-1]], [round(float(c), 4) for c in C.den / C.den[-1]]
    ([0.866, 1.0], [0.2887, 1.0])
    """
    _interval("phase_lead_deg", phase_lead_deg, 0, 90, closed=False)
    _positive(frequency=frequency)
    alpha = _lead_alpha(phase_lead_deg)
    T = 1 / (frequency * np.sqrt(alpha))
    return transfer_function([T, 1.0], [alpha * T, 1.0])


def lag_compensator(beta, zero_frequency):
    r"""Phase-lag compensator :math:`(Ts + 1)/(\beta Ts + 1)` that cuts the high-frequency gain by :math:`\beta`.

    Parameters
    ----------
    beta : float
        Pole-zero ratio :math:`\beta` [-], greater than 1.
    zero_frequency : float
        Zero :math:`1/T` [rad/s], positive; the pole is at
        :math:`1/(\beta T)`.

    Returns
    -------
    scipy.signal.TransferFunction
        :math:`C(s)` with unit gain at DC and :math:`1/\beta` at high
        frequency.

    Notes
    -----
    Placed a decade or so below the crossover, the lag lowers the loop gain
    there by :math:`\beta` while adding only a few degrees of lag.
    Equivalently, it raises the low-frequency gain, and so the error
    constants, by :math:`\beta` relative to the crossover gain. Inputs are
    raw SI numbers with no dimensional metadata; see
    :doc:`/guides/control_systems`.

    Examples
    --------
    >>> C = lag_compensator(10.0, 0.1)
    >>> C.num.tolist(), C.den.tolist()
    ([0.1, 0.01], [1.0, 0.01])
    """
    _positive(beta=beta, zero_frequency=zero_frequency)
    if beta <= 1:
        raise ValueError("beta must be greater than 1")
    T = 1 / zero_frequency
    return transfer_function([T, 1.0], [beta * T, 1.0])


def lead_lag_compensator(phase_lead_deg, lead_frequency, beta, lag_zero_frequency):
    r"""Series lead-lag compensator, the product of :func:`lead_compensator` and :func:`lag_compensator`.

    Parameters
    ----------
    phase_lead_deg : float
        Maximum phase lead of the lead part [deg], strictly between 0 and 90.
    lead_frequency : float
        Frequency of that maximum lead [rad/s], positive.
    beta : float
        Pole-zero ratio of the lag part [-], greater than 1.
    lag_zero_frequency : float
        Zero of the lag part [rad/s], positive.

    Returns
    -------
    scipy.signal.TransferFunction
        :math:`C(s)` with unit gain at DC.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/control_systems`.

    Examples
    --------
    >>> C = lead_lag_compensator(30.0, 2.0, 10.0, 0.1)
    >>> len(C.den) - 1
    2
    """
    return _series(lead_compensator(phase_lead_deg, lead_frequency), lag_compensator(beta, lag_zero_frequency))


def _result(kind, compensator, loop, w):
    compensated = _series(compensator, loop)
    m = stability_margins(compensated, w)
    return CompensatorResult(kind, compensator, loop, compensated, m["phase_margin_deg"], m["gain_cross_frequency"], m["gain_margin"])


def _grid(w):
    w = np.logspace(-4, 4, 20000) if w is None else np.sort(np.asarray(w, float))
    _positive(w=w)
    return w


def _check_design_inputs(phase_margin_deg, extra_deg):
    _interval("phase_margin_deg", phase_margin_deg, 0, 90, closed=False)
    _interval("extra_deg", extra_deg, 0, 45)


def design_lead(loop, phase_margin_deg, extra_deg=5.0, w=None):
    r"""Lead compensator that raises a loop's phase margin, by Ogata's frequency-response procedure.

    1. Read the present phase margin :math:`\gamma_0` of :math:`L`.
    2. The lead must supply :math:`\phi_m = \gamma - \gamma_0 + \delta`,
       where :math:`\delta` (``extra_deg``) pays for the extra lag at the
       crossover, which moves up.
    3. Put the peak lead at the frequency where :math:`|L| = \sqrt\alpha`,
       so that the compensator's gain :math:`1/\sqrt\alpha` makes it the new
       crossover.

    Parameters
    ----------
    loop : scipy.signal.lti
        Open-loop model :math:`L(s)`, with its gain already set for the
        required error constant.
    phase_margin_deg : float
        Required phase margin :math:`\gamma` [deg], strictly between 0 and 90.
    extra_deg : float, optional
        Safety allowance :math:`\delta` [deg], between 0 and 45.
    w : array_like, optional
        Positive angular frequencies [rad/s] for reading the response;
        defaults to 20000 points from :math:`10^{-4}` to :math:`10^4`.

    Returns
    -------
    CompensatorResult
        The compensator, both loops and the margins actually achieved.

    Raises
    ------
    ValueError
        If the loop has no gain crossover, already meets the margin, or
        needs 90° or more of lead from one stage.

    Notes
    -----
    A single lead stage gives less than 90° and, in practice, about 60°;
    beyond that use two stages or :func:`design_lead_lag`. The achieved
    margin differs from :math:`\gamma` by however far :math:`\delta`
    misjudges the extra lag. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/control_systems`.

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> r = design_lead(transfer_function([40], [1, 2, 0]), 50.0)
    >>> round(r.phase_margin_deg, 1), round(r.gain_cross_frequency, 2)
    (49.8, 8.85)
    """
    _check_design_inputs(phase_margin_deg, extra_deg)
    w = _grid(w)
    pm0 = stability_margins(loop, w)["phase_margin_deg"]
    if np.isnan(pm0):
        raise ValueError("the loop has no gain crossover on the frequency grid")
    lead = phase_margin_deg - pm0 + extra_deg
    if lead <= 0:
        raise ValueError(f"the loop already has a {pm0:.1f} deg phase margin")
    if lead >= 90:
        raise ValueError(f"one lead stage cannot supply {lead:.1f} deg; use two stages or design_lead_lag")
    alpha = _lead_alpha(lead)
    log_mag = np.log(np.abs(nyquist(loop, w)))
    wm, _, _ = _crossings(w, log_mag - 0.5 * np.log(alpha))
    return _result("lead", lead_compensator(lead, float(wm[-1])), loop, w)


def design_lag(loop, phase_margin_deg, extra_deg=5.0, zero_ratio=10.0, w=None):
    r"""Lag compensator that lowers a loop's crossover to where it has the required phase margin.

    1. Find the frequency :math:`\omega_c` where :math:`\angle L =
       -180^\circ + \gamma + \delta`, the allowance :math:`\delta` paying
       for the lag's residual phase there.
    2. Choose :math:`\beta = |L(j\omega_c)|`, so that the lag's
       high-frequency gain :math:`1/\beta` makes :math:`\omega_c` the new
       crossover.
    3. Put the lag zero at :math:`\omega_c/` ``zero_ratio``.

    Parameters
    ----------
    loop : scipy.signal.lti
        Open-loop model :math:`L(s)`, with its gain already set for the
        required error constant.
    phase_margin_deg : float
        Required phase margin :math:`\gamma` [deg], strictly between 0 and 90.
    extra_deg : float, optional
        Allowance :math:`\delta` for the lag's phase at :math:`\omega_c`
        [deg], between 0 and 45.
    zero_ratio : float, optional
        Ratio of :math:`\omega_c` to the lag zero [-], greater than 1.
    w : array_like, optional
        Positive angular frequencies [rad/s] for reading the response.

    Returns
    -------
    CompensatorResult
        The compensator, both loops and the margins actually achieved.

    Raises
    ------
    ValueError
        If the phase never reaches the target, or the gain there is already
        below one (a lag cannot raise it).

    Notes
    -----
    A zero a decade below :math:`\omega_c` costs about
    :math:`\arctan 0.1 \approx 5.7^\circ` there, which is why
    :math:`\delta = 5^\circ` is the usual allowance. The price is a slower
    response: the crossover, and so the bandwidth, goes down. Inputs are raw
    SI numbers with no dimensional metadata; see
    :doc:`/guides/control_systems`.

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> r = design_lag(transfer_function([5], [0.5, 1.5, 1, 0]), 40.0)
    >>> round(r.phase_margin_deg, 1), round(r.gain_cross_frequency, 3)
    (39.9, 0.564)
    """
    _check_design_inputs(phase_margin_deg, extra_deg)
    _positive(zero_ratio=zero_ratio)
    if zero_ratio <= 1:
        raise ValueError("zero_ratio must be greater than 1")
    w = _grid(w)
    h = nyquist(loop, w)
    phase = np.unwrap(np.angle(h)) * 180 / np.pi
    wc, _, _ = _crossings(w, phase + 180 - phase_margin_deg - extra_deg, 360.0)
    if not wc.size:
        raise ValueError("the loop phase never reaches the target on the frequency grid")
    wc = float(wc[0])
    beta = float(np.abs(nyquist(loop, wc)))
    if beta <= 1:
        raise ValueError("the loop gain is already below one at the target crossover; a lag cannot raise it")
    return _result("lag", lag_compensator(beta, wc / zero_ratio), loop, w)


def design_lead_lag(loop, phase_margin_deg, crossover_frequency, extra_deg=5.0, zero_ratio=10.0, w=None):
    r"""Lead-lag compensator giving a chosen crossover frequency and phase margin.

    The lead puts its peak at the chosen crossover :math:`\omega_c` and
    supplies :math:`\phi_m = \gamma + \delta - 180^\circ - \angle
    L(j\omega_c)`; the lag then sets :math:`\beta = |L(j\omega_c)|/\sqrt\alpha`
    so that the total gain there is one.

    Parameters
    ----------
    loop : scipy.signal.lti
        Open-loop model :math:`L(s)`, with its gain already set for the
        required error constant.
    phase_margin_deg : float
        Required phase margin :math:`\gamma` [deg], strictly between 0 and 90.
    crossover_frequency : float
        Required gain-crossover frequency :math:`\omega_c` [rad/s], positive.
    extra_deg : float, optional
        Allowance :math:`\delta` for the lag's phase at :math:`\omega_c`
        [deg], between 0 and 45.
    zero_ratio : float, optional
        Ratio of :math:`\omega_c` to the lag zero [-], greater than 1.
    w : array_like, optional
        Positive angular frequencies [rad/s] for reading the achieved
        margins.

    Returns
    -------
    CompensatorResult
        The compensator, both loops and the margins actually achieved.

    Raises
    ------
    ValueError
        If the required lead is not between 0 and 90°, or the loop gain at
        :math:`\omega_c` is too low for a lag to be needed (then use
        :func:`design_lead`).

    Notes
    -----
    The lead buys phase at the crossover, the lag gives back the excess gain
    the lead adds, so the error constants are unchanged and the crossover
    lands where it was asked. The phase of :math:`L` is read in
    :math:`[-360^\circ, 0^\circ)`. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/control_systems`.

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> r = design_lead_lag(transfer_function([20], [1, 3, 2, 0]), 45.0, 1.0)
    >>> round(r.phase_margin_deg, 1), round(r.gain_cross_frequency, 3)
    (44.6, 1.004)
    """
    _check_design_inputs(phase_margin_deg, extra_deg)
    _positive(crossover_frequency=crossover_frequency, zero_ratio=zero_ratio)
    if zero_ratio <= 1:
        raise ValueError("zero_ratio must be greater than 1")
    h = complex(nyquist(loop, crossover_frequency))
    angle = np.degrees(np.angle(h)) % 360 - 360
    lead = phase_margin_deg + extra_deg - 180 - angle
    if not 0 < lead < 90:
        raise ValueError(f"the crossover needs {lead:.1f} deg of lead; one stage gives between 0 and 90")
    beta = abs(h) / np.sqrt(_lead_alpha(lead))
    if beta <= 1:
        raise ValueError("the loop gain at the crossover is too low for a lag; use design_lead")
    C = lead_lag_compensator(lead, crossover_frequency, beta, crossover_frequency / zero_ratio)
    return _result("lead-lag", C, loop, _grid(w))
