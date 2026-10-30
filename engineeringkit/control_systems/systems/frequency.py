r"""Frequency response: Bode plots, the Nyquist criterion, stability margins and the Nichols chart.

Evaluating :math:`L(j\omega)` along the imaginary axis tells how a loop
amplifies and delays sinusoids. Nyquist's criterion counts encirclements of
:math:`-1` to find closed-loop instability; Bode's gain and phase margins
measure how far the loop is from it. Nichols plotted the open-loop gain in
decibels against its phase, where the closed-loop gain of a unity-feedback
loop is a fixed grid of contours.
"""

from __future__ import annotations

import numpy as np
from scipy import signal

from engineeringkit._validation import positive as _positive
from engineeringkit.control_systems.utils.polynomials import num_den

__all__ = ["bode", "nyquist", "nyquist_encirclements", "stability_margins", "nichols", "nichols_closed_loop"]


def bode(sys, w=None):
    """Bode magnitude and phase of a continuous-time model.

    Parameters
    ----------
    sys : scipy.signal.lti
        Continuous-time model.
    w : array_like, optional
        Angular frequencies [rad/s]; chosen automatically if omitted.

    Returns
    -------
    tuple of ndarray
        ``(w, magnitude, phase)``: angular frequency [rad/s], gain [dB] and
        phase [deg].

    Notes
    -----
    Wraps :func:`scipy.signal.bode`. See :doc:`/guides/control_systems`.

    References
    ----------
    H. W. Bode, *Network Analysis and Feedback Amplifier Design* (Van
    Nostrand, 1945).

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> w, mag, phase = bode(transfer_function([1], [1, 1]), [1.0])
    >>> round(float(mag[0]), 3), round(float(phase[0]), 1)
    (-3.01, -45.0)
    """
    if w is not None:
        _positive(w=w)
    return signal.bode(sys, w=w)


def nyquist(sys, w):
    r"""Complex frequency response :math:`G(j\omega)`.

    Parameters
    ----------
    sys : scipy.signal.lti
        SISO continuous-time model.
    w : float or array_like
        Angular frequencies [rad/s]; negative values give the mirror image.

    Returns
    -------
    complex or ndarray of complex
        Output/input gain at each frequency.

    Notes
    -----
    See :doc:`/guides/control_systems`.

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> complex(np.round(nyquist(transfer_function([2, 2], [1, 3, 2]), 2.0), 12))
    (0.5-0.5j)
    """
    num, den = num_den(sys)
    w = np.asarray(w, float)
    if not np.all(np.isfinite(w)):
        raise ValueError("w must contain finite values")
    s = 1j * w
    return np.polyval(num, s) / np.polyval(den, s)


def nyquist_encirclements(sys, w=None):
    r"""Net clockwise encirclements of :math:`-1` by the Nyquist plot of a loop.

    Nyquist's criterion: the closed loop :math:`L/(1+L)` has
    :math:`Z = N + P` right-half-plane poles, where :math:`P` is the number of
    open-loop right-half-plane poles and :math:`N` the value returned here.

    Parameters
    ----------
    sys : scipy.signal.lti
        Proper open-loop model :math:`L(s)` with no poles on the imaginary axis.
    w : array_like, optional
        Positive angular frequencies [rad/s] spanning all features of
        :math:`L`; defaults to 40001 points from :math:`10^{-4}` to
        :math:`10^4`.

    Returns
    -------
    int
        Clockwise encirclements of :math:`-1` (negative if counter-clockwise).

    Raises
    ------
    ValueError
        If :math:`L` is improper, has imaginary-axis poles, or passes through
        :math:`-1` on the grid.

    Notes
    -----
    The plot is traced for :math:`\omega` from :math:`-\infty` to
    :math:`\infty`; for a proper loop the infinite semicircle maps to a single
    point. A grid that misses a fast phase swing can miscount. See
    :doc:`/guides/control_systems`.

    References
    ----------
    H. Nyquist, "Regeneration theory," *Bell Syst. Tech. J.* 11 (1932)
    126-147.

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> nyquist_encirclements(transfer_function([12], [1, 3, 3, 1]))
    2
    """
    num, den = num_den(sys)
    if len(num) > len(den):
        raise ValueError("the loop transfer function must be proper")
    roots = np.roots(den)
    if np.any(np.abs(roots.real) <= 1e-10 * max(1.0, float(np.max(np.abs(roots), initial=0.0)))):
        raise ValueError("poles on the imaginary axis need an indented contour, which is not supported")
    w = np.logspace(-4, 4, 40001) if w is None else np.sort(np.asarray(w, float))
    _positive(w=w)
    s = 1j * np.concatenate((-w[::-1], w))
    return_difference = 1 + np.polyval(num, s) / np.polyval(den, s)
    if np.min(np.abs(return_difference)) < 1e-9:
        raise ValueError("the Nyquist plot passes through -1: the closed loop is marginally stable")
    phase = np.unwrap(np.angle(return_difference))
    return int(round(-(phase[-1] - phase[0]) / (2 * np.pi)))


def _crossings(w, value, period=None):
    """Log-frequencies where ``value`` crosses zero, or any multiple of ``period``.

    Returns the interpolated frequencies (linear in log-frequency between
    grid points) and the index of the grid segment holding each crossing.
    """
    k = np.floor(value / period) if period else (value >= 0).astype(float)
    i = np.flatnonzero(np.diff(k) != 0)
    target = np.maximum(k[i], k[i + 1]) * period if period else np.zeros(i.size)
    frac = (target - value[i]) / (value[i + 1] - value[i])
    return np.exp(np.log(w[i]) + frac * (np.log(w[i + 1]) - np.log(w[i]))), i, frac


def stability_margins(sys, w=None):
    r"""Worst gain and phase margins over all crossovers of a loop's frequency response.

    Parameters
    ----------
    sys : scipy.signal.lti
        SISO open-loop model :math:`L(s)`.
    w : array_like, optional
        Positive angular frequencies [rad/s]; defaults to 20000 points from
        :math:`10^{-4}` to :math:`10^4`.

    Returns
    -------
    dict
        ``gain_margin`` [-], the margin at the phase crossover nearest
        instability (``inf`` without a phase crossover);
        ``phase_margin_deg`` [deg], the margin of smallest magnitude over
        the gain crossovers (``nan`` without one); ``gain_cross_frequency`` and
        ``phase_cross_frequency`` [rad/s], the crossovers at which those worst
        margins are read (``nan`` if absent); and the arrays
        ``phase_cross_frequencies`` [rad/s] with their ``gain_margins`` [-],
        and ``gain_cross_frequencies`` [rad/s] with their
        ``phase_margins_deg`` [deg], over every crossover in increasing
        frequency.

    Notes
    -----
    A phase crossover is any frequency where the phase passes
    :math:`-180^\circ + 360^\circ k`. Its gain margin
    :math:`1/|L(j\omega)|` is below one when *reducing* the gain by that
    factor reaches instability, as in a conditionally stable loop; the worst
    gain margin is the one with the smallest :math:`|\log \mathrm{GM}|`.
    Phase margins are wrapped into :math:`[-180^\circ, 180^\circ)`; a
    negative one means the plot crosses the unit circle on the far side of
    :math:`-1`, and a phase *lead* of that size reaches instability. With one
    crossover and a stable open loop that signals an unstable closed loop,
    but an open-loop unstable loop that needs counter-clockwise encirclements
    can cross on both sides, so the worst phase margin is the one of
    smallest magnitude.
    Crossings are interpolated linearly in log-frequency between grid points,
    so a grid that misses a crossover misses its margin. Time delays need a
    rational approximation such as :func:`pade_delay`. See
    :doc:`/guides/control_systems`.

    References
    ----------
    H. W. Bode, "Relations between attenuation and phase in feedback
    amplifier design," *Bell Syst. Tech. J.* 19 (1940) 421-454.

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> m = stability_margins(transfer_function([2], [1, 3, 2, 0]))
    >>> round(m["gain_margin"], 2), round(m["phase_cross_frequency"], 3)
    (3.0, 1.414)
    """
    w = np.logspace(-4, 4, 20000) if w is None else np.sort(np.asarray(w, float))
    _positive(w=w)
    h = nyquist(sys, w)
    log_mag = np.log(np.abs(h))
    phase = np.unwrap(np.angle(h)) * 180 / np.pi
    wg, ig, fg = _crossings(w, log_mag)
    wp, ip, fp = _crossings(w, phase + 180, 360.0)
    pm = (phase[ig] + fg * (phase[ig + 1] - phase[ig]) + 360) % 360 - 180
    gm = np.exp(-(log_mag[ip] + fp * (log_mag[ip + 1] - log_mag[ip])))
    jg = int(np.argmin(np.abs(pm))) if pm.size else None
    jp = int(np.argmin(np.abs(np.log(gm)))) if gm.size else None
    return {
        "gain_margin": np.inf if jp is None else float(gm[jp]),
        "phase_margin_deg": np.nan if jg is None else float(pm[jg]),
        "gain_cross_frequency": np.nan if jg is None else float(wg[jg]),
        "phase_cross_frequency": np.nan if jp is None else float(wp[jp]),
        "gain_margins": gm,
        "phase_cross_frequencies": wp,
        "phase_margins_deg": pm,
        "gain_cross_frequencies": wg,
    }


def nichols_closed_loop(gain_db, phase_deg):
    r"""Closed-loop gain and phase of a unity-feedback loop from its open-loop Nichols coordinates.

    Parameters
    ----------
    gain_db : float or array_like
        Open-loop gain :math:`20\log_{10}|L|` [dB].
    phase_deg : float or array_like
        Open-loop phase :math:`\angle L` [deg].

    Returns
    -------
    tuple of ndarray
        ``(closed_gain_db, closed_phase_deg)`` of :math:`T = L/(1 + L)`, the
        phase in :math:`(-180^\circ, 180^\circ]`. The gain grows without
        bound towards the critical point :math:`L = -1`.

    Notes
    -----
    The contours of constant :math:`|T|` drawn in these coordinates are the
    M-contours of the Nichols chart; they repeat every :math:`360^\circ` of
    open-loop phase. Inputs are raw SI numbers with no dimensional metadata;
    see :doc:`/guides/control_systems`.

    References
    ----------
    H. M. James, N. B. Nichols and R. S. Phillips, *Theory of
    Servomechanisms*, MIT Radiation Laboratory Series 25 (McGraw-Hill, 1947).

    Examples
    --------
    >>> m, p = nichols_closed_loop(0.0, -90.0)
    >>> round(float(m), 3), round(float(p), 1)
    (-3.01, -45.0)
    """
    gain_db, phase_deg = np.asarray(gain_db, float), np.asarray(phase_deg, float)
    if not (np.all(np.isfinite(gain_db)) and np.all(np.isfinite(phase_deg))):
        raise ValueError("gain_db and phase_deg must be finite")
    L = 10 ** (gain_db / 20) * np.exp(1j * np.radians(phase_deg))
    T = L / (1 + L)
    return 20 * np.log10(np.abs(T)), np.degrees(np.angle(T))


def nichols(sys, w=None):
    r"""Open-loop Nichols locus of a loop, with the closed-loop response and resonant peak it implies.

    Parameters
    ----------
    sys : scipy.signal.lti
        SISO open-loop model :math:`L(s)`.
    w : array_like, optional
        Positive angular frequencies [rad/s]; defaults to 4000 points from
        :math:`10^{-3}` to :math:`10^3`.

    Returns
    -------
    dict
        ``w`` [rad/s]; the open-loop ``gain_db`` [dB] and continuous
        ``phase_deg`` [deg]; the closed-loop ``closed_gain_db`` [dB] and
        ``closed_phase_deg`` [deg] of :math:`L/(1 + L)`; and the resonant
        peak ``resonant_peak_db`` [dB] with its ``resonant_frequency``
        [rad/s], the largest closed-loop gain on the grid.

    Notes
    -----
    The resonant peak is the highest M-contour the locus touches. The
    open-loop phase is unwrapped from the lowest frequency, so it can leave
    :math:`[-360^\circ, 0^\circ]`. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/control_systems`.

    References
    ----------
    H. M. James, N. B. Nichols and R. S. Phillips, *Theory of
    Servomechanisms*, MIT Radiation Laboratory Series 25 (McGraw-Hill, 1947).

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> n = nichols(transfer_function([1], [1, 1, 0]))
    >>> round(n["resonant_peak_db"], 2), round(n["resonant_frequency"], 3)
    (1.25, 0.707)
    """
    w = np.logspace(-3, 3, 4000) if w is None else np.sort(np.asarray(w, float))
    _positive(w=w)
    h = nyquist(sys, w)
    gain_db = 20 * np.log10(np.abs(h))
    phase_deg = np.degrees(np.unwrap(np.angle(h)))
    T = h / (1 + h)
    closed_gain_db = 20 * np.log10(np.abs(T))
    i = int(np.argmax(closed_gain_db))
    return {
        "w": w,
        "gain_db": gain_db,
        "phase_deg": phase_deg,
        "closed_gain_db": closed_gain_db,
        "closed_phase_deg": np.degrees(np.unwrap(np.angle(T))),
        "resonant_peak_db": float(closed_gain_db[i]),
        "resonant_frequency": float(w[i]),
    }
