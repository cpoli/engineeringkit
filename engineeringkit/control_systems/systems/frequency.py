r"""Frequency response: Bode plots, the Nyquist criterion and stability margins.

Evaluating :math:`L(j\omega)` along the imaginary axis tells how a loop
amplifies and delays sinusoids. Nyquist's criterion counts encirclements of
:math:`-1` to find closed-loop instability; Bode's gain and phase margins
measure how far the loop is from it.
"""

from __future__ import annotations

import numpy as np
from scipy import signal

from engineeringkit._validation import positive as _positive
from engineeringkit.control_systems.utils.polynomials import num_den

__all__ = ["bode", "nyquist", "nyquist_encirclements", "stability_margins"]


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


def _crossing(w, value, target):
    """Log-frequency of the first crossing of ``value`` through ``target``, by linear interpolation."""
    i = np.flatnonzero(np.diff(np.sign(value - target)) != 0)
    if not i.size:
        return np.nan, None
    i = i[0]
    frac = (target - value[i]) / (value[i + 1] - value[i])
    return float(np.exp(np.log(w[i]) + frac * (np.log(w[i + 1]) - np.log(w[i])))), i


def stability_margins(sys, w=None):
    """Gain and phase margins from the first crossovers of a loop's frequency response.

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
        ``gain_margin`` [-] (``inf`` without a phase crossover),
        ``phase_margin_deg`` [deg] (``nan`` without a gain crossover),
        ``gain_cross_frequency`` and ``phase_cross_frequency`` [rad/s].

    Notes
    -----
    Crossings are interpolated linearly in log-frequency between grid
    points. Only the first crossover of each kind is reported, so loops with
    multiple crossovers, time delays or negative gain need separate analysis.
    See :doc:`/guides/control_systems`.

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
    w = np.logspace(-4, 4, 20000) if w is None else np.asarray(w, float)
    _positive(w=w)
    h = nyquist(sys, w)
    mag = np.abs(h)
    phase = np.unwrap(np.angle(h)) * 180 / np.pi
    wg, ig = _crossing(w, np.log(mag), 0.0)
    wp, ip = _crossing(w, phase, -180.0)
    pm = np.nan if ig is None else 180 + float(np.interp(np.log(wg), np.log(w[ig : ig + 2]), phase[ig : ig + 2]))
    gm = np.inf if ip is None else 1 / float(np.exp(np.interp(np.log(wp), np.log(w[ip : ip + 2]), np.log(mag[ip : ip + 2]))))
    return {"gain_margin": float(gm), "phase_margin_deg": float(pm), "gain_cross_frequency": float(wg), "phase_cross_frequency": float(wp)}
