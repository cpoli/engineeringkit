r"""Pulse shaping for zero intersymbol interference: Nyquist's criterion and the raised cosine.

Nyquist (1928) showed that a pulse is free of intersymbol interference at
rate :math:`1/T` exactly when its spectrum, folded at multiples of
:math:`1/T`, adds to a constant. The raised-cosine family meets the
criterion with a spectrum that rolls off over a band :math:`\alpha/T` above
the minimum :math:`1/2T`; splitting it into two root-raised-cosine filters
puts a matched filter at the receiver.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["raised_cosine_pulse", "raised_cosine_spectrum", "root_raised_cosine_pulse", "raised_cosine_bandwidth"]

_SINGULAR = 1e-10


def _inputs(name, value, symbol_period, rolloff):
    _finite(**{name: value})
    _positive(symbol_period=symbol_period)
    _interval("rolloff", rolloff, 0, 1)
    return np.asarray(value, dtype=float), np.asarray(symbol_period, dtype=float), np.asarray(rolloff, dtype=float)


def _scalar(value):
    return float(value) if value.ndim == 0 else value


def raised_cosine_pulse(t, symbol_period, rolloff):
    r"""Raised-cosine pulse :math:`p(t) = \operatorname{sinc}(t/T)\,\dfrac{\cos(\pi\alpha t/T)}{1 - (2\alpha t/T)^2}`.

    Parameters
    ----------
    t : float or array_like
        Time [s].
    symbol_period : float or array_like
        Symbol period :math:`T` [s], positive.
    rolloff : float or array_like
        Roll-off (excess-bandwidth) factor :math:`\alpha` [-] in ``[0, 1]``;
        zero gives Nyquist's ideal sinc pulse.

    Returns
    -------
    float or ndarray
        Pulse value [-], one at :math:`t = 0` and zero at every other
        multiple of :math:`T`.

    Notes
    -----
    The removable singularity at :math:`|t| = T/2\alpha` is replaced by its
    limit :math:`(\pi/4)\operatorname{sinc}(1/2\alpha)`. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/communications`.

    References
    ----------
    H. Nyquist, "Certain topics in telegraph transmission theory," *Trans.
    AIEE* 47 (1928) 617-644.

    Examples
    --------
    >>> np.abs(raised_cosine_pulse([0.0, 1.0, 2.0], 1.0, 0.35)).round(12).tolist()
    [1.0, 0.0, 0.0]
    """
    t, T, a = _inputs("t", t, symbol_period, rolloff)
    x = t / T
    denominator = 1 - (2 * a * x) ** 2
    singular = np.abs(denominator) < _SINGULAR
    limit = np.where(a > 0, np.pi / 4 * np.sinc(1 / (2 * np.where(a > 0, a, 1.0))), 0.0)
    return _scalar(np.where(singular, limit, np.sinc(x) * np.cos(np.pi * a * x) / np.where(singular, 1.0, denominator)))


def raised_cosine_spectrum(f, symbol_period, rolloff):
    r"""Fourier transform of :func:`raised_cosine_pulse`: flat, then a half-cosine roll-off.

    .. math::

        P(f) = \begin{cases}
        T, & |f| \le \dfrac{1-\alpha}{2T} \\
        \dfrac{T}{2}\left[1 + \cos\dfrac{\pi T}{\alpha}\left(|f| - \dfrac{1-\alpha}{2T}\right)\right],
        & \dfrac{1-\alpha}{2T} < |f| \le \dfrac{1+\alpha}{2T} \\
        0, & \text{otherwise.}
        \end{cases}

    Parameters
    ----------
    f : float or array_like
        Frequency [Hz].
    symbol_period : float or array_like
        Symbol period :math:`T` [s], positive.
    rolloff : float or array_like
        Roll-off factor :math:`\alpha` [-] in ``[0, 1]``.

    Returns
    -------
    float or ndarray
        Spectrum [s]. Copies shifted by every multiple of :math:`1/T` add up
        to the constant :math:`T`, which is Nyquist's criterion. With
        :math:`\alpha = 0` the band edge takes the midpoint value
        :math:`T/2`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(raised_cosine_spectrum(0.25, 1.0, 0.5) + raised_cosine_spectrum(0.25 - 1.0, 1.0, 0.5))
    1.0
    """
    f, T, a = _inputs("f", f, symbol_period, rolloff)
    edge = (1 - a) / (2 * T)
    nu = np.abs(f)
    rolled = T / 2 * (1 + np.cos(np.pi * T / np.where(a > 0, a, 1.0) * (nu - edge)))
    band = np.where(a > 0, rolled, np.where(nu == edge, T / 2, 0.0))
    return _scalar(np.where(nu < edge, T, np.where(nu <= (1 + a) / (2 * T), band, 0.0)))


def root_raised_cosine_pulse(t, symbol_period, rolloff):
    r"""Unit-energy root-raised-cosine pulse, whose spectrum is :math:`\sqrt{P(f)}`.

    .. math::

        h(t) = \frac{1}{\sqrt T}\,
        \frac{\sin[\pi(1-\alpha)t/T] + 4\alpha (t/T)\cos[\pi(1+\alpha)t/T]}
        {\pi (t/T)\,[1 - (4\alpha t/T)^2]}

    Parameters
    ----------
    t : float or array_like
        Time [s].
    symbol_period : float or array_like
        Symbol period :math:`T` [s], positive.
    rolloff : float or array_like
        Roll-off factor :math:`\alpha` [-] in ``[0, 1]``.

    Returns
    -------
    float or ndarray
        Pulse value [s^-1/2]. The pulse convolved with itself is the
        raised-cosine pulse, so a root-raised-cosine transmit filter and the
        same filter as the receiver's matched filter give zero intersymbol
        interference together.

    Notes
    -----
    The singular points :math:`t = 0` and :math:`|t| = T/4\alpha` take
    their limits. Inputs are raw SI numbers with no dimensional metadata;
    see :doc:`/guides/communications`.

    Examples
    --------
    >>> round(float(root_raised_cosine_pulse(0.0, 1.0, 0.0)), 12)
    1.0
    """
    t, T, a = _inputs("t", t, symbol_period, rolloff)
    x = t / T
    at_zero = np.abs(x) < _SINGULAR
    at_quarter = np.abs(1 - (4 * a * x) ** 2) < _SINGULAR
    regular = ~(at_zero | at_quarter)
    xr = np.where(regular, x, 0.1)
    h = (np.sin(np.pi * (1 - a) * xr) + 4 * a * xr * np.cos(np.pi * (1 + a) * xr)) / (np.pi * xr * (1 - (4 * a * xr) ** 2))
    angle = np.pi / (4 * np.where(a > 0, a, 1.0))
    quarter = (a / np.sqrt(2)) * ((1 + 2 / np.pi) * np.sin(angle) + (1 - 2 / np.pi) * np.cos(angle))
    return _scalar(np.where(at_zero, 1 - a + 4 * a / np.pi, np.where(at_quarter, quarter, h)) / np.sqrt(T))


def raised_cosine_bandwidth(symbol_rate, rolloff):
    r"""Occupied (absolute) bandwidth of raised-cosine signalling, :math:`B = (1 + \alpha)R_s/2`.

    Parameters
    ----------
    symbol_rate : float or array_like
        Symbol rate :math:`R_s = 1/T` [baud], nonnegative.
    rolloff : float or array_like
        Roll-off factor :math:`\alpha` [-] in ``[0, 1]``.

    Returns
    -------
    float or ndarray
        Baseband bandwidth [Hz]; a passband (double-sideband) signal
        occupies twice this.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(raised_cosine_bandwidth(2400.0, 0.25))
    1500.0
    """
    _nonnegative(symbol_rate=symbol_rate)
    _interval("rolloff", rolloff, 0, 1)
    return (1 + np.asarray(rolloff, dtype=float)) * np.asarray(symbol_rate, dtype=float) / 2
