r"""Transfer functions, feedback interconnection and time responses.

A SISO linear time-invariant system is described by
:math:`G(s) = N(s)/D(s)`. Its poles (roots of :math:`D`) set the natural
modes and its zeros (roots of :math:`N`) shape how inputs excite them.
Black's negative-feedback loop gives the closed-loop map
:math:`G/(1 + GH)`.
"""

from __future__ import annotations

import numpy as np
from scipy import signal

from engineeringkit._validation import finite as _finite
from engineeringkit.control_systems.utils.polynomials import as_polynomial, num_den

__all__ = ["transfer_function", "feedback", "poles", "zeros", "dc_gain", "step_response", "impulse_response"]


def transfer_function(num, den):
    """Continuous-time transfer function :math:`N(s)/D(s)`.

    Parameters
    ----------
    num, den : array_like
        Numerator and denominator coefficients in descending powers of
        :math:`s`. The denominator must have degree at least zero and must not
        be identically zero.

    Returns
    -------
    scipy.signal.TransferFunction
        Continuous-time SISO model.

    Notes
    -----
    Units of :math:`G` are output units per input unit. See
    :doc:`/guides/control_systems`.

    Examples
    --------
    >>> G = transfer_function([1], [1, 1])
    >>> G.den.tolist()
    [1.0, 1.0]
    """
    return signal.TransferFunction(as_polynomial(num, "num"), as_polynomial(den, "den"))


def feedback(G, H=1.0, sign=-1):
    r"""Closed-loop transfer function of a feedback loop, :math:`G/(1 - \sigma GH)`.

    Parameters
    ----------
    G : scipy.signal.lti
        Forward-path model.
    H : scipy.signal.lti or float, optional
        Feedback-path model or constant gain (unity feedback by default).
    sign : {-1, 1}, optional
        Feedback sign :math:`\sigma`; ``-1`` is negative feedback.

    Returns
    -------
    scipy.signal.TransferFunction
        Closed-loop model from reference to output.

    Notes
    -----
    Pole-zero cancellations are not removed. See
    :doc:`/guides/control_systems`.

    References
    ----------
    H. S. Black, "Stabilized feedback amplifiers," *Bell Syst. Tech. J.* 13
    (1934) 1-18.

    Examples
    --------
    >>> T = feedback(transfer_function([1000], [1]), 0.01)
    >>> round(dc_gain(T), 3)
    90.909
    """
    if sign not in (-1, 1):
        raise ValueError("sign must be -1 (negative feedback) or +1")
    gn, gd = num_den(G)
    if np.ndim(H) == 0 and not hasattr(H, "num"):
        _finite(H=H)
        hn, hd = np.array([float(H)]), np.array([1.0])
    else:
        hn, hd = num_den(H)
    num = np.polymul(gn, hd)
    den = np.polysub(np.polymul(gd, hd), sign * np.polymul(gn, hn))
    return transfer_function(num, den)


def poles(sys):
    """Poles of a SISO model, the roots of its denominator.

    Parameters
    ----------
    sys : scipy.signal.lti
        SISO continuous-time model.

    Returns
    -------
    ndarray of complex
        Poles [1/s].

    Notes
    -----
    The system is asymptotically stable when every pole has a negative real
    part. See :doc:`/guides/control_systems`.

    Examples
    --------
    >>> np.sort(poles(transfer_function([1], [1, 3, 2]))).real.tolist()
    [-2.0, -1.0]
    """
    return np.roots(num_den(sys)[1])


def zeros(sys):
    """Zeros of a SISO model, the roots of its numerator.

    Parameters
    ----------
    sys : scipy.signal.lti
        SISO continuous-time model.

    Returns
    -------
    ndarray of complex
        Zeros [1/s].

    Notes
    -----
    See :doc:`/guides/control_systems`.

    Examples
    --------
    >>> zeros(transfer_function([2, 2], [1, 3, 2])).real.tolist()
    [-1.0]
    """
    return np.roots(num_den(sys)[0])


def dc_gain(sys):
    """Steady-state gain :math:`G(0)` for a constant input.

    Parameters
    ----------
    sys : scipy.signal.lti
        SISO continuous-time model.

    Returns
    -------
    float
        Output per unit input at zero frequency.

    Raises
    ------
    ValueError
        If the model has a pole at :math:`s = 0` (an integrator), so the gain
        is unbounded.

    Notes
    -----
    The steady state is reached only if the system is stable. See
    :doc:`/guides/control_systems`.

    Examples
    --------
    >>> dc_gain(transfer_function([2, 2], [1, 3, 2]))
    1.0
    """
    num, den = num_den(sys)
    if den[-1] == 0:
        raise ValueError("the model has a pole at s = 0; its DC gain is unbounded")
    return float(num[-1] / den[-1])


def step_response(sys, T=None):
    """Response to a unit step input from rest.

    Parameters
    ----------
    sys : scipy.signal.lti
        Continuous-time model.
    T : array_like, optional
        Time samples [s]; chosen automatically if omitted.

    Returns
    -------
    tuple of ndarray
        ``(t, y)``: time [s] and output per unit step.

    Notes
    -----
    Wraps :func:`scipy.signal.step`. See :doc:`/guides/control_systems`.

    Examples
    --------
    >>> t, y = step_response(transfer_function([1], [1, 1]), [0.0, 1.0])
    >>> round(float(y[-1]), 4)
    0.6321
    """
    return signal.step(sys, T=T)


def impulse_response(sys, T=None):
    """Response to a unit impulse input from rest.

    Parameters
    ----------
    sys : scipy.signal.lti
        Continuous-time model.
    T : array_like, optional
        Time samples [s]; chosen automatically if omitted.

    Returns
    -------
    tuple of ndarray
        ``(t, h)``: time [s] and impulse response [output/(input s)].

    Notes
    -----
    Wraps :func:`scipy.signal.impulse`. See :doc:`/guides/control_systems`.

    Examples
    --------
    >>> t, h = impulse_response(transfer_function([1], [1, 1]), [0.0, 1.0])
    >>> round(float(h[-1]), 4)
    0.3679
    """
    return signal.impulse(sys, T=T)
