r"""Nonlinear loops: describing functions, limit-cycle prediction and the phase plane.

Kochenburger (1950) analysed relay ("contactor") servomechanisms by
replacing the nonlinearity with its *describing function*
:math:`N(A)`, the complex gain from a sinusoid :math:`A\sin\omega t` to the
fundamental of its output. A loop :math:`L(s)` closed through it can
sustain an oscillation of amplitude :math:`A` and frequency :math:`\omega`
where :math:`L(j\omega) = -1/N(A)`. Andronov, Vitt and Khaikin (1937)
studied such self-sustained oscillations in the phase plane
:math:`(x, \dot x)`, where they are closed orbits.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.control_systems.core.base import PhasePortraitResult
from engineeringkit.control_systems.systems.frequency import stability_margins
from engineeringkit.integrators import integrate_ode

__all__ = [
    "describing_function",
    "relay_describing_function",
    "saturation_describing_function",
    "dead_zone_describing_function",
    "limit_cycle",
    "phase_portrait",
]


def describing_function(nonlinearity, amplitude, samples=4096):
    r"""Describing function of a memoryless nonlinearity, from the first Fourier harmonic of its output.

    Parameters
    ----------
    nonlinearity : callable
        ``nonlinearity(x)`` evaluated elementwise on an array of inputs.
    amplitude : float or array_like
        Input amplitude :math:`A` [input units], positive.
    samples : int, optional
        Points per period of the quadrature, at least 16.

    Returns
    -------
    complex or ndarray of complex
        :math:`N(A) = (b_1 + j a_1)/A`, where :math:`b_1` and :math:`a_1`
        are the sine and cosine coefficients of the output's fundamental
        [output units per input unit].

    Notes
    -----
    The rectangle rule on a uniform grid is spectrally accurate for smooth
    nonlinearities and first-order accurate across a jump. An odd,
    single-valued nonlinearity has a real :math:`N`. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/control_systems`.

    Examples
    --------
    >>> round(float(describing_function(lambda x: x**3, 2.0).real), 6)
    3.0
    """
    _positive(amplitude=amplitude)
    _integer("samples", samples, minimum=16)
    A = np.asarray(amplitude, float)
    theta = 2 * np.pi * np.arange(samples) / samples
    y = np.asarray(nonlinearity(A[..., None] * np.sin(theta)), float)
    b1 = 2 * np.mean(y * np.sin(theta), axis=-1)
    a1 = 2 * np.mean(y * np.cos(theta), axis=-1)
    N = (b1 + 1j * a1) / A
    return complex(N) if N.ndim == 0 else N


def relay_describing_function(amplitude, output_level, hysteresis=0.0):
    r"""Describing function of a relay, ideal or with hysteresis.

    Parameters
    ----------
    amplitude : float or array_like
        Input amplitude :math:`A` [input units], greater than ``hysteresis``.
    output_level : float
        Relay output :math:`\pm M` [output units], positive.
    hysteresis : float, optional
        Switching threshold :math:`h` [input units], nonnegative: the relay
        switches up when the input rises past :math:`h` and down when it
        falls past :math:`-h`.

    Returns
    -------
    float, complex or ndarray
        :math:`N(A) = \frac{4M}{\pi A}e^{-j\arcsin(h/A)}`; real when
        :math:`h = 0` [output units per input unit].

    Notes
    -----
    Hysteresis makes the relay lag the input by :math:`\arcsin(h/A)`.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/control_systems`.

    References
    ----------
    R. J. Kochenburger, "A frequency response method for analyzing and
    synthesizing contactor servomechanisms," *Trans. AIEE* 69 (1950)
    270-284.

    Examples
    --------
    >>> round(float(relay_describing_function(1.0, np.pi / 4)), 12)
    1.0
    """
    _positive(amplitude=amplitude, output_level=output_level)
    _nonnegative(hysteresis=hysteresis)
    A = np.asarray(amplitude, float)
    if np.any(A <= hysteresis):
        raise ValueError("amplitude must exceed the hysteresis")
    N = 4 * output_level / (np.pi * A)
    if hysteresis:
        r = hysteresis / A
        N = N * (np.sqrt(1 - r * r) - 1j * r)
    return N[()]


def _saturation_gain(A, limit, slope):
    r = np.minimum(limit / A, 1.0)
    return 2 * slope / np.pi * (np.arcsin(r) + r * np.sqrt(1 - r * r))


def saturation_describing_function(amplitude, limit, slope=1.0):
    r"""Describing function of a saturation: linear with gain ``slope`` for inputs up to ``limit``.

    Parameters
    ----------
    amplitude : float or array_like
        Input amplitude :math:`A` [input units], positive.
    limit : float
        Input :math:`a` at which the output saturates [input units],
        positive; the output saturates at :math:`\pm k a`.
    slope : float, optional
        Linear gain :math:`k` [output units per input unit], positive.

    Returns
    -------
    float or ndarray
        :math:`N = k` for :math:`A \le a`, else
        :math:`\frac{2k}{\pi}\left[\arcsin\frac aA + \frac aA
        \sqrt{1 - (a/A)^2}\right]`.

    Notes
    -----
    :math:`N` falls from :math:`k` towards :math:`4ka/(\pi A)`, the relay
    value, as the amplitude grows. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/control_systems`.

    Examples
    --------
    >>> float(saturation_describing_function(0.5, 1.0, 2.0))
    2.0
    """
    _positive(amplitude=amplitude, limit=limit, slope=slope)
    return _saturation_gain(np.asarray(amplitude, float), limit, slope)[()]


def dead_zone_describing_function(amplitude, width, slope=1.0):
    r"""Describing function of a dead zone: zero output for inputs below ``width``, then gain ``slope``.

    Parameters
    ----------
    amplitude : float or array_like
        Input amplitude :math:`A` [input units], positive.
    width : float
        Half-width :math:`\delta` of the dead zone [input units], positive.
    slope : float, optional
        Gain :math:`k` outside the dead zone [output units per input unit],
        positive.

    Returns
    -------
    float or ndarray
        :math:`N = 0` for :math:`A \le \delta`, else :math:`k` minus the
        saturation describing function with limit :math:`\delta`.

    Notes
    -----
    A dead zone is a linear gain minus a saturation, and the describing
    function is linear in the nonlinearity. Inputs are raw SI numbers with
    no dimensional metadata; see :doc:`/guides/control_systems`.

    Examples
    --------
    >>> float(dead_zone_describing_function(0.5, 1.0))
    0.0
    """
    _positive(amplitude=amplitude, width=width, slope=slope)
    return (slope - _saturation_gain(np.asarray(amplitude, float), width, slope))[()]


def limit_cycle(loop, describing, amplitude_bounds, w=None, samples=400):
    r"""Limit cycles predicted by the describing function for a loop closed through a real nonlinearity.

    With a real :math:`N(A)`, :math:`L(j\omega) = -1/N(A)` needs
    :math:`\angle L = -180^\circ` (a phase crossover :math:`\omega_{180}`)
    and :math:`N(A) = 1/|L(j\omega_{180})|`, the gain margin there.

    Parameters
    ----------
    loop : scipy.signal.lti
        Linear part :math:`L(s)` of the loop, open-loop stable.
    describing : callable
        ``describing(A)`` returning the real describing function, vectorised
        over amplitudes.
    amplitude_bounds : tuple of float
        ``(A_min, A_max)`` [input units], the positive amplitude range
        searched.
    w : array_like, optional
        Positive angular frequencies [rad/s] for finding the phase
        crossovers, as in :func:`stability_margins`.
    samples : int, optional
        Log-spaced amplitudes scanned for sign changes, at least 2.

    Returns
    -------
    list of dict
        One dict per predicted limit cycle, in increasing frequency then
        amplitude: ``amplitude`` [input units], ``frequency`` [rad/s],
        ``period`` [s] and ``stable`` (bool). Empty if none is predicted.

    Raises
    ------
    ValueError
        If the bounds are not increasing and positive, or the describing
        function is not real.

    Notes
    -----
    By Loeb's criterion, a cycle is stable when :math:`N` falls through the
    gain margin as :math:`A` grows: a larger oscillation then sees a loop
    gain below the critical one and decays back. The method assumes
    :math:`L` filters out the harmonics, so it gives estimates, not proofs.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/control_systems`.

    References
    ----------
    R. J. Kochenburger, "A frequency response method for analyzing and
    synthesizing contactor servomechanisms," *Trans. AIEE* 69 (1950)
    270-284.

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> relay = lambda A: relay_describing_function(A, 1.0)
    >>> (cycle,) = limit_cycle(transfer_function([1], [1, 3, 3, 1]), relay, (0.01, 10))
    >>> round(cycle["amplitude"], 4), round(cycle["frequency"], 3)
    (0.1592, 1.732)
    """
    a_min, a_max = amplitude_bounds
    _positive(a_min=a_min, a_max=a_max)
    _integer("samples", samples, minimum=2)
    if a_min >= a_max:
        raise ValueError("amplitude_bounds must be increasing")
    amplitudes = np.logspace(np.log10(a_min), np.log10(a_max), samples)
    values = np.asarray(describing(amplitudes))
    if np.iscomplexobj(values) and np.any(values.imag != 0):
        raise ValueError("only real describing functions (odd, single-valued nonlinearities) are supported")
    values = values.real.astype(float)
    margins = stability_margins(loop, w)
    cycles = []
    for frequency, gain_margin in zip(margins["phase_cross_frequencies"], margins["gain_margins"]):
        excess = values - gain_margin
        for i in np.flatnonzero(np.sign(excess[:-1]) * np.sign(excess[1:]) < 0):
            amplitude = optimize.brentq(lambda a: float(np.real(describing(a))) - gain_margin, amplitudes[i], amplitudes[i + 1])
            cycles.append(
                {
                    "amplitude": float(amplitude),
                    "frequency": float(frequency),
                    "period": float(2 * np.pi / frequency),
                    "stable": bool(excess[i] > 0),
                }
            )
    return cycles


def phase_portrait(acceleration, x0, v0, t):
    r"""Trajectory of :math:`\ddot x = f(x, \dot x)` in the phase plane.

    Parameters
    ----------
    acceleration : callable
        ``acceleration(x, v)`` returning :math:`\ddot x` [position units/s²];
        it is also evaluated on arrays to draw the direction field.
    x0, v0 : float
        Initial position [position units] and velocity [position units/s].
    t : array_like
        At least two strictly increasing sample times [s], fine enough for a
        fixed-step RK4 integration.

    Returns
    -------
    PhasePortraitResult
        Time [s], position and velocity, and the right-hand side for plotting
        the direction field.

    Notes
    -----
    Fixed RK4 steps tolerate a discontinuous :math:`f`, such as a relay, at
    first-order accuracy across each switch. Time is eliminated in the
    plane: the slope of a trajectory is :math:`dv/dx = f(x, v)/v`. Inputs
    are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/control_systems`.

    References
    ----------
    A. A. Andronov and S. E. Khaikin, *Theory of Oscillations* (ONTI,
    Moscow, 1937; with A. A. Vitt, 2nd ed., 1959).

    Examples
    --------
    >>> r = phase_portrait(lambda x, v: -x, 1.0, 0.0, np.linspace(0, np.pi, 2001))
    >>> round(r.summary()["final_x"], 6)
    -1.0
    """
    _finite(x0=x0, v0=v0)

    def rhs(_, state):
        return np.array([state[1], acceleration(state[0], state[1])], float)

    sol = integrate_ode(rhs, [x0, v0], t, method="rk4")
    return PhasePortraitResult(sol.t, sol.history[:, 0], sol.history[:, 1], acceleration)
