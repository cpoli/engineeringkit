r"""Analog and digital modulation.

Amplitude modulation (Fessenden), suppressed-carrier modulation (Carson),
wideband frequency modulation (Armstrong) with Carson's bandwidth rule, and
binary/quadrature phase-shift keying, and the Gray-labelled M-PSK and
square M-QAM constellations.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import nonnegative as _nonnegative

__all__ = ["am_signal", "dsb_sc", "fm_signal", "carson_bandwidth", "bpsk", "qpsk", "psk_constellation", "qam_constellation"]


def _binary(bits, name="bits"):
    b = np.asarray(bits)
    if not b.size or not np.all((b == 0) | (b == 1)):
        raise ValueError(f"{name} must contain only zero or one")
    return b.astype(float)


def am_signal(message, t, fc, mu=1.0, carrier_amplitude=1.0):
    r"""Conventional (full-carrier) AM, :math:`s(t) = A_c[1 + \mu m(t)]\cos 2\pi f_c t`.

    Parameters
    ----------
    message : array_like
        Normalised message :math:`m(t)` with :math:`|m| \le 1` [-].
    t : array_like
        Sample times [s], same shape as ``message``.
    fc : float
        Carrier frequency [Hz], nonnegative.
    mu : float, optional
        Modulation index [-], nonnegative; ``mu > 1`` overmodulates and
        distorts an envelope detector's output.
    carrier_amplitude : float, optional
        Carrier amplitude :math:`A_c` [V].

    Returns
    -------
    ndarray
        Modulated waveform [V].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> am_signal([0.5], [0.0], 1.0, mu=0.5).tolist()
    [1.25]
    """
    m, t = np.asarray(message, dtype=float), np.asarray(t, dtype=float)
    _finite(message=m, t=t, carrier_amplitude=carrier_amplitude)
    _nonnegative(fc=fc, mu=mu)
    return carrier_amplitude * (1 + mu * m) * np.cos(2 * np.pi * fc * t)


def dsb_sc(message, t, fc):
    r"""Double-sideband suppressed-carrier modulation, :math:`s(t) = m(t)\cos 2\pi f_c t`.

    Parameters
    ----------
    message : array_like
        Message samples [V].
    t : array_like
        Sample times [s].
    fc : float
        Carrier frequency [Hz], nonnegative.

    Returns
    -------
    ndarray
        Modulated waveform [V]; demodulate by multiplying with a coherent
        carrier and low-pass filtering.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> dsb_sc([2.0], [0.0], 5.0).tolist()
    [2.0]
    """
    m, t = np.asarray(message, dtype=float), np.asarray(t, dtype=float)
    _finite(message=m, t=t)
    _nonnegative(fc=fc)
    return m * np.cos(2 * np.pi * fc * t)


def fm_signal(message, t, fc, kf):
    r"""Frequency modulation, :math:`\cos(2\pi f_c t + 2\pi k_f\int_0^t m\,d\tau)`.

    The phase integral uses the trapezoidal rule, so nonuniform time grids
    are allowed; the modulation phase is zero at the first sample.

    Parameters
    ----------
    message : array_like
        Message samples :math:`m(t)`.
    t : array_like
        Strictly increasing sample times [s], same shape as ``message``.
    fc : float
        Carrier frequency [Hz], nonnegative.
    kf : float
        Frequency sensitivity [Hz per message unit].

    Returns
    -------
    ndarray
        Unit-amplitude FM waveform.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    E. H. Armstrong, "A method of reducing disturbances in radio signaling by
    a system of frequency modulation," *Proc. IRE* 24 (1936) 689-740.

    Examples
    --------
    >>> fm_signal([1.0, 1.0], [0.0, 0.25], 0.0, 1.0).round(12).tolist()
    [1.0, 0.0]
    """
    t = np.asarray(t, dtype=float)
    m = np.asarray(message, dtype=float)
    _finite(t=t, message=m, kf=kf)
    _nonnegative(fc=fc)
    if t.ndim != 1 or m.shape != t.shape or len(t) < 2 or np.any(np.diff(t) <= 0):
        raise ValueError("message and t must be matching vectors with increasing times")
    phase = 2 * np.pi * fc * t + 2 * np.pi * kf * cumulative_trapezoid(m, t, initial=0)
    return np.cos(phase)


def carson_bandwidth(peak_deviation, message_bandwidth):
    r"""Carson's rule for the 98%-power bandwidth of FM, :math:`B \approx 2(\Delta f + W)`.

    Parameters
    ----------
    peak_deviation : float or array_like
        Peak frequency deviation :math:`\Delta f` [Hz], nonnegative.
    message_bandwidth : float or array_like
        Highest message frequency :math:`W` [Hz], nonnegative.

    Returns
    -------
    float or ndarray
        Transmission bandwidth [Hz].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    J. R. Carson, "Notes on the theory of modulation," *Proc. IRE* 10 (1922)
    57-64.

    Examples
    --------
    >>> float(carson_bandwidth(75e3, 15e3))
    180000.0
    """
    _nonnegative(peak_deviation=peak_deviation, message_bandwidth=message_bandwidth)
    return 2 * (np.asarray(peak_deviation, dtype=float) + message_bandwidth)


def bpsk(bits):
    r"""Binary phase-shift keying: map bits ``0, 1`` to antipodal symbols ``-1, +1``.

    Parameters
    ----------
    bits : array_like
        Binary data (any integer or float dtype holding only 0 and 1).

    Returns
    -------
    ndarray
        Real unit-energy symbols.

    Notes
    -----
    Inputs are converted to float before arithmetic, so unsigned zeros do
    not wrap; see :doc:`/guides/communications`.

    Examples
    --------
    >>> bpsk([0, 1, 1]).tolist()
    [-1.0, 1.0, 1.0]
    """
    return 2 * _binary(bits) - 1


def qpsk(bit_pairs):
    r"""Quadrature phase-shift keying: two bits per unit-energy complex symbol.

    The first bit of each pair sets the sign of the in-phase component and
    the second the quadrature component (Gray labelling), so each symbol is
    :math:`(\pm 1 \pm j)/\sqrt2`.

    Parameters
    ----------
    bit_pairs : array_like
        An even number of binary values.

    Returns
    -------
    ndarray
        Complex symbols, half as many as bits.

    Notes
    -----
    See :doc:`/guides/communications`.

    References
    ----------
    C. R. Cahn, "Combined digital phase and amplitude modulation communication
    systems," *IRE Trans. Commun. Syst.* 8 (1960) 150-155.

    Examples
    --------
    >>> np.round(qpsk([1, 1]) * np.sqrt(2), 12).tolist()
    [(1+1j)]
    """
    raw = np.asarray(bit_pairs)
    if raw.size % 2:
        raise ValueError("QPSK requires an even number of binary bits")
    b = _binary(raw, "bit_pairs").reshape(-1, 2)
    return (2 * b[:, 0] - 1 + 1j * (2 * b[:, 1] - 1)) / np.sqrt(2)


def _gray(n):
    return n ^ (n >> 1)


def psk_constellation(M, phase_offset=0.0):
    r"""Gray-labelled M-ary phase-shift-keying constellation on the unit circle.

    Parameters
    ----------
    M : int
        Number of phases, a power of two, at least 2.
    phase_offset : float, optional
        Phase of the first point [rad].

    Returns
    -------
    ndarray of complex
        ``M`` unit-energy points, indexed by label: entry ``k`` is the symbol
        sent for the bit pattern ``k``. Neighbouring phases carry labels that
        differ in one bit, so a symbol error to a neighbour costs one bit.

    Notes
    -----
    The point at phase :math:`2\pi i/M + \phi_0` carries the Gray label
    :math:`i \oplus \lfloor i/2 \rfloor`. See :doc:`/guides/communications`.

    References
    ----------
    C. R. Cahn, "Performance of digital phase-modulation communication
    systems," *IRE Trans. Commun. Syst.* 7 (1959) 3-6.

    Examples
    --------
    >>> np.round(np.angle(psk_constellation(4), deg=True)).tolist()
    [0.0, 90.0, -90.0, 180.0]
    """
    _integer("M", M, 2)
    if M & (M - 1):
        raise ValueError("M must be a power of two")
    _finite(phase_offset=phase_offset)
    i = np.arange(M)
    points = np.empty(M, dtype=complex)
    points[_gray(i)] = np.exp(1j * (2 * np.pi * i / M + phase_offset))
    return points


def qam_constellation(M):
    r"""Gray-labelled square M-ary quadrature-amplitude-modulation constellation with unit average energy.

    The in-phase and quadrature rails each carry :math:`\sqrt M` equally
    spaced levels :math:`\pm 1, \pm 3, \ldots`, scaled by
    :math:`\sqrt{3/2(M-1)}` so that the average symbol energy is one.

    Parameters
    ----------
    M : int
        Number of points, a power of four, at least 4 (4, 16, 64, ...).

    Returns
    -------
    ndarray of complex
        ``M`` points indexed by label; the high half of the label's bits is
        the Gray code of the in-phase level and the low half that of the
        quadrature level.

    Notes
    -----
    See :doc:`/guides/communications`.

    References
    ----------
    C. N. Campopiano and B. G. Glazer, "A coherent digital amplitude and
    phase modulation scheme," *IRE Trans. Commun. Syst.* 10 (1962) 90-95.

    Examples
    --------
    >>> points = qam_constellation(16)
    >>> round(float(np.mean(np.abs(points) ** 2)), 12)
    1.0
    """
    _integer("M", M, 4)
    k = int(M).bit_length() - 1
    if M & (M - 1) or k % 2:
        raise ValueError("M must be a power of four")
    side = 1 << (k // 2)
    i = np.arange(side)
    levels = np.empty(side)
    levels[_gray(i)] = 2 * i - (side - 1)
    labels = np.arange(M)
    points = levels[labels >> (k // 2)] + 1j * levels[labels & (side - 1)]
    return points * np.sqrt(3 / (2 * (M - 1)))
