r"""Digital filters: Butterworth, Chebyshev and elliptic IIR, the bilinear transform, window and equiripple FIR, frequency response.

Hurewicz (1947) described a sampled system by its z-transform
:math:`H(z) = B(z)/A(z)`; its frequency response is :math:`H` on the unit
circle, :math:`z = e^{j2\pi f/f_s}`. Butterworth's (1930) maximally flat
low-pass, carried into discrete time by the bilinear transform, has
:math:`|H|^2 = 1/(1 + (\tan(\pi f/f_s)/\tan(\pi f_c/f_s))^{2n})`.
Chebyshev's (1854) polynomials spread the error evenly over the passband
or the stopband instead, and Cauer's (1931) elliptic filters over both,
for the lowest order that meets a specification. Golden and Kaiser (1964)
designed digital filters from analog ones with the bilinear transform
:math:`s = K(1 - z^{-1})/(1 + z^{-1})`, prewarping the band edges.
Kaiser (1966) designed FIR filters by windowing the ideal impulse response,
with formulas for the length and window shape that meet a given
attenuation; Parks and McClellan (1972) found the linear-phase FIR filter
whose weighted error ripples evenly, the best in the minimax sense.
Jackson (1970) showed why high-order recursive filters are built as
cascades of second-order sections: the roots of one high-degree polynomial
are too sensitive to rounding in its coefficients.
"""

from __future__ import annotations

import numpy as np
from scipy import signal, special

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive

__all__ = [
    "butterworth_lowpass",
    "butterworth_highpass",
    "chebyshev1_lowpass",
    "chebyshev2_lowpass",
    "elliptic_lowpass",
    "iir_order",
    "prewarp_frequency",
    "bilinear_transform",
    "fir_lowpass",
    "kaiser_order",
    "equiripple_lowpass",
    "equiripple_order",
    "filter_signal",
    "frequency_response",
]


def _iir_design(design, order, cutoff, fs, btype, output, *specification):
    _integer("order", order, 1)
    _positive(cutoff=cutoff, fs=fs)
    if cutoff >= fs / 2:
        raise ValueError("cutoff must lie below the Nyquist frequency fs/2")
    if output not in ("ba", "sos"):
        raise ValueError('output must be "ba" or "sos"')
    return design(order, *specification, cutoff, btype=btype, fs=fs, output=output)


def _butterworth(order, cutoff, fs, btype, output):
    return _iir_design(signal.butter, order, cutoff, fs, btype, output)


def butterworth_lowpass(order, cutoff, fs, output="ba"):
    r"""Digital Butterworth low-pass filter by the bilinear transform.

    Parameters
    ----------
    order : int
        Filter order :math:`n`, at least 1.
    cutoff : float
        -3 dB frequency :math:`f_c` [Hz], between 0 and :math:`f_s/2`.
    fs : float
        Sampling rate [Hz], positive.
    output : {"ba", "sos"}, optional
        ``"ba"`` for numerator and denominator polynomials in
        :math:`z^{-1}`; ``"sos"`` for a cascade of second-order sections.

    Returns
    -------
    tuple of ndarray or ndarray
        ``(b, a)`` [-], or an ``(n_sections, 6)`` array of
        ``[b0, b1, b2, 1, a1, a2]`` rows [-].

    Notes
    -----
    :math:`|H(f)|^2 = 1/(1 + (\tan(\pi f/f_s)/\tan(\pi f_c/f_s))^{2n})`: the
    first :math:`2n-1` derivatives vanish at :math:`f = 0`, and the gain is
    :math:`1/\sqrt2` at :math:`f_c` for every order. Above about order 8,
    or for a cutoff far below :math:`f_s/2`, the ``"ba"`` polynomials lose
    precision and the filter can come out wrong or unstable; use ``"sos"``.
    See :doc:`/guides/signal_processing`.

    References
    ----------
    S. Butterworth, "On the theory of filter amplifiers," *Wireless
    Engineer* 7 (1930) 536-541; L. B. Jackson, "On the interaction of
    roundoff noise and dynamic range in digital filters," *Bell Syst. Tech.
    J.* 49 (1970) 159-184.

    Examples
    --------
    >>> b, a = butterworth_lowpass(1, 250.0, 1000.0)
    >>> np.round(b, 12).tolist()
    [0.5, 0.5]
    """
    return _butterworth(order, cutoff, fs, "lowpass", output)


def butterworth_highpass(order, cutoff, fs, output="ba"):
    r"""Digital Butterworth high-pass filter by the bilinear transform.

    Parameters
    ----------
    order : int
        Filter order :math:`n`, at least 1.
    cutoff : float
        -3 dB frequency :math:`f_c` [Hz], between 0 and :math:`f_s/2`.
    fs : float
        Sampling rate [Hz], positive.
    output : {"ba", "sos"}, optional
        Polynomials or second-order sections, as in
        :func:`butterworth_lowpass`.

    Returns
    -------
    tuple of ndarray or ndarray
        ``(b, a)`` or second-order sections [-].

    Notes
    -----
    :math:`|H(f)|^2 = 1/(1 + (\tan(\pi f_c/f_s)/\tan(\pi f/f_s))^{2n})`. See
    :func:`butterworth_lowpass` on precision, and
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> b, a = butterworth_highpass(1, 250.0, 1000.0)
    >>> b.tolist()
    [0.5, -0.5]
    """
    return _butterworth(order, cutoff, fs, "highpass", output)


def chebyshev1_lowpass(order, ripple_db, cutoff, fs, output="ba"):
    r"""Digital Chebyshev type I low-pass filter (equiripple passband) by the bilinear transform.

    Parameters
    ----------
    order : int
        Filter order :math:`n`, at least 1.
    ripple_db : float
        Peak-to-peak passband ripple :math:`R` [dB], positive.
    cutoff : float
        Passband edge [Hz], where the gain last equals :math:`-R` dB,
        between 0 and :math:`f_s/2`.
    fs : float
        Sampling rate [Hz], positive.
    output : {"ba", "sos"}, optional
        Polynomials or second-order sections, as in
        :func:`butterworth_lowpass`.

    Returns
    -------
    tuple of ndarray or ndarray
        ``(b, a)`` or second-order sections [-].

    Notes
    -----
    The analog prototype has :math:`|H|^2 = 1/(1 + \varepsilon^2 T_n^2(\Omega/\Omega_p))`,
    with :math:`T_n` the Chebyshev polynomial and
    :math:`\varepsilon^2 = 10^{R/10} - 1`. :math:`|T_n| \le 1` on
    :math:`[-1, 1]`, so the passband gain swings :math:`n` times between 1
    and :math:`1/\sqrt{1+\varepsilon^2}`; outside it :math:`T_n` grows as
    fast as any polynomial of degree :math:`n` can, which is why the
    transition is steeper than Butterworth's for the same order. Frequencies
    are mapped by the bilinear transform with the edge prewarped. See
    :doc:`/guides/signal_processing`.

    References
    ----------
    P. L. Chebyshev, "Théorie des mécanismes connus sous le nom de
    parallélogrammes," *Mém. Acad. Sci. St.-Pétersbourg, Savants étrangers*
    7 (1854) 539-568; A. V. Oppenheim and R. W. Schafer, *Discrete-Time
    Signal Processing*, 3rd ed., Pearson (2010), sec. 7.3.

    Examples
    --------
    >>> b, a = chebyshev1_lowpass(4, 1.0, 100.0, 1000.0)
    >>> f, H = frequency_response(b, a, worN=[0.0, 100.0], fs=1000.0)
    >>> np.round(20 * np.log10(np.abs(H)), 6).tolist()
    [-1.0, -1.0]
    """
    _positive(ripple_db=ripple_db)
    return _iir_design(signal.cheby1, order, cutoff, fs, "lowpass", output, ripple_db)


def chebyshev2_lowpass(order, attenuation_db, cutoff, fs, output="ba"):
    r"""Digital Chebyshev type II (inverse Chebyshev) low-pass filter: flat passband, equiripple stopband.

    Parameters
    ----------
    order : int
        Filter order :math:`n`, at least 1.
    attenuation_db : float
        Minimum stopband attenuation :math:`A` [dB], positive.
    cutoff : float
        Stopband edge [Hz], where the gain first reaches :math:`-A` dB,
        between 0 and :math:`f_s/2`.
    fs : float
        Sampling rate [Hz], positive.
    output : {"ba", "sos"}, optional
        Polynomials or second-order sections, as in
        :func:`butterworth_lowpass`.

    Returns
    -------
    tuple of ndarray or ndarray
        ``(b, a)`` or second-order sections [-].

    Notes
    -----
    The prototype :math:`|H|^2 = 1/(1 + 1/(\varepsilon^2 T_n^2(\Omega_s/\Omega)))`
    is maximally flat at DC and has zeros on the unit circle in the
    stopband, where the gain ripples between 0 and :math:`10^{-A/20}`. See
    :func:`chebyshev1_lowpass` and :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> b, a = chebyshev2_lowpass(4, 40.0, 150.0, 1000.0)
    >>> f, H = frequency_response(b, a, worN=[0.0, 150.0], fs=1000.0)
    >>> np.round(20 * np.log10(np.abs(H)), 6).tolist()
    [0.0, -40.0]
    """
    _positive(attenuation_db=attenuation_db)
    return _iir_design(signal.cheby2, order, cutoff, fs, "lowpass", output, attenuation_db)


def elliptic_lowpass(order, ripple_db, attenuation_db, cutoff, fs, output="ba"):
    r"""Digital elliptic (Cauer) low-pass filter: equiripple in both the passband and the stopband.

    Parameters
    ----------
    order : int
        Filter order :math:`n`, at least 1.
    ripple_db : float
        Peak-to-peak passband ripple :math:`R` [dB], positive.
    attenuation_db : float
        Minimum stopband attenuation :math:`A` [dB], above ``ripple_db``.
    cutoff : float
        Passband edge [Hz], between 0 and :math:`f_s/2`.
    fs : float
        Sampling rate [Hz], positive.
    output : {"ba", "sos"}, optional
        Polynomials or second-order sections, as in
        :func:`butterworth_lowpass`.

    Returns
    -------
    tuple of ndarray or ndarray
        ``(b, a)`` or second-order sections [-].

    Notes
    -----
    The prototype :math:`|H|^2 = 1/(1 + \varepsilon^2 R_n^2(\Omega/\Omega_p))`
    uses Cauer's elliptic rational function :math:`R_n`, which ripples
    between :math:`\pm 1` in the passband and stays above
    :math:`\sqrt{(10^{A/10}-1)/\varepsilon^2}` in the stopband. Spreading
    the error over both bands gives the narrowest transition of any filter
    of order :math:`n` with that ripple and attenuation, at the cost of the
    most nonlinear phase. See :func:`iir_order` and
    :doc:`/guides/signal_processing`.

    References
    ----------
    W. Cauer, *Siebschaltungen*, VDI-Verlag, Berlin (1931); W. Cauer,
    *Synthesis of Linear Communication Networks*, McGraw-Hill (1958).

    Examples
    --------
    >>> b, a = elliptic_lowpass(4, 0.5, 60.0, 100.0, 1000.0)
    >>> f, H = frequency_response(b, a, worN=[100.0], fs=1000.0)
    >>> round(float(20 * np.log10(np.abs(H[0]))), 6)
    -0.5
    """
    _positive(ripple_db=ripple_db, attenuation_db=attenuation_db)
    if attenuation_db <= ripple_db:
        raise ValueError("attenuation_db must exceed ripple_db")
    return _iir_design(signal.ellip, order, cutoff, fs, "lowpass", output, ripple_db, attenuation_db)


def iir_order(pass_edge, stop_edge, ripple_db, attenuation_db, fs, family="butterworth"):
    r"""Minimum order of a digital low-pass IIR filter that meets a specification.

    Parameters
    ----------
    pass_edge : float
        Passband edge :math:`f_p` [Hz].
    stop_edge : float
        Stopband edge :math:`f_s'` [Hz], with :math:`0 < f_p < f_s' < f_s/2`.
    ripple_db : float
        Largest passband loss :math:`R` [dB], positive.
    attenuation_db : float
        Smallest stopband loss :math:`A` [dB], above ``ripple_db``.
    fs : float
        Sampling rate [Hz], positive.
    family : {"butterworth", "chebyshev", "elliptic"}, optional
        Approximation; ``"chebyshev"`` holds for both types.

    Returns
    -------
    int
        The smallest order :math:`n`.

    Notes
    -----
    With the prewarped selectivity :math:`k = \tan(\pi f_p/f_s)/\tan(\pi f_s'/f_s)`
    and discrimination :math:`k_1 = \sqrt{(10^{R/10}-1)/(10^{A/10}-1)}`,
    :math:`n` is the ceiling of :math:`\log k_1/\log k` (Butterworth),
    :math:`\cosh^{-1}(1/k_1)/\cosh^{-1}(1/k)` (Chebyshev), or
    :math:`K(k)K'(k_1)/(K'(k)K(k_1))` (elliptic), with :math:`K` the
    complete elliptic integral of the first kind and :math:`K'(k) = K(\sqrt{1-k^2})`.
    The orders agree with SciPy's ``buttord``, ``cheb1ord`` and ``ellipord``.
    See :doc:`/guides/signal_processing`.

    References
    ----------
    T. W. Parks and C. S. Burrus, *Digital Filter Design*, Wiley (1987),
    ch. 7; A. V. Oppenheim and R. W. Schafer, *Discrete-Time Signal
    Processing*, 3rd ed., Pearson (2010), sec. 7.3.

    Examples
    --------
    >>> [iir_order(100.0, 150.0, 1.0, 60.0, 1000.0, family) for family in ("butterworth", "chebyshev", "elliptic")]
    [17, 9, 6]
    """
    _positive(pass_edge=pass_edge, stop_edge=stop_edge, ripple_db=ripple_db, attenuation_db=attenuation_db, fs=fs)
    if not pass_edge < stop_edge < fs / 2:
        raise ValueError("need 0 < pass_edge < stop_edge < fs/2")
    if attenuation_db <= ripple_db:
        raise ValueError("attenuation_db must exceed ripple_db")
    k = np.tan(np.pi * pass_edge / fs) / np.tan(np.pi * stop_edge / fs)
    k1 = np.sqrt(np.expm1(np.log(10) * ripple_db / 10) / np.expm1(np.log(10) * attenuation_db / 10))
    if family == "butterworth":
        n = np.log(k1) / np.log(k)
    elif family == "chebyshev":
        n = np.arccosh(1 / k1) / np.arccosh(1 / k)
    elif family == "elliptic":
        n = special.ellipk(k**2) * special.ellipkm1(k1**2) / (special.ellipkm1(k**2) * special.ellipk(k1**2))
    else:
        raise ValueError('family must be "butterworth", "chebyshev" or "elliptic"')
    return int(np.ceil(n))


def prewarp_frequency(frequency, fs):
    r"""Analog frequency that the bilinear transform maps onto a given digital frequency.

    Parameters
    ----------
    frequency : float or array_like
        Digital frequency :math:`f` [Hz], from 0 to below :math:`f_s/2`.
    fs : float
        Sampling rate [Hz], positive.

    Returns
    -------
    float or ndarray
        Analog frequency :math:`f_a = (f_s/\pi)\tan(\pi f/f_s)` [Hz].

    Notes
    -----
    The bilinear transform with :math:`K = 2f_s` takes the analog axis
    :math:`s = j2\pi f_a` onto the unit circle at
    :math:`f = (f_s/\pi)\arctan(\pi f_a/f_s)`, compressing all of
    :math:`[0, \infty)` into :math:`[0, f_s/2)`. Designing the analog
    prototype at :math:`f_a` puts its edge at :math:`f` after the
    transform. :math:`f_a \approx f` only for :math:`f \ll f_s`. See
    :func:`bilinear_transform` and :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> round(float(prewarp_frequency(250.0, 1000.0)), 6)
    318.309886
    """
    _positive(fs=fs)
    _finite(frequency=frequency)
    f = np.asarray(frequency, dtype=float)
    if np.any(f < 0) or np.any(f >= fs / 2):
        raise ValueError("frequency must lie in [0, fs/2)")
    return (fs / np.pi * np.tan(np.pi * f / fs))[()]


def bilinear_transform(b, a, fs, prewarp=None):
    r"""Digital filter :math:`H(z)` from an analog :math:`H(s)` by the bilinear transform, written out.

    Parameters
    ----------
    b, a : array_like
        Analog numerator and denominator coefficients in descending powers
        of :math:`s` [-], with ``a`` of degree at least that of ``b``.
    fs : float
        Sampling rate [Hz], positive.
    prewarp : float, optional
        Frequency [Hz], below :math:`f_s/2`, at which the digital response
        is made to equal the analog one exactly. Without it :math:`K = 2f_s`.

    Returns
    -------
    tuple of ndarray
        ``(b, a)``: coefficients in :math:`z^{-1}` [-], with ``a[0] = 1``.

    Notes
    -----
    Substituting :math:`s = K(1 - z^{-1})/(1 + z^{-1})` and multiplying
    through by :math:`(1 + z^{-1})^N`, each term :math:`c_m s^m` of degree
    :math:`N` polynomials becomes :math:`c_m K^m (1 - z^{-1})^m (1 + z^{-1})^{N-m}`.
    The left half :math:`s`-plane maps inside the unit circle, so stable
    analog filters give stable digital ones, and the :math:`j\Omega` axis maps
    onto the circle with :math:`\Omega = K\tan(\omega/2)`. Prewarping at
    :math:`f_0` sets :math:`K = 2\pi f_0/\tan(\pi f_0/f_s)`. The result
    equals ``scipy.signal.bilinear`` with ``fs = K/2``. See
    :func:`prewarp_frequency` and :doc:`/guides/signal_processing`.

    References
    ----------
    R. M. Golden and J. F. Kaiser, "Design of wideband sampled-data
    filters," *Bell Syst. Tech. J.* 43 (1964) 1533-1546.

    Examples
    --------
    >>> b, a = bilinear_transform([1.0], [1.0, 1.0], 0.5)
    >>> b.tolist(), a.tolist()
    ([0.5, 0.5], [1.0, 0.0])
    """
    _positive(fs=fs)
    b = np.atleast_1d(np.asarray(b, dtype=float))
    a = np.atleast_1d(np.asarray(a, dtype=float))
    _finite(b=b, a=a)
    if b.ndim != 1 or a.ndim != 1 or a[0] == 0 or b.size > a.size:
        raise ValueError("b and a must be 1D, a[0] nonzero and deg(b) <= deg(a)")
    if prewarp is None:
        K = 2.0 * fs
    else:
        _positive(prewarp=prewarp)
        if prewarp >= fs / 2:
            raise ValueError("prewarp must lie below fs/2")
        K = 2 * np.pi * prewarp / np.tan(np.pi * prewarp / fs)
    n = a.size - 1

    def substitute(c):
        out = np.zeros(n + 1)
        for m, coefficient in enumerate(c[::-1]):
            out += coefficient * K**m * np.polymul(_binomial_power([1.0, -1.0], m), _binomial_power([1.0, 1.0], n - m))
        return out

    num, den = substitute(b), substitute(a)
    return num / den[0], den / den[0]


def _binomial_power(p, m):
    """Coefficients of ``p`` raised to the power ``m``."""
    out = np.array([1.0])
    for _ in range(m):
        out = np.polymul(out, p)
    return out


def fir_lowpass(numtaps, cutoff, fs, window="hamming"):
    r"""Linear-phase FIR low-pass filter by the window method (SciPy's ``firwin``).

    Parameters
    ----------
    numtaps : int
        Number of coefficients :math:`N`, at least 1.
    cutoff : float
        Cutoff frequency [Hz], between 0 and :math:`f_s/2`; the gain there
        is about one half (-6 dB).
    fs : float
        Sampling rate [Hz], positive.
    window : str or tuple, optional
        Window, e.g. ``"hamming"`` or ``("kaiser", beta)``; see
        :func:`kaiser_order`.

    Returns
    -------
    ndarray
        Symmetric impulse response [-], scaled to unit gain at DC.

    Notes
    -----
    The ideal low-pass impulse response
    :math:`2(f_c/f_s)\operatorname{sinc}(2 f_c (n - (N-1)/2)/f_s)` is
    truncated to :math:`N` samples and tapered by the window. Symmetry gives
    a constant group delay of :math:`(N-1)/2` samples. See
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> round(float(fir_lowpass(31, 100.0, 1000.0).sum()), 12)
    1.0
    """
    _integer("numtaps", numtaps, 1)
    _positive(cutoff=cutoff, fs=fs)
    return signal.firwin(numtaps, cutoff, fs=fs, window=window)


def kaiser_order(attenuation_db, transition_width, fs):
    r"""Kaiser's formulas for the FIR length and window shape that meet a stopband attenuation.

    Parameters
    ----------
    attenuation_db : float
        Required stopband attenuation (and passband ripple)
        :math:`A = -20\log_{10}\delta` [dB], positive.
    transition_width : float
        Width of the transition band :math:`\Delta f` [Hz], positive and
        below :math:`f_s/2`.
    fs : float
        Sampling rate [Hz], positive.

    Returns
    -------
    tuple
        ``(numtaps, beta)``: number of coefficients (int) and Kaiser window
        parameter :math:`\beta` [-], for ``window=("kaiser", beta)``.

    Notes
    -----
    All parameters are scalar-only: the result sizes one filter.

    With :math:`\Delta\omega = 2\pi\Delta f/f_s`,
    :math:`N - 1 = \lceil (A - 7.95)/(2.285\,\Delta\omega) \rceil` and
    :math:`\beta = 0.1102(A - 8.7)` for :math:`A > 50`,
    :math:`0.5842(A-21)^{0.4} + 0.07886(A-21)` for :math:`21 \le A \le 50`,
    and 0 below. The formulas are empirical and usually land within about
    a decibel of :math:`A`. See :doc:`/guides/signal_processing`.

    References
    ----------
    J. F. Kaiser, "Digital filters," ch. 7 in F. F. Kuo and J. F. Kaiser
    (eds.), *System Analysis by Digital Computer*, Wiley (1966); J. F.
    Kaiser, "Nonrecursive digital filter design using the I0-sinh window
    function," *Proc. IEEE Int. Symp. Circuits Syst.* (1974) 20-23.

    Examples
    --------
    >>> numtaps, beta = kaiser_order(60.0, 50.0, 1000.0)
    >>> numtaps, round(beta, 4)
    (74, 5.6533)
    """
    _positive(attenuation_db=attenuation_db, transition_width=transition_width, fs=fs)
    if transition_width >= fs / 2:
        raise ValueError("transition_width must lie below fs/2")
    a = float(attenuation_db)
    if a > 50:
        beta = 0.1102 * (a - 8.7)
    elif a >= 21:
        beta = 0.5842 * (a - 21) ** 0.4 + 0.07886 * (a - 21)
    else:
        beta = 0.0
    numtaps = int(np.ceil((a - 7.95) / (2.285 * 2 * np.pi * transition_width / fs))) + 1
    return max(numtaps, 1), beta


def equiripple_lowpass(numtaps, pass_edge, stop_edge, fs, weight=1.0):
    r"""Parks-McClellan equiripple linear-phase FIR low-pass filter (SciPy's ``remez``).

    Parameters
    ----------
    numtaps : int
        Number of coefficients :math:`N`, at least 3.
    pass_edge : float
        Passband edge [Hz].
    stop_edge : float
        Stopband edge [Hz], with :math:`0 < f_p < f_s' < f_s/2`.
    fs : float
        Sampling rate [Hz], positive.
    weight : float, optional
        Stopband error weight relative to the passband [-], positive; the
        designed deviations satisfy :math:`\delta_p = \text{weight}\cdot\delta_s`.

    Returns
    -------
    ndarray
        Symmetric impulse response [-].

    Notes
    -----
    Among all symmetric FIR filters of length :math:`N`, the result
    minimizes the largest weighted error
    :math:`\max |W(f)(H(f) - D(f))|` over both bands. By the alternation
    theorem the optimum's error touches its peak, with alternating sign, at
    least :math:`(N+1)/2 + 1` frequencies, and Remez's exchange algorithm
    finds them. For the same length it beats any window design, whose
    largest ripple sits next to the band edge. The transition band is a
    "don't care" region. See :func:`equiripple_order` and
    :doc:`/guides/signal_processing`.

    References
    ----------
    T. W. Parks and J. H. McClellan, "Chebyshev approximation for
    nonrecursive digital filters with linear phase," *IEEE Trans. Circuit
    Theory* 19 (1972) 189-194, doi:10.1109/TCT.1972.1083419.

    Examples
    --------
    >>> h = equiripple_lowpass(41, 100.0, 150.0, 1000.0)
    >>> h.size, bool(np.allclose(h, h[::-1]))
    (41, True)
    """
    _integer("numtaps", numtaps, 3)
    _positive(pass_edge=pass_edge, stop_edge=stop_edge, fs=fs, weight=weight)
    if not pass_edge < stop_edge < fs / 2:
        raise ValueError("need 0 < pass_edge < stop_edge < fs/2")
    return signal.remez(numtaps, [0.0, pass_edge, stop_edge, fs / 2], [1.0, 0.0], weight=[1.0, weight], fs=fs)


def equiripple_order(passband_deviation, stopband_deviation, transition_width, fs):
    r"""Kaiser's estimate of the Parks-McClellan filter length for given deviations.

    Parameters
    ----------
    passband_deviation : float
        Passband deviation :math:`\delta_p` [-], in (0, 1).
    stopband_deviation : float
        Stopband deviation :math:`\delta_s` [-], in (0, 1).
    transition_width : float
        Transition band :math:`\Delta f` [Hz], positive and below :math:`f_s/2`.
    fs : float
        Sampling rate [Hz], positive.

    Returns
    -------
    int
        Number of coefficients :math:`N`.

    Notes
    -----
    :math:`N - 1 \approx (-20\log_{10}\sqrt{\delta_p\delta_s} - 13)/(14.6\,\Delta f/f_s)`,
    rounded up. The estimate is empirical and can be a few taps short;
    check the design and add taps if needed. Pass
    ``weight = passband_deviation / stopband_deviation`` to
    :func:`equiripple_lowpass`. See :doc:`/guides/signal_processing`.

    References
    ----------
    J. F. Kaiser, "Nonrecursive digital filter design using the I0-sinh
    window function," *Proc. IEEE Int. Symp. Circuits Syst.* (1974) 20-23;
    Oppenheim and Schafer, *Discrete-Time Signal Processing*, 3rd ed.
    (2010), sec. 7.7.

    Examples
    --------
    >>> equiripple_order(0.01, 0.001, 50.0, 1000.0)
    52
    """
    _positive(transition_width=transition_width, fs=fs)
    _interval("passband_deviation", passband_deviation, 0.0, 1.0, closed=False)
    _interval("stopband_deviation", stopband_deviation, 0.0, 1.0, closed=False)
    if transition_width >= fs / 2:
        raise ValueError("transition_width must lie below fs/2")
    attenuation = -10 * np.log10(passband_deviation * stopband_deviation)
    return int(np.ceil((attenuation - 13) / (14.6 * transition_width / fs))) + 1


def _sections(sos):
    sos = np.asarray(sos, dtype=float)
    if sos.ndim != 2 or sos.shape[1] != 6 or not sos.shape[0]:
        raise ValueError("second-order sections must be an (n_sections, 6) array; pass a=None only with sections")
    _finite(sos=sos)
    return sos


def filter_signal(b, a, x):
    r"""Apply a digital filter, given as polynomials or as second-order sections.

    Parameters
    ----------
    b : array_like
        Numerator coefficients in :math:`z^{-1}` [-], or, when ``a`` is
        ``None``, an ``(n_sections, 6)`` array of second-order sections.
    a : array_like or None
        Denominator coefficients [-], or ``None`` for sections.
    x : array_like
        Input samples [signal units].

    Returns
    -------
    ndarray
        Output :math:`y`, with :math:`\sum_k a_k y_{n-k} = \sum_k b_k x_{n-k}`
        and zero initial state [signal units].

    Notes
    -----
    The filter runs forward once, so its start-up transient and phase delay
    remain. Sections are applied one after another (``scipy.signal.sosfilt``);
    see :func:`butterworth_lowpass` and :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> filter_signal([1, 1], [1], [1, 0, 0]).tolist()
    [1.0, 1.0, 0.0]
    """
    if a is None:
        return signal.sosfilt(_sections(b), np.asarray(x, dtype=np.result_type(x, np.float64)))
    b, a, x = np.asarray(b), np.asarray(a), np.asarray(x)
    dtype = np.result_type(b, a, x, np.float64)
    return signal.lfilter(b.astype(dtype), a.astype(dtype), x.astype(dtype))


def frequency_response(b, a=1, worN=512, fs=2 * np.pi):
    r"""Frequency response :math:`H(e^{j2\pi f/f_s})` of a digital filter.

    Parameters
    ----------
    b : array_like
        Numerator coefficients [-], or second-order sections when ``a`` is
        ``None``.
    a : array_like or None, optional
        Denominator coefficients [-]; 1 for an FIR filter, ``None`` for
        sections.
    worN : int or array_like, optional
        Number of frequencies from 0 to :math:`f_s/2`, or the frequencies
        themselves [units of ``fs``].
    fs : float, optional
        Sampling rate [Hz]; the default :math:`2\pi` gives frequencies in
        rad/sample.

    Returns
    -------
    tuple of ndarray
        ``(frequency, H)``: frequencies [units of ``fs``] and complex gain [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> f, H = frequency_response([0.5, 0.5], 1, worN=[0.0, 250.0], fs=1000.0)
    >>> np.round(np.abs(H), 6).tolist()
    [1.0, 0.707107]
    """
    if a is None:
        return signal.sosfreqz(_sections(b), worN=worN, fs=fs)
    return signal.freqz(b, a, worN=worN, fs=fs)
