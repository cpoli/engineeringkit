r"""The Hilbert transform and the analytic signal: envelope, instantaneous frequency and the FIR Hilbert transformer.

Hilbert (1905) paired a function on the circle with its conjugate, the
boundary values of the real and imaginary parts of one analytic function.
On the frequency axis the pairing is a filter,
:math:`H(f) = -j\operatorname{sgn}(f)`, which shifts every component by
-90°. Gabor (1946) added the transform as an imaginary part to a real
signal, :math:`z = x + j\mathcal{H}x`, which has no negative frequencies:
:math:`|z|` is the envelope and the rate of change of :math:`\arg z` the
instantaneous frequency.
"""

from __future__ import annotations

import numpy as np
from scipy import signal

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive

__all__ = ["analytic_signal", "hilbert_transform", "envelope", "instantaneous_frequency", "hilbert_fir"]


def _real(x):
    x = np.asarray(x)
    _finite(x=x)
    if x.ndim != 1 or np.iscomplexobj(x) or x.size < 2:
        raise ValueError("x must be a real 1D signal with at least 2 samples")
    return x.astype(float)


def analytic_signal(x):
    r"""Analytic signal :math:`z = x + j\mathcal{H}x` by suppressing the negative-frequency half of the DFT.

    Parameters
    ----------
    x : array_like
        Real 1D samples [signal units], at least 2.

    Returns
    -------
    ndarray of complex
        :math:`z`, with :math:`\operatorname{Re} z = x` [signal units].

    Notes
    -----
    The DFT bins are multiplied by 1 at DC (and at Nyquist for even
    :math:`N`), 2 at positive frequencies and 0 at negative ones, then
    inverted. The record is treated as one period of a periodic signal, so
    a record that does not hold a whole number of cycles leaks at its ends.
    Equal to ``scipy.signal.hilbert``. See :doc:`/guides/signal_processing`.

    References
    ----------
    D. Gabor, "Theory of communication," *J. IEE* 93, part III (1946)
    429-457; S. L. Marple Jr., "Computing the discrete-time 'analytic'
    signal via FFT," *IEEE Trans. Signal Process.* 47 (1999) 2600-2603.

    Examples
    --------
    >>> z = analytic_signal(np.cos(2 * np.pi * np.arange(8) / 8))
    >>> (np.round(z.imag, 12) + 0.0).tolist()[:3]
    [0.0, 0.707106781187, 1.0]
    """
    x = _real(x)
    n = x.size
    weights = np.zeros(n)
    weights[0] = 1
    weights[1 : (n + 1) // 2] = 2
    if n % 2 == 0:
        weights[n // 2] = 1
    return np.fft.ifft(np.fft.fft(x) * weights)


def hilbert_transform(x):
    r"""Discrete Hilbert transform :math:`\mathcal{H}x`: every component shifted by -90°.

    Parameters
    ----------
    x : array_like
        Real 1D samples [signal units], at least 2.

    Returns
    -------
    ndarray
        :math:`\operatorname{Im} z` of :func:`analytic_signal` [signal units];
        :math:`\cos\omega n` becomes :math:`\sin\omega n`.

    Notes
    -----
    :math:`\mathcal{H}(\mathcal{H}x) = -x` for a record with no DC or
    Nyquist content, and :math:`x` and :math:`\mathcal{H}x` are orthogonal.
    See :doc:`/guides/signal_processing`.

    References
    ----------
    D. Hilbert, "Grundzüge einer allgemeinen Theorie der linearen
    Integralgleichungen (Dritte Mitteilung)," *Nachr. Ges. Wiss. Göttingen,
    Math.-Phys. Kl.* (1905) 307-338.

    Examples
    --------
    >>> (np.round(hilbert_transform(np.cos(2 * np.pi * np.arange(4) / 4)), 12) + 0.0).tolist()
    [0.0, 1.0, 0.0, -1.0]
    """
    return analytic_signal(x).imag


def envelope(x):
    r"""Envelope :math:`|z|` of a real signal, the magnitude of its analytic signal.

    Parameters
    ----------
    x : array_like
        Real 1D samples [signal units], at least 2.

    Returns
    -------
    ndarray
        Envelope [signal units].

    Notes
    -----
    For :math:`x = a_n\cos(\omega_c n + \phi_n)` with :math:`a_n` and
    :math:`\phi_n` varying slowly compared with the carrier (Bedrosian's
    condition), :math:`|z| = a_n`: an AM detector without a rectifier. See
    :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> np.round(envelope(3 * np.cos(2 * np.pi * np.arange(8) / 8)), 12).tolist()
    [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
    """
    return np.abs(analytic_signal(x))


def instantaneous_frequency(x, fs):
    r"""Instantaneous frequency :math:`\frac{1}{2\pi}\frac{d}{dt}\arg z` from the unwrapped analytic phase.

    Parameters
    ----------
    x : array_like
        Real 1D samples [signal units], at least 2.
    fs : float
        Sampling rate [Hz], positive.

    Returns
    -------
    ndarray
        ``len(x) - 1`` frequencies [Hz], each between two samples.

    Notes
    -----
    The phase difference of consecutive samples, times :math:`f_s/2\pi`.
    It is meaningful for a single component; two tones give a frequency
    that oscillates between them. See :doc:`/guides/signal_processing`.

    Examples
    --------
    >>> x = np.cos(2 * np.pi * 50 * np.arange(100) / 1000)
    >>> np.round(instantaneous_frequency(x, 1000.0)[40:42], 6).tolist()
    [50.0, 50.0]
    """
    _positive(fs=fs)
    return np.diff(np.unwrap(np.angle(analytic_signal(x)))) * fs / (2 * np.pi)


def hilbert_fir(numtaps, window="hamming"):
    r"""Windowed FIR Hilbert transformer, :math:`h_m = 2/(\pi m)` for odd :math:`m`, 0 for even.

    Parameters
    ----------
    numtaps : int
        Odd number of coefficients :math:`N`, at least 3; the taps sit at
        :math:`m = n - (N-1)/2`.
    window : str or tuple, optional
        Window, as in :func:`~engineeringkit.signal_processing.fir_lowpass`.

    Returns
    -------
    ndarray
        Antisymmetric impulse response [-], delay :math:`(N-1)/2` samples.

    Notes
    -----
    The ideal :math:`-j\operatorname{sgn}(\omega)` response has the
    impulse response :math:`(1 - \cos\pi m)/(\pi m)`, infinitely long. The
    truncated, windowed version is a type III linear-phase filter: its gain
    is zero at DC and Nyquist and close to 1 between. Delay :math:`x` by
    :math:`(N-1)/2` samples to pair it with the output. See
    :doc:`/guides/signal_processing`.

    References
    ----------
    A. V. Oppenheim and R. W. Schafer, *Discrete-Time Signal Processing*,
    3rd ed., Pearson (2010), sec. 12.4.

    Examples
    --------
    >>> np.round(hilbert_fir(5, window="boxcar") * np.pi, 12).tolist()
    [0.0, -2.0, 0.0, 2.0, 0.0]
    """
    _integer("numtaps", numtaps, 3)
    if numtaps % 2 == 0:
        raise ValueError("numtaps must be odd")
    m = np.arange(numtaps) - (numtaps - 1) // 2
    h = np.zeros(numtaps)
    odd = m % 2 == 1
    h[odd] = 2 / (np.pi * m[odd])
    return h * signal.get_window(window, numtaps, fftbins=False)
