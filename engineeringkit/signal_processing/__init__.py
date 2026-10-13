"""Signals and systems from Fourier analysis to digital filtering and sampling."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive


def sampling_period(fs):
    """Sampling period.

    Parameters
    ----------
    fs : Hz
        Input in the units stated above.

    Returns
    -------
    result
        s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _positive(fs=fs)
    return 1 / fs


def nyquist_frequency(fs):
    """Nyquist frequency.

    Parameters
    ----------
    fs : Hz
        Input in the units stated above.

    Returns
    -------
    result
        Hz.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _positive(fs=fs)
    return fs / 2


def alias_frequency(f, fs):
    """Alias frequency.

    Parameters
    ----------
    f : Hz
        Input in the units stated above.
    fs : Hz
        Input in the units stated above.

    Returns
    -------
    result
        Hz.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _positive(fs=fs)
    return abs(((f + fs / 2) % fs) - fs / 2)


def db20(x):
    """Db20.

    Parameters
    ----------
    x : signal amplitude units
        Input in the units stated above.

    Returns
    -------
    result
        dB.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return 20 * np.log10(np.abs(x))


def rms(x):
    """Rms.

    Parameters
    ----------
    x : signal amplitude units
        Input in the units stated above.

    Returns
    -------
    result
        signal amplitude units.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    x = np.asarray(x)
    return np.sqrt(np.mean(np.abs(x) ** 2))


def energy(x):
    """Energy.

    Parameters
    ----------
    x : signal amplitude units
        Input in the units stated above.

    Returns
    -------
    result
        squared amplitude units (no dt weighting).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    x = np.asarray(x)
    return float(np.sum(np.abs(x) ** 2))


def power(x):
    """Power.

    Parameters
    ----------
    x : signal amplitude units
        Input in the units stated above.

    Returns
    -------
    result
        squared amplitude units (no dt weighting).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    x = np.asarray(x)
    return float(np.mean(np.abs(x) ** 2))


def convolution(x, h, mode="full"):
    """Convolution.

    Parameters
    ----------
    x : signal amplitude units
        Input in the units stated above.
    h : filter coefficient sequences (dimensionless)
        Input in the units stated above.
    mode : NumPy/SciPy option string
        Input in the units stated above.

    Returns
    -------
    result
        signal amplitude units.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return np.convolve(x, h, mode)


def correlation(x, y, mode="full"):
    """Correlation.

    Parameters
    ----------
    x : signal amplitude units
        Input in the units stated above.
    y : signal amplitude units
        Input in the units stated above.
    mode : NumPy/SciPy option string
        Input in the units stated above.

    Returns
    -------
    result
        signal amplitude units.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return signal.correlate(x, y, mode=mode)


def moving_average(x, n):
    """Moving average.

    Parameters
    ----------
    x : signal amplitude units
        Input in the units stated above.
    n : integer count
        Input in the units stated above.

    Returns
    -------
    result
        signal amplitude units.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _integer("n", n, 1)
    x = np.asarray(x)
    _finite(x=x)
    if x.ndim != 1 or n > len(x):
        raise ValueError("window n must not exceed the 1D signal length")
    return np.convolve(x, np.ones(n) / n, "valid")


def dft(x):
    """Dft.

    Parameters
    ----------
    x : signal amplitude units
        Input in the units stated above.

    Returns
    -------
    result
        signal amplitude units.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return np.fft.fft(x)


def idft(X):
    """Idft.

    Parameters
    ----------
    X : complex DFT coefficients, same amplitude units
        Input in the units stated above.

    Returns
    -------
    result
        signal amplitude units.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return np.fft.ifft(X)


def fft_frequencies(n, fs):
    """Fft frequencies.

    Parameters
    ----------
    n : integer count
        Input in the units stated above.
    fs : Hz
        Input in the units stated above.

    Returns
    -------
    result
        Hz.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _positive(fs=fs)
    return np.fft.fftfreq(n, 1 / fs)


def power_spectral_density(x, fs):
    """Power spectral density.

    Parameters
    ----------
    x : signal amplitude units
        Input in the units stated above.
    fs : Hz
        Input in the units stated above.

    Returns
    -------
    result
        (frequency [Hz], PSD [amplitude^2/Hz]).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _positive(fs=fs)
    return signal.periodogram(x, fs=fs)


def windowed_spectrum(x, fs, window="hann"):
    """Windowed spectrum.

    Parameters
    ----------
    x : signal amplitude units
        Input in the units stated above.
    fs : Hz
        Input in the units stated above.
    window : NumPy/SciPy option string
        Input in the units stated above.

    Returns
    -------
    result
        (frequency [Hz], unnormalized FFT magnitude).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _positive(fs=fs)
    x = np.asarray(x)
    w = signal.get_window(window, len(x))
    return np.fft.rfftfreq(len(x), 1 / fs), np.abs(np.fft.rfft(x * w))


def butterworth_lowpass(order, cutoff, fs):
    """Butterworth lowpass.

    Parameters
    ----------
    order : integer count
        Input in the units stated above.
    cutoff : Hz
        Input in the units stated above.
    fs : Hz
        Input in the units stated above.

    Returns
    -------
    result
        (numerator b, denominator a), dimensionless.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return signal.butter(order, cutoff, fs=fs, output="ba")


def butterworth_highpass(order, cutoff, fs):
    """Butterworth highpass.

    Parameters
    ----------
    order : integer count
        Input in the units stated above.
    cutoff : Hz
        Input in the units stated above.
    fs : Hz
        Input in the units stated above.

    Returns
    -------
    result
        (numerator b, denominator a), dimensionless.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return signal.butter(order, cutoff, btype="highpass", fs=fs, output="ba")


def fir_lowpass(numtaps, cutoff, fs, window="hamming"):
    """Fir lowpass.

    Parameters
    ----------
    numtaps : integer count
        Input in the units stated above.
    cutoff : Hz
        Input in the units stated above.
    fs : Hz
        Input in the units stated above.
    window : NumPy/SciPy option string
        Input in the units stated above.

    Returns
    -------
    result
        dimensionless FIR coefficients.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return signal.firwin(numtaps, cutoff, fs=fs, window=window)


def filter_signal(b, a, x):
    """Filter signal.

    Parameters
    ----------
    b : filter coefficient sequences (dimensionless)
        Input in the units stated above.
    a : filter coefficient sequences (dimensionless)
        Input in the units stated above.
    x : signal amplitude units
        Input in the units stated above.

    Returns
    -------
    result
        signal amplitude units.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    b, a, x = np.asarray(b), np.asarray(a), np.asarray(x)
    dtype = np.result_type(b, a, x, np.float64)
    return signal.lfilter(b.astype(dtype), a.astype(dtype), x.astype(dtype))


def frequency_response(b, a=1, worN=512, fs=2 * np.pi):
    """Frequency response.

    Parameters
    ----------
    b : filter coefficient sequences (dimensionless)
        Input in the units stated above.
    a : filter coefficient sequences (dimensionless)
        Input in the units stated above.
    worN : integer count or frequency array in fs units
        Input in the units stated above.
    fs : Hz
        Input in the units stated above.

    Returns
    -------
    result
        (frequencies in fs units, complex gain).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return signal.freqz(b, a, worN=worN, fs=fs)


def snr_db(signal_values, noise):
    """Snr db.

    Parameters
    ----------
    signal_values : signal amplitude units
        Input in the units stated above.
    noise : signal amplitude units
        Input in the units stated above.

    Returns
    -------
    result
        dB.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    values, noise = np.asarray(signal_values), np.asarray(noise)
    _finite(signal_values=values, noise=noise)
    if values.shape != noise.shape:
        raise ValueError("signal and noise must have matching shapes")
    sp, npower = np.mean(abs(values) ** 2), np.mean(abs(noise) ** 2)
    if sp == 0 and npower == 0:
        raise ValueError("SNR is undefined for zero signal and zero noise")
    if npower == 0:
        return float("inf")
    if sp == 0:
        return float("-inf")
    return 10 * np.log10(sp / npower)


def quantization_step(vmin, vmax, bits):
    """Step of a uniform mid-rise ADC with 2**bits bins over [vmin, vmax].

    Parameters
    ----------
    vmin : V (or consistent amplitude unit)
        Input in the units stated above.
    vmax : V (or consistent amplitude unit)
        Input in the units stated above.
    bits : integer count
        Input in the units stated above.

    Returns
    -------
    result
        signal amplitude units.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    if not isinstance(bits, (int, np.integer)) or isinstance(bits, bool) or not 1 <= bits <= 52:
        raise ValueError("bits must be an integer between 1 and 52")
    if not np.all(np.isfinite([vmin, vmax])) or vmax <= vmin:
        raise ValueError("finite vmax must exceed vmin")
    return (vmax - vmin) / (2**bits)


def quantize(x, vmin, vmax, bits):
    """Return mid-rise reconstruction levels, saturating outside ADC range.

        There are exactly 2**bits levels located at bin centers. Inputs and output
        share units (typically volts). This models ideal uniform quantization.

    Parameters
    ----------
    x : signal amplitude units
        Input in the units stated above.
    vmin : V (or consistent amplitude unit)
        Input in the units stated above.
    vmax : V (or consistent amplitude unit)
        Input in the units stated above.
    bits : integer count
        Input in the units stated above.

    Returns
    -------
    result
        signal amplitude units.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    q = quantization_step(vmin, vmax, bits)
    x = np.asarray(x, float)
    if not np.all(np.isfinite(x)):
        raise ValueError("samples must be finite")
    codes = np.clip(np.floor((x - vmin) / q), 0, 2**bits - 1)
    return vmin + (codes + 0.5) * q


def sinc_reconstruct(samples, fs, t):
    """Sinc reconstruct.

    Parameters
    ----------
    samples : signal amplitude units
        Input in the units stated above.
    fs : Hz
        Input in the units stated above.
    t : s
        Input in the units stated above.

    Returns
    -------
    result
        signal amplitude units.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/signal_processing`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _positive(fs=fs)
    samples = np.asarray(samples)
    tn = np.arange(len(samples)) / fs
    return np.sinc(fs * (np.asarray(t)[:, None] - tn)) @ samples


@dataclass(frozen=True)
class Spectrum:
    frequency: np.ndarray
    magnitude: np.ndarray

    @classmethod
    def from_signal(cls, x, fs):
        _positive(fs=fs)
        _finite(x=x)
        x = np.asarray(x)
        return cls(np.fft.rfftfreq(len(x), 1 / fs), np.abs(np.fft.rfft(x)))

    def plot(self):
        fig, ax = plt.subplots()
        ax.plot(self.frequency, self.magnitude)
        ax.set(xlabel="frequency [Hz]", ylabel="magnitude")
        return fig, ax
