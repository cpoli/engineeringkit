"""engineeringkit.signal_processing: sampling, Fourier analysis, spectra, filters, adaptive filters and wavelets.

Sampling, aliasing and sinc reconstruction; uniform quantization and SNR;
the DFT, a hand-written radix-2 FFT and Fourier partial sums; convolution
and correlation; the periodogram, windowed spectra, Welch's averaged
periodogram, the short-time Fourier transform and Goertzel's algorithm;
Butterworth IIR filters as polynomials or second-order sections,
window-method FIR filters with Kaiser's design formulas; the FIR Wiener
filter and the LMS adaptive filter; and the Haar wavelet transform. NumPy's
FFT and SciPy's filter design and filtering are used where they are the
standard tool; the algorithms named after their inventors are written out.
"""

from engineeringkit import __version__
from engineeringkit.signal_processing.core.base import LMSResult, Spectrum, STFTResult
from engineeringkit.signal_processing.systems.convolution import convolution, correlation, moving_average
from engineeringkit.signal_processing.systems.filters import (
    butterworth_highpass,
    butterworth_lowpass,
    filter_signal,
    fir_lowpass,
    frequency_response,
    kaiser_order,
)
from engineeringkit.signal_processing.systems.fourier import dft, fft_frequencies, fourier_partial_sum, idft, radix2_fft
from engineeringkit.signal_processing.systems.optimal import lms_filter, wiener_filter
from engineeringkit.signal_processing.systems.sampling import (
    alias_frequency,
    nyquist_frequency,
    quantization_step,
    quantize,
    sampling_period,
    sinc_reconstruct,
)
from engineeringkit.signal_processing.systems.spectral import goertzel, power_spectral_density, stft, welch_psd, windowed_spectrum
from engineeringkit.signal_processing.systems.wavelets import haar_transform, inverse_haar_transform
from engineeringkit.signal_processing.utils.measures import db20, energy, power, rms, snr_db

__all__ = [
    "__version__",
    "Spectrum",
    "STFTResult",
    "LMSResult",
    "db20",
    "rms",
    "energy",
    "power",
    "snr_db",
    "sampling_period",
    "nyquist_frequency",
    "alias_frequency",
    "sinc_reconstruct",
    "quantization_step",
    "quantize",
    "dft",
    "idft",
    "fft_frequencies",
    "radix2_fft",
    "fourier_partial_sum",
    "convolution",
    "correlation",
    "moving_average",
    "power_spectral_density",
    "windowed_spectrum",
    "welch_psd",
    "stft",
    "goertzel",
    "butterworth_lowpass",
    "butterworth_highpass",
    "fir_lowpass",
    "kaiser_order",
    "filter_signal",
    "frequency_response",
    "wiener_filter",
    "lms_filter",
    "haar_transform",
    "inverse_haar_transform",
]
