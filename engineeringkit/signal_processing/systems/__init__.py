"""Concrete signal-processing models: sampling, Fourier analysis, convolution, spectra, filters, optimal filters, wavelets."""

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

__all__ = [
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
