"""engineeringkit.signal_processing: sampling, Fourier analysis, spectra, filters, adaptive filters and wavelets.

Sampling, aliasing and sinc reconstruction; uniform quantization and SNR;
decimation and interpolation; the DFT, a hand-written radix-2 FFT, Fourier
partial sums and the DCT; linear and circular convolution, overlap-add and
correlation; the periodogram, windowed spectra, Welch's averaged
periodogram, the short-time Fourier transform, Goertzel's algorithm and the
cepstrum; autoregressive models by the Levinson-Durbin recursion; the
Hilbert transform and analytic signal; Butterworth, Chebyshev and elliptic
IIR filters as polynomials or second-order sections, the bilinear transform
with prewarping, window-method FIR filters with Kaiser's design formulas
and Parks-McClellan equiripple filters; the FIR Wiener filter and the LMS
adaptive filter; and the Haar wavelet transform. NumPy's FFT and SciPy's
filter design and filtering are used where they are the standard tool; the
algorithms named after their inventors are written out.
"""

from engineeringkit import __version__
from engineeringkit.signal_processing.core.base import ARModel, LMSResult, Spectrum, STFTResult
from engineeringkit.signal_processing.systems.autoregressive import levinson_durbin, yule_walker
from engineeringkit.signal_processing.systems.convolution import circular_convolution, convolution, correlation, moving_average, overlap_add
from engineeringkit.signal_processing.systems.filters import (
    bilinear_transform,
    butterworth_highpass,
    butterworth_lowpass,
    chebyshev1_lowpass,
    chebyshev2_lowpass,
    elliptic_lowpass,
    equiripple_lowpass,
    equiripple_order,
    filter_signal,
    fir_lowpass,
    frequency_response,
    iir_order,
    kaiser_order,
    prewarp_frequency,
)
from engineeringkit.signal_processing.systems.fourier import dct, dft, fft_frequencies, fourier_partial_sum, idct, idft, radix2_fft
from engineeringkit.signal_processing.systems.hilbert import (
    analytic_signal,
    envelope,
    hilbert_fir,
    hilbert_transform,
    instantaneous_frequency,
)
from engineeringkit.signal_processing.systems.optimal import lms_filter, wiener_filter
from engineeringkit.signal_processing.systems.sampling import (
    alias_frequency,
    decimate,
    downsample,
    interpolate,
    nyquist_frequency,
    quantization_step,
    quantize,
    sampling_period,
    sinc_reconstruct,
    upsample,
)
from engineeringkit.signal_processing.systems.spectral import (
    goertzel,
    power_spectral_density,
    real_cepstrum,
    stft,
    welch_psd,
    windowed_spectrum,
)
from engineeringkit.signal_processing.systems.wavelets import haar_transform, inverse_haar_transform
from engineeringkit.signal_processing.utils.measures import db20, energy, power, rms, snr_db

__all__ = [
    "__version__",
    "Spectrum",
    "STFTResult",
    "LMSResult",
    "ARModel",
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
    "downsample",
    "upsample",
    "decimate",
    "interpolate",
    "dft",
    "idft",
    "fft_frequencies",
    "radix2_fft",
    "fourier_partial_sum",
    "dct",
    "idct",
    "convolution",
    "circular_convolution",
    "overlap_add",
    "correlation",
    "moving_average",
    "power_spectral_density",
    "windowed_spectrum",
    "welch_psd",
    "stft",
    "goertzel",
    "real_cepstrum",
    "levinson_durbin",
    "yule_walker",
    "analytic_signal",
    "hilbert_transform",
    "envelope",
    "instantaneous_frequency",
    "hilbert_fir",
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
    "wiener_filter",
    "lms_filter",
    "haar_transform",
    "inverse_haar_transform",
]
