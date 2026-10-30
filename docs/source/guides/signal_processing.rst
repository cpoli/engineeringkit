Signal processing: sampled energy, spectra and quantization
===========================================================

Sample rate is Hz and time is seconds. Signals may be complex for DFT, RMS,
energy, power and SNR. Real FFT helpers require real data. energy is sum(abs(x)^2)
without sample-period weighting; it is a discrete signal measure, not joules.
For an approximation to a continuous integral multiply by 1/fs. power is the
mean squared amplitude and rms is its square root.

Spectral conventions
--------------------

DFT is the NumPy unnormalized forward transform; inverse includes 1/N.
windowed_spectrum and Spectrum return raw magnitudes, not calibrated one-sided
amplitudes. PSD uses SciPy periodogram density scaling, units amplitude^2/Hz.
frequency_response returns frequencies in the supplied fs units (default
2*pi means rad/sample). Digital cutoff parameters use Hz with fs in Hz.
The BA IIR representation loses precision above about order 8, or for a
cutoff far below fs/2; pass output="sos" to the Butterworth designs and give
the sections to filter_signal or frequency_response with a=None.

.. code-block:: python

    import numpy as np
    from engineeringkit.signal_processing import dft, energy
    x = np.exp(2j*np.pi*np.arange(16)/16)
    assert np.isclose(energy(x), energy(dft(x))/len(x))

Tests check Parseval, PSD integration, Butterworth cutoff gain 1/sqrt(2), FIR DC
gain and sinc reconstruction at the sample nodes. Truncated sinc reconstruction
is not an exact infinite bandlimited interpolator away from sample nodes.

ADC and SNR
-----------

Quantization is uniform mid-rise with 2**bits bins on [vmin,vmax]. Levels are
at bin centers, endpoints saturate, and the maximum in-range error is half a
step. Bits are restricted to 1..52 so floating-point level indexing is meaningful.
SNR uses mean squared magnitudes for complex values and matching signal/noise
shapes. Nonzero signal with zero noise returns +infinity; zero signal with
nonzero noise returns -infinity; zero/zero is rejected as undefined.

Fourier series, FFT and wavelets
--------------------------------

fourier_partial_sum evaluates a0 + sum(a_k cos(k w t) + b_k sin(k w t)) for the
coefficients given; it does not compute coefficients from a waveform.
radix2_fft is a recursive decimation-in-time FFT written for reading, not
speed, and needs a power-of-two length; its multiplication count includes the
trivial twiddle W^0. haar_transform is orthonormal (scaling 1/sqrt(2)), returns
[a_J, d_J, ..., d_1] like PyWavelets' wavedec, and needs a length divisible by
2**levels; no boundary extension is applied.

Spectral estimation
-------------------

welch_psd uses no detrending and matches scipy.signal.welch(..., detrend=False);
segments that do not fit after the last full one are dropped. stft does not pad
the record, so frame centres start L/2 samples in, and its coefficients are not
scaled. goertzel accepts any frequency, on or between bins, and returns the
unnormalized DFT sum at that frequency.

Kaiser design, Wiener and LMS
-----------------------------

kaiser_order uses Kaiser's empirical formulas, which usually land within about
a decibel of the requested attenuation; check the designed response.
wiener_filter solves the normal equations from biased sample correlations of
the whole record and assumes stationarity. lms_filter starts from zero weights
unless given others; the mean weights converge for 0 < mu < 2/lambda_max(R),
and a larger step trades excess error for speed.

IIR design: Chebyshev, elliptic and the bilinear transform
-----------------------------------------------------------

chebyshev1_lowpass and elliptic_lowpass take the passband edge as cutoff, where
the gain last equals -ripple_db; chebyshev2_lowpass takes the stopband edge,
where the gain first reaches -attenuation_db. Ripple is peak-to-peak in dB.
iir_order returns the smallest order meeting a low-pass specification for
"butterworth", "chebyshev" (both types) or "elliptic", with edges prewarped
as the bilinear designs need; it agrees with SciPy's buttord, cheb1ord and
ellipord. Butterworth's own cutoff is the -3 dB point, not the passband edge,
so place it between the edges. bilinear_transform takes analog polynomials in
descending powers of s and returns coefficients in z^-1 with a[0] = 1. With
K = 2 fs the analog axis is warped, f_a = (fs/pi) tan(pi f/fs); either design
the prototype at prewarp_frequency(f, fs) or pass prewarp=f so that f maps
exactly. Above order 8 or so use output="sos", as for Butterworth.

FIR design: Parks-McClellan
---------------------------

equiripple_lowpass calls scipy.signal.remez with bands [0, f_pass] and
[f_stop, fs/2]; the transition band is unconstrained. weight is the stopband
weight relative to the passband, so the designed deviations satisfy
delta_p = weight * delta_s. equiripple_order is Kaiser's empirical estimate
and can be a few taps short: check the design, as for kaiser_order. Remez can
fail to converge for very long filters or very narrow bands.

Convolution, zero padding and overlap-add
-----------------------------------------

circular_convolution zero-pads both inputs to n and wraps the linear result
modulo n; it is evaluated as an n x n product, for reading, not speed.
Padding to n >= len(x) + len(h) - 1 makes it the linear convolution.
overlap_add returns exactly convolution(x, h) ("full"), real for real inputs,
with FFT blocks of the next power of two above block_length + len(h) - 1.
Zero padding a record before a DFT (n_fft in windowed_spectrum and
real_cepstrum) samples the same spectrum more finely; it does not separate
components closer than about fs/N for an N-sample record.

Sample-rate conversion
----------------------

downsample and upsample do no filtering. decimate applies a Hamming FIR
low-pass with cutoff fs/(2M) (20M + 1 taps by default), removes its delay
and keeps every M-th sample; interpolate inserts zeros and applies a
windowed-sinc filter of gain L that leaves the original samples unchanged.
Both treat the record as zero outside its ends, so about numtaps/2 samples
at each end are transients. Both use integer factors; a rational factor L/M
is interpolate by L then decimate by M.

Hilbert transform, autoregressive models and the cepstrum
---------------------------------------------------------

analytic_signal matches scipy.signal.hilbert: the record is one period of a
periodic signal, so a record without a whole number of cycles leaks at its
ends. envelope and instantaneous_frequency are meaningful for a single
narrowband component. hilbert_fir is a windowed type III FIR with delay
(numtaps-1)/2 and zero gain at DC and Nyquist. levinson_durbin takes
autocorrelation lags r_0..r_p and raises ValueError if they are not positive
definite; ARModel coefficients follow A(z) = 1 + a_1 z^-1 + ... + a_p z^-p,
and ARModel.psd is one-sided, integrating to the model variance. yule_walker
uses the biased autocorrelation and removes the mean unless demean=False. The
order must be chosen: watch prediction_errors flatten. dct and idct are the
orthonormal DCT-II and DCT-III, written as matrix products, equal to
scipy.fft.dct(norm="ortho"). real_cepstrum raises ValueError when the
spectrum has an exact zero; its cepstrum is periodic in n_fft, so pad to
keep echoes from wrapping.

Tests and gallery examples check the FFT against NumPy, the Gibbs limit
Si(pi)/pi, Haar energy and inversion, Welch against SciPy, Goertzel against
DFT bins, the Wiener and LMS solutions against a known FIR system, Kaiser
designs against their specification, and second-order sections where the
polynomial form fails. They also check Chebyshev and elliptic ripple bounds
and minimum orders against SciPy, the bilinear transform against
scipy.signal.bilinear and the prewarped Butterworth, Parks-McClellan
alternation and weighting, circular convolution's time aliasing, overlap-add
against direct convolution, interpolation through the original samples and
decimation against aliasing, the Hilbert transform of a cosine and the
envelope of AM, the DCT against SciPy and its orthonormality, Levinson-Durbin
against a Toeplitz solve and Yule-Walker on a known AR(2) process, and the
cepstral echo peak alpha/2.
