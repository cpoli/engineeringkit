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

Tests and gallery examples check the FFT against NumPy, the Gibbs limit
Si(pi)/pi, Haar energy and inversion, Welch against SciPy, Goertzel against
DFT bins, the Wiener and LMS solutions against a known FIR system, Kaiser
designs against their specification, and second-order sections where the
polynomial form fails.
