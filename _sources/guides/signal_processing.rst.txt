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
The BA IIR representation can be numerically fragile for high order; this
package does not yet expose second-order-section filtering.

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
