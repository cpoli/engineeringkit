Communications: linear ratios, modulation phase and channel scope
=================================================================

Frequency/bandwidth are Hz, transmit/receive/noise power W, distance m and
absolute noise temperature K. Gains, SNR and Eb/N0 inputs are linear ratios,
not dB. db is 10*log10 of a positive power ratio; dbm uses 1 mW reference.
Zero logarithm inputs are rejected; use explicit -infinity for an ideal zero
power dB limit. Shannon capacity assumes AWGN and unrestricted ideal coding,
not an achievable throughput guarantee. Hartley/Nyquist helpers assume an
ideal bandlimited real baseband channel.

Link budgets and noise
----------------------

Friis assumes free space, far-field antennas, matched polarization/impedance
and consistent dimensionless gains. Positive distance alone does not establish
the far-field condition; no antenna aperture is provided for that check.
Cascaded noise factors use one positive gain and F>=1 per stage, with Friis's
noise formula Ftotal=F1+(F2-1)/G1+... . Thermal noise is k*T*B.

Modulation corrections
----------------------

BPSK/QPSK accept only binary data. Values are converted before arithmetic,
preventing unsigned 0 from wrapping to 255. QPSK returns unit-energy symbols
with real and imaginary signs determined by each bit pair; BER assumes coherent
Gray-labelled QPSK/BPSK in AWGN.

FM phase is 2*pi*fc*t+2*pi*kf*integral(message dt), with modulation phase zero
at the first supplied sample. Trapezoidal integration supports nonuniform time
grids. The old rectangular cumsum introduced a one-sample phase offset. AM
expects normalized message input; preventing overmodulation remains the
caller's responsibility. Neither model checks carrier Nyquist sampling.

.. code-block:: python

    import numpy as np
    from engineeringkit.communications import fm_signal
    t = np.linspace(0, .1, 101)
    waveform = fm_signal(np.ones_like(t)*2, t, fc=10, kf=3)
    # Analytic comparison: cos(2*pi*16*t).

Antennas
--------

effective_aperture is G*lambda^2/(4*pi) for a linear gain (not dBi);
aperture_gain is its inverse with an aperture efficiency in (0, 1]. Both use
the free-space wavelength c/f. directivity integrates an azimuthally
symmetric power pattern U(theta) sampled on increasing angles from 0 to pi
by the trapezoidal rule, so its accuracy is that of the grid. The
half-wave dipole helpers assume a thin wire with a sinusoidal current: the
pattern is cos(pi/2*cos(theta))/sin(theta), the directivity 4/Cin(2*pi) =
1.641 and the radiation resistance eta/(4*pi)*Cin(2*pi) = 73.08 ohm with
eta = 376.73 ohm. The feed reactance (about +42.5 ohm at exactly half a
wavelength) and wire thickness are not modelled.

Analog noise performance
------------------------

baseband_snr is gamma = P_R/(N0*W) with a one-sided N0. am_output_snr is
mu^2*Pm/(1+mu^2*Pm)*gamma for envelope detection and fm_output_snr is
3*beta^2*Pm*gamma for a limiter-discriminator, with beta = peak deviation / W
and Pm the mean square of a message normalised to abs(m) <= 1 (0.5 for a
tone). Both hold only above threshold: the envelope detector needs the
carrier well above the noise, and FM needs gamma above fm_threshold_snr =
2*(beta+1)*CNR_th (CNR_th = 10 by default). Pre-emphasis is not included.
DSB-SC and SSB with coherent detection give gamma itself.

Pulse shaping
-------------

raised_cosine_pulse has unit peak and zeros at every other multiple of the
symbol period T; raised_cosine_spectrum is its Fourier transform in seconds
(T in the flat band, T/2 at the edge when the roll-off is zero), and
root_raised_cosine_pulse has unit energy (units s^-1/2) so that it
convolved with itself is the raised cosine. Removable singularities take
their limits. The pulses are infinitely long: truncating them, as any FIR
filter must, leaves a small residual ISI. raised_cosine_bandwidth is the
baseband width (1+alpha)*Rs/2; a passband signal needs twice that.

Detection and constellations
----------------------------

matched_filter returns the time-reversed conjugate pulse; matched_filter_snr
is the peak output SNR 2E/N0 with a one-sided N0, valid for white Gaussian
noise and sampling at the pulse's end. psk_constellation (M a power of two)
and qam_constellation (M a power of four, square grids only) return
unit-average-energy points indexed by their Gray label. ser_mpsk_awgn
evaluates Craig's exact integral by adaptive quadrature, one value at a
time; ser_mqam_awgn is the exact square-QAM formula. Both take Es/N0 (not
Eb/N0; Es = Eb*log2 M) and assume coherent minimum-distance detection, as
minimum_distance_detect performs. Gray labelling makes the bit-error rate
close to SER/log2 M at moderate and high SNR.

Digitization
------------

pcm_encode is a mid-rise quantizer over [-V, V] with 2^N levels, saturating
outside the range, writing natural-binary level numbers MSB first;
pcm_decode returns interval centres, identical to
engineeringkit.signal_processing.quantize. pcm_sqnr_db is the large-N value
6.02N+1.76 dB for a full-scale sine; a simulation falls short by about
2.3*2^(-N/2) dB. mu_law_compress and mu_law_expand use the continuous
mu-law curve (mu = 255 by default), not the segmented G.711 table, and
reject samples beyond full scale. delta_modulate is linear (fixed-step)
delta modulation without leakage or step adaptation; it is free of slope
overload when the signal slope stays below delta_slope_limit = step*fs.

Source and channel coding
-------------------------

entropy and huffman_code take probabilities that sum to one (within 1e-9);
zero-probability symbols still get codewords. Huffman ties are broken by
creation order, so codewords are deterministic but other optimal codes with
the same lengths exist. CRC routines do plain polynomial division: register
preset to zero, no reflection, no final XOR (the XMODEM variant of
CRC-16-CCITT, check value 0x31C3); generators are integers that include the
leading term. convolutional_encode is a feedforward rate-1/n encoder whose
generators are tap integers, MSB for the current bit (octal (7, 5) is the
standard K = 3 code); it appends K-1 zero tail bits. viterbi_decode
requires that termination, uses the Hamming metric on hard bits or the
squared Euclidean metric on soft values with bit 1 sent as +1, and keeps
the whole survivor table (memory grows with message length times states).
convolutional_free_distance is found by Dijkstra's search over the trellis.

Fading and diversity
--------------------

rayleigh_pdf takes the mean-square envelope Omega. The BER formulas assume
flat, slow Rayleigh fading known at the receiver, coherent BPSK, and for
ber_bpsk_mrc_rayleigh and selection_outage_probability independent branches
of equal mean SNR. maximal_ratio_combine divides by the summed branch power,
so its output has unit gain on the symbol; a column whose gains are all zero
is rejected. Correlated branches, Rician fading and Doppler spectra are not
modelled.

Tests cover initial phase, constant/linear-message phase integration,
linear/log link-budget consistency, distance-squared attenuation, cascade
noise, AWGN BER and binary-symmetric-channel limits, the matched filter's
Cauchy-Schwarz bound, Nyquist's folded-spectrum criterion, PCM against the
signal_processing quantizer, Huffman optimality, the CRC-16 check value,
Viterbi against exhaustive maximum-likelihood search, M-PSK/QAM, Rayleigh
and diversity error rates by simulation, and the dipole's numerically
integrated directivity and radiation resistance.

Reference: `SciPy cumulative trapezoid <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.cumulative_trapezoid.html>`_.
