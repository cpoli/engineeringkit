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

Tests cover initial phase, constant/linear-message phase integration,
linear/log link-budget consistency, distance-squared attenuation, cascade
noise, AWGN BER and binary-symmetric-channel limits.

Reference: `SciPy cumulative trapezoid <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.cumulative_trapezoid.html>`_.
