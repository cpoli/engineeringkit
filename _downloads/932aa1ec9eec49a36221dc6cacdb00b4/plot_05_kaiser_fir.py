r"""
Kaiser's window method for FIR design (1966)
============================================

An ideal low-pass filter has an infinitely long sinc impulse response.
Truncating it gives a finite (FIR) filter whose response ripples near the
band edge. Tapering the truncated response with a window trades ripple
for a wider transition band. James Kaiser of Bell Labs (1966) used the
prolate-spheroidal-like :math:`I_0`-Bessel window, whose single parameter
:math:`\beta` sets that trade continuously. He later (1974) gave empirical
formulas for the length and :math:`\beta` that meet a required stopband
attenuation and transition width
(:func:`~engineeringkit.signal_processing.kaiser_order`,
:func:`~engineeringkit.signal_processing.fir_lowpass`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.signal_processing import db20, fir_lowpass, frequency_response, kaiser_order

# %%
# Design to a specification
# -------------------------
#
# Pass up to 100 Hz, stop from 150 Hz, at 1 kHz sampling, with stopband
# attenuations from 30 to 90 dB. The cutoff is placed mid-transition.

fs, f_pass, f_stop = 1000.0, 100.0, 150.0
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
achieved = {}
for attenuation in (30.0, 50.0, 70.0, 90.0):
    numtaps, beta = kaiser_order(attenuation, f_stop - f_pass, fs)
    h = fir_lowpass(numtaps, (f_pass + f_stop) / 2, fs, window=("kaiser", beta))
    f, H = frequency_response(h, 1, worN=8192, fs=fs)
    stop = -db20(abs(H[f >= f_stop]).max())
    ripple = abs(H[f <= f_pass]).max() - 1
    achieved[attenuation] = stop
    print(f"A = {attenuation:2.0f} dB: {numtaps:3d} taps, beta = {beta:5.2f}, stopband {stop:5.1f} dB, passband ripple {ripple:.1e}")
    axes[0].plot(f, db20(np.maximum(abs(H), 1e-8)), lw=0.8, label=f"A = {attenuation:.0f} dB, N = {numtaps}")
    # Kaiser's formulas are empirical: they land within about a decibel.
    assert stop > attenuation - 1.5
axes[0].axvspan(f_pass, f_stop, color="0.9")
axes[0].set(xlabel="frequency [Hz]", ylabel="gain [dB]", ylim=(-120, 5), title="Kaiser designs")
axes[0].legend(fontsize=8)

# %%
# Linear phase
# ------------
#
# The coefficients are symmetric, so every frequency is delayed by
# :math:`(N-1)/2` samples and the waveform's shape survives.

numtaps, beta = kaiser_order(60.0, f_stop - f_pass, fs)
h = fir_lowpass(numtaps, 125.0, fs, window=("kaiser", beta))
assert np.allclose(h, h[::-1])
f, H = frequency_response(h, 1, worN=np.linspace(1, 90, 50), fs=fs)
delay = -np.gradient(np.unwrap(np.angle(H)), 2 * np.pi * f / fs)
assert np.allclose(delay, (numtaps - 1) / 2)
print(f"group delay {delay.mean():.1f} samples = (N - 1)/2 for N = {numtaps}")

# %%
# Against a fixed window
# ----------------------
#
# A Hamming window always gives about 53 dB, whatever the length; the
# Kaiser window reaches any attenuation by choosing :math:`\beta`.

for numtaps in (41, 81):
    f, H = frequency_response(fir_lowpass(numtaps, 125.0, fs, window="hamming"), 1, worN=8192, fs=fs)
    width = 3.3 * fs / numtaps
    print(f"Hamming, {numtaps} taps: stopband {-db20(abs(H[f >= 125.0 + width / 2]).max()):.1f} dB")
axes[1].stem(h, basefmt=" ")
axes[1].set(xlabel="tap", ylabel="coefficient", title=f"60 dB design: {h.size} taps")
fig.tight_layout()

plt.show()
