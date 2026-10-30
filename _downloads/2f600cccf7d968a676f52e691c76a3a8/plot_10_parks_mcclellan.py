r"""
Parks and McClellan's equiripple FIR filter (1972)
==================================================

A window design is easy but wasteful: its largest ripple sits next to the
band edge and the ripple dies away from it, so most of the band is better
than it needs to be. Thomas Parks and his student James McClellan at Rice
University posed linear-phase FIR design as Chebyshev approximation:
minimize the largest weighted error over the passband and stopband. The
alternation theorem says the optimum's error peaks with alternating signs
at least :math:`L + 2` times, and Remez's exchange algorithm finds those
frequencies (:func:`~engineeringkit.signal_processing.equiripple_lowpass`,
:func:`~engineeringkit.signal_processing.equiripple_order`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

from engineeringkit.signal_processing import db20, equiripple_lowpass, equiripple_order, fir_lowpass, frequency_response, kaiser_order

# %%
# Equiripple against Kaiser's window, same length
# -----------------------------------------------
#
# Pass 0-100 Hz, stop from 150 Hz, at 1 kHz, with 41 taps each.

fs, f_pass, f_stop, numtaps = 1000.0, 100.0, 150.0, 41
h_pm = equiripple_lowpass(numtaps, f_pass, f_stop, fs)
h_kaiser = fir_lowpass(numtaps, (f_pass + f_stop) / 2, fs, window=("kaiser", kaiser_order(50.0, 50.0, fs)[1]))
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
stopband = {}
for name, h in (("Parks-McClellan", h_pm), ("Kaiser window", h_kaiser)):
    f, H = frequency_response(h, 1, worN=8192, fs=fs)
    axes[0].plot(f, db20(np.maximum(np.abs(H), 1e-8)), label=name, lw=1)
    stopband[name] = -db20(np.abs(H[f >= f_stop]).max())
    deviation = np.abs(np.abs(H[f <= f_pass]) - 1).max()
    print(f"{name:16s}: stopband {stopband[name]:5.1f} dB, passband deviation {deviation:.4f}")
axes[0].axvspan(f_pass, f_stop, color="0.9")
axes[0].set(xlabel="frequency [Hz]", ylabel="gain [dB]", ylim=(-90, 5), title=f"{numtaps} taps")
axes[0].legend()

# %%
# The alternation theorem
# -----------------------
#
# For a symmetric filter of :math:`N = 2L + 1` taps the error
# :math:`E(f) = W(f)(|H| - D)` reaches its peak, alternating in sign, at
# least :math:`L + 2` times over the two bands.

grid = np.concatenate([np.linspace(0, f_pass, 4000), np.linspace(f_stop, fs / 2, 14000)])
f, H = frequency_response(h_pm, 1, worN=grid, fs=fs)
amplitude = np.real(H * np.exp(1j * 2 * np.pi * f / fs * (numtaps - 1) / 2))
error = np.where(f <= f_pass, amplitude - 1, amplitude)
pieces = [np.flatnonzero(f <= f_pass), np.flatnonzero(f >= f_stop)]
extrema = []
for idx in pieces:
    e = np.abs(error[idx])
    inner = 1 + np.flatnonzero((e[1:-1] >= e[:-2]) & (e[1:-1] >= e[2:]))
    extrema += [idx[i] for i in sorted({0, *inner, e.size - 1})]
peak = np.abs(error).max()
extrema = [i for i in extrema if np.abs(error[i]) > 0.99 * peak]
alternations = 1 + np.count_nonzero(np.diff(np.sign(error[extrema])))
for idx in pieces:
    axes[1].plot(f[idx], error[idx], color="C0")
axes[1].plot(f[extrema], error[extrema], "o", ms=3, color="C3")
axes[1].axhline(peak, color="0.5", lw=0.8)
axes[1].axhline(-peak, color="0.5", lw=0.8)
axes[1].set(xlabel="frequency [Hz]", ylabel="error", title=f"{alternations} alternations, L + 2 = {(numtaps - 1) // 2 + 2}")
fig.tight_layout()
print(f"alternations: {alternations}")

# %%
# Weights and length
# ------------------
#
# Weighting the stopband by :math:`\delta_p/\delta_s` sets the ratio of the
# deviations, and Kaiser's estimate gives the length for both at once. The
# estimate is empirical: here it lands about 16% above the target
# deviations, and a few more taps would close the gap.

delta_p, delta_s = 0.01, 0.001
n = equiripple_order(delta_p, delta_s, f_stop - f_pass, fs)
h = equiripple_lowpass(n, f_pass, f_stop, fs, weight=delta_p / delta_s)
f, H = frequency_response(h, 1, worN=8192, fs=fs)
dp, ds = np.abs(np.abs(H[f <= f_pass]) - 1).max(), np.abs(H[f >= f_stop]).max()
print(f"estimated {n} taps: deviations {dp:.4f} and {ds:.5f}")
assert np.isclose(dp / ds, delta_p / delta_s, rtol=0.02)
assert np.allclose(h, signal.remez(n, [0, f_pass, f_stop, fs / 2], [1, 0], weight=[1, 10], fs=fs))
assert alternations >= (numtaps - 1) // 2 + 2
assert stopband["Parks-McClellan"] > stopband["Kaiser window"]

plt.show()
