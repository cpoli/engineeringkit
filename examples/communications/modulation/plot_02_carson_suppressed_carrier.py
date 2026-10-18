r"""
Carson's suppressed-carrier modulation (1915)
=============================================

In 1915 John Renshaw Carson of AT&T showed that the carrier in AM carries no
information. Multiplying the message by the carrier alone,
:math:`s(t) = m(t)\cos 2\pi f_c t`
(:func:`~engineeringkit.communications.dsb_sc`), puts all the transmitted
power into the sidebands. The receiver must supply a coherent carrier: mix,
low-pass filter, and the message returns, scaled by :math:`\tfrac12\cos\phi`
for a phase error :math:`\phi`. Carson also showed that one sideband is
enough, the basis of single-sideband telephony.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfiltfilt

from engineeringkit.communications import dsb_sc

# %%
# Coherent demodulation and phase error
# -------------------------------------

fs, fc = 200e3, 10e3
t = np.arange(0, 0.02, 1 / fs)
message = np.cos(2 * np.pi * 300 * t) + 0.5 * np.sin(2 * np.pi * 700 * t)
s = dsb_sc(message, t, fc)
lowpass = butter(6, 2e3, fs=fs, output="sos")
phases = np.radians([0, 45, 90])
fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
for phi in phases:
    recovered = 2 * sosfiltfilt(lowpass, s * np.cos(2 * np.pi * fc * t + phi))
    gain = np.dot(recovered, message) / np.dot(message, message)
    print(f"phase error {np.degrees(phi):3.0f} deg: recovered gain {gain:+.3f} (cos phi = {np.cos(phi):+.3f})")
    assert np.isclose(gain, np.cos(phi), atol=0.01)
    axes[0].plot(t * 1e3, recovered, label=f"φ = {np.degrees(phi):.0f}°")
axes[0].plot(t * 1e3, message, "k:", label="message")
axes[0].set(xlim=(0, 10), xlabel="time [ms]", title="Coherent detection")
axes[0].legend(fontsize=8)

# %%
# No carrier line in the spectrum
# -------------------------------

freq = np.fft.rfftfreq(t.size, 1 / fs)
axes[1].plot(freq / 1e3, np.abs(np.fft.rfft(s)) / t.size)
axes[1].set(xlim=(8.5, 11.5), xlabel="frequency [kHz]", title="DSB-SC spectrum: sidebands only")
carrier_bin = np.argmin(np.abs(freq - fc))
assert np.abs(np.fft.rfft(s))[carrier_bin] / t.size < 1e-10
fig.tight_layout()

plt.show()
