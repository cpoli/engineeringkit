r"""
Bennett's quantization noise (1948)
===================================

Pulse-code modulation, worked out at Bell Laboratories during the Second
World War, rounds every sample to one of :math:`2^N` levels
(:func:`~engineeringkit.signal_processing.quantize`). W. R. Bennett asked
what that rounding does to the signal. For a signal that is busy compared
with the step :math:`q`
(:func:`~engineeringkit.signal_processing.quantization_step`), the error is
spread evenly over :math:`[-q/2, q/2]`, so it behaves like added white noise
of power :math:`q^2/12`. For a full-scale sine this gives the rule every ADC
data sheet quotes, :math:`\mathrm{SNR} = 6.02N + 1.76` dB: each extra bit
buys 6 dB.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.signal_processing import quantization_step, quantize, rms, snr_db

# %%
# The error is uniform with power q^2/12
# --------------------------------------
#
# A full-scale sine at an irrational fraction of the sampling rate, so the
# samples do not repeat, quantized to 8 bits over [-1, 1] V. The model is
# an approximation: here it is about 1.5% low in rms (0.14 dB).

n = np.arange(200_000)
x = np.sin(2 * np.pi * n * (np.sqrt(2) - 1) / 10) * (1 - 1e-9)
bits = 8
q = quantization_step(-1.0, 1.0, bits)
error = quantize(x, -1.0, 1.0, bits) - x
print(f"q = {q * 1e3:.3f} mV; rms error {rms(error) * 1e3:.4f} mV, q/sqrt(12) = {q / np.sqrt(12) * 1e3:.4f} mV")
assert np.abs(error).max() <= q / 2 + 1e-15
assert np.isclose(rms(error), q / np.sqrt(12), rtol=0.03)

# %%
# Six decibels per bit
# --------------------
#
# The rule is a large-N result: with only a few levels the error is not
# yet uniform, and the measured SNR falls short.

resolutions = np.arange(2, 17)
measured = [snr_db(x, quantize(x, -1.0, 1.0, int(b)) - x) for b in resolutions]
rule = 6.02 * resolutions + 1.76
for b, m in zip(resolutions[::4], measured[::4]):
    print(f"{b:2d} bits: SNR {m:6.2f} dB, 6.02N + 1.76 = {6.02 * b + 1.76:6.2f} dB")
assert np.allclose(measured[4:], rule[4:], atol=0.3)
assert measured[0] < rule[0] - 0.5

# %%
# When the signal is not busy
# ---------------------------
#
# A slow ramp spanning only a few steps gives a sawtooth error, correlated
# with the signal: Bennett's noise model needs many levels crossed.

ramp = np.linspace(-0.1, 0.1, 2000)
ramp_error = quantize(ramp, -1.0, 1.0, 5) - ramp

fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
axes[0].hist(error / q, bins=40, density=True)
axes[0].axhline(1.0, color="k", ls="--")
axes[0].set(xlabel="error / q", ylabel="density", title="8-bit error is uniform")
axes[1].plot(resolutions, measured, "o", label="measured")
axes[1].plot(resolutions, rule, label="6.02 N + 1.76 dB")
axes[1].set(xlabel="bits N", ylabel="SNR [dB]", title="Full-scale sine")
axes[1].legend()
axes[2].plot(ramp, ramp_error / quantization_step(-1.0, 1.0, 5))
axes[2].set(xlabel="input [V]", ylabel="error / q", title="5-bit ramp: a sawtooth, not noise")
fig.tight_layout()

plt.show()
