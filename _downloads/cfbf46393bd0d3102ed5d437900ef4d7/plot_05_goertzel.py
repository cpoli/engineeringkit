r"""
Goertzel's algorithm
====================

Gerald Goertzel (1958) evaluated a single term of a Fourier series with a
second-order recursion that needs one real multiplication per sample
(:func:`~engineeringkit.signal_processing.goertzel`). When only a handful
of frequencies matter, it is cheaper than a full FFT and needs no buffer
of samples. Its classic use is the touch-tone (DTMF) telephone keypad,
where each key sends two tones, one from a low group and one from a high
group, and the receiver has to measure only eight frequencies.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.signal_processing import dft, goertzel

# %%
# One bin of the DFT
# ------------------
#
# At a frequency :math:`kf_s/N` the recursion returns exactly DFT bin
# :math:`k`.

rng = np.random.default_rng(1958)
fs, n = 8000.0, 205  # the classic DTMF block length
x = rng.normal(size=n)
for k in (3, 18, 40):
    assert np.isclose(goertzel(x, k * fs / n, fs), dft(x)[k])
print("Goertzel agrees with the DFT at bins 3, 18 and 40")

# %%
# Decoding a phone number
# -----------------------

low, high = [697.0, 770.0, 852.0, 941.0], [1209.0, 1336.0, 1477.0]
keys = ["123", "456", "789", "*0#"]
number = "5550123"
t = np.arange(n) / fs
decoded = ""
energies = []
for key in number:
    row = next(i for i, r in enumerate(keys) if key in r)
    column = keys[row].index(key)
    tone = np.sin(2 * np.pi * low[row] * t) + np.sin(2 * np.pi * high[column] * t) + rng.normal(0, 0.5, n)
    power = abs(goertzel(tone, low + high, fs)) ** 2
    energies.append(power)
    decoded += keys[int(np.argmax(power[:4]))][int(np.argmax(power[4:]))]
print(f"sent {number}, decoded {decoded}")
assert decoded == number

# %%
# Cost: eight Goertzel filters against one FFT
# --------------------------------------------
#
# Eight filters over 205 samples take about :math:`8 \times 205` real
# multiplications. A 256-point radix-2 FFT takes
# :math:`(N/2)\log_2 N = 1024` complex multiplications, about 4096 real ones,
# to compute every bin when only eight are needed.

goertzel_cost = 8 * n
fft_cost = 4 * 256 // 2 * 8
print(f"Goertzel: {goertzel_cost} real multiplications; FFT: {fft_cost}")
assert goertzel_cost < fft_cost

fig, ax = plt.subplots(figsize=(7, 3.5))
ax.imshow(np.array(energies).T, aspect="auto", cmap="viridis")
ax.set(xticks=range(len(number)), xticklabels=list(number), yticks=range(7), yticklabels=[f"{f:.0f} Hz" for f in low + high])
ax.set(xlabel="key pressed", title="Goertzel power at the eight DTMF frequencies")
fig.tight_layout()

plt.show()
