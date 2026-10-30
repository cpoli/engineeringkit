r"""
Nyquist's ISI criterion: raised-cosine pulses (1928)
====================================================

The same 1928 paper that set the limit of :math:`2B` pulses per second also
said which pulse shapes reach it without intersymbol interference: those
whose spectrum, copied at every multiple of the symbol rate, adds up to a
constant. The ideal sinc pulse meets the criterion with the least bandwidth
but decays as :math:`1/t`, so any timing error lets many neighbours leak in.
The raised-cosine family (:func:`~engineeringkit.communications.raised_cosine_pulse`,
:func:`~engineeringkit.communications.raised_cosine_spectrum`) spends an
excess bandwidth :math:`\alpha/2T` on a cosine roll-off that keeps the
criterion and makes the tails fall as :math:`1/t^3`. Split into two
root-raised-cosine filters (:func:`~engineeringkit.communications.root_raised_cosine_pulse`),
one at each end of the link, it also gives the receiver a matched filter.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.communications import (
    bpsk,
    raised_cosine_bandwidth,
    raised_cosine_pulse,
    raised_cosine_spectrum,
    root_raised_cosine_pulse,
)
from engineeringkit.communications.visualizers import plot_eye_diagram

# %%
# Zero crossings at every other symbol, and the folded spectrum
# --------------------------------------------------------------

T = 1.0
t = np.linspace(-6, 6, 2401)
f = np.linspace(-1.2, 1.2, 1201)
fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
for a in (0.0, 0.35, 1.0):
    axes[0].plot(t, raised_cosine_pulse(t, T, a), label=rf"$\alpha$ = {a}")
    axes[1].plot(f, raised_cosine_spectrum(f, T, a))
    k = np.arange(-6, 7)
    assert np.allclose(raised_cosine_pulse(k * T, T, a), k == 0, atol=1e-12)
    folded = sum(raised_cosine_spectrum(f + n / T, T, a) for n in range(-3, 4))
    assert np.allclose(folded[np.abs(np.abs(f) - 0.5) > 1e-9], T)
    print(f"alpha = {a}: bandwidth at 2400 baud = {float(raised_cosine_bandwidth(2400, a)):.0f} Hz")
axes[0].plot(np.arange(-6, 7), np.zeros(13), "k.", ms=4)
axes[0].set(xlabel="t / T", ylabel="p(t)", title="Raised-cosine pulses")
axes[0].legend()
axes[1].set(xlabel="f T", ylabel="P(f) / T", title="Spectra")

# %%
# Eye diagrams through a root-raised-cosine pair
# ----------------------------------------------
#
# Transmit and receive filters are root raised cosines, so the end-to-end
# pulse is a raised cosine and the eye closes to two points at the
# sampling instant, however small the roll-off. The small residue left is
# the cost of truncating each filter to eight symbols either side.

rng = np.random.default_rng(1928)
sps, span = 16, 8
taps = root_raised_cosine_pulse(np.arange(-span * sps, span * sps + 1) / sps, T, 0.35) / np.sqrt(sps)
symbols = bpsk(rng.integers(0, 2, 400))
upsampled = np.zeros(symbols.size * sps)
upsampled[::sps] = symbols
received = np.convolve(np.convolve(upsampled, taps), taps)
delay = 2 * span * sps
plot_eye_diagram(received, sps, ax=axes[2], offset=delay + sps // 2 + sps * 10)
samples = received[delay : delay + symbols.size * sps : sps]
isi = np.max(np.abs(samples[10:-10] - symbols[10:-10]))
print(f"worst intersymbol interference at the sampling instants: {isi:.1e}")
assert isi < 1e-2
axes[2].set(title=r"Eye after RRC $\times$ RRC, $\alpha$ = 0.35")
fig.tight_layout()

plt.show()
