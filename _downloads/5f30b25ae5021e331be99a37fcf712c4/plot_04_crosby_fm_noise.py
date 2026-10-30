r"""
Crosby's noise analysis of FM against AM (1937)
===============================================

Armstrong demonstrated wideband FM's quiet reception in 1935. Murray
Crosby of RCA explained it. Measured against the reference SNR
:math:`\gamma = P_R/N_0W` (:func:`~engineeringkit.communications.baseband_snr`),
an envelope-detected AM receiver gives at best :math:`\gamma/3`, because two
thirds of the power sits in the carrier
(:func:`~engineeringkit.communications.am_output_snr`). An FM discriminator
gives :math:`3\beta^2 P_m\gamma` (:func:`~engineeringkit.communications.fm_output_snr`):
every doubling of the deviation buys 6 dB. Crosby also found the catch.
Below a carrier-to-noise ratio of about 10 dB in the FM bandwidth, the
improvement collapses (:func:`~engineeringkit.communications.fm_threshold_snr`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfiltfilt

from engineeringkit.communications import am_output_snr, am_signal, carson_bandwidth, fm_output_snr, fm_threshold_snr

# %%
# Simulated receivers
# -------------------
#
# A 1 kHz tone in a message band W = 5 kHz, at complex baseband: the AM
# receiver takes the envelope, the FM receiver differentiates the phase.
# Each first filters the received signal to its channel (:math:`\pm W` for
# AM, twice the Carson band for FM, since Carson's band holds only 98% of
# the power and clipping the rest distorts), as an IF filter would, and low-pass filters
# the detector output to W. The FM peak deviation is :math:`\beta W`.

rng = np.random.default_rng(1937)
fs, W, fm_tone = 400e3, 5e3, 1e3
t = np.arange(int(0.2 * fs)) / fs
message = np.cos(2 * np.pi * fm_tone * t)
lowpass = butter(8, W, fs=fs, output="sos")
beta = 5.0
if_am = lowpass
if_fm = butter(8, carson_bandwidth(beta * W, W), fs=fs, output="sos")  # twice Carson: no sideband clipping


def output_snr(recovered, reference):
    recovered, reference = recovered[2000:-2000], reference[2000:-2000]
    recovered = recovered - recovered.mean()  # the AM detector also passes the carrier's DC
    gain = np.dot(recovered, reference) / np.dot(reference, reference)
    return np.mean((gain * reference) ** 2) / np.mean((recovered - gain * reference) ** 2)


def receive(signal, gamma):
    power = np.mean(np.abs(signal) ** 2)
    N0 = power / (gamma * W)
    noise = (rng.standard_normal(t.size) + 1j * rng.standard_normal(t.size)) * np.sqrt(N0 * fs / 2)
    return signal + noise


gammas_db = np.arange(10, 46, 5.0)
print(f"FM bandwidth by Carson's rule: {carson_bandwidth(beta * W, W) / 1e3:.0f} kHz for a {2 * W / 1e3:.0f} kHz AM channel")
am_measured, fm_measured = [], []
for g in 10 ** (gammas_db / 10):
    am_rx = receive(am_signal(message, t, 0.0) + 0j, g)
    am_measured.append(output_snr(sosfiltfilt(lowpass, np.abs(sosfiltfilt(if_am, am_rx))), message))
    phase = 2 * np.pi * beta * W * np.cumsum(message) / fs  # peak deviation beta * W
    fm_rx = receive(np.exp(1j * phase), g)
    frequency = np.diff(np.unwrap(np.angle(sosfiltfilt(if_fm, fm_rx))), prepend=0.0) * fs / (2 * np.pi)
    fm_measured.append(output_snr(sosfiltfilt(lowpass, frequency), message))
am_measured, fm_measured = 10 * np.log10(am_measured), 10 * np.log10(fm_measured)
am_theory = 10 * np.log10(am_output_snr(10 ** (gammas_db / 10)))
fm_theory = 10 * np.log10(fm_output_snr(10 ** (gammas_db / 10), beta))
threshold_db = 10 * np.log10(fm_threshold_snr(beta))
for g, a, f, ft in zip(gammas_db, am_measured, fm_measured, fm_theory):
    print(f"gamma {g:4.0f} dB: AM {a:5.1f} dB, FM {f:5.1f} dB (theory {ft:5.1f})")
above = gammas_db > threshold_db + 3
assert np.allclose(am_measured, am_theory, atol=1.0)
assert np.allclose(fm_measured[above], fm_theory[above], atol=1.5)
assert fm_measured[0] < fm_theory[0] - 5  # below threshold the FM advantage collapses

fig, ax = plt.subplots(figsize=(6.5, 4.2))
ax.plot(gammas_db, gammas_db, "k:", label=r"baseband, DSB-SC, SSB: $\gamma$")
ax.plot(gammas_db, am_theory, "C0", label=r"AM, $\mu = 1$: $\gamma/3$")
ax.plot(gammas_db, am_measured, "C0o")
ax.plot(gammas_db, fm_theory, "C1", label=rf"FM, $\beta$ = {beta:.0f}: $3\beta^2\gamma/2$")
ax.plot(gammas_db, fm_measured, "C1s")
ax.axvline(threshold_db, color="C3", ls="--", label="FM threshold (CNR 10 dB)")
ax.set(xlabel=r"reference SNR $\gamma$ [dB]", ylabel="output SNR [dB]")
ax.legend(fontsize=8)
fig.tight_layout()

plt.show()
