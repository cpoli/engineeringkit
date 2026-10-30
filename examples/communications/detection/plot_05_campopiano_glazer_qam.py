r"""
Campopiano and Glazer's square QAM (1962)
=========================================

Phase-shift keying wastes the inside of the circle. C. N. Campopiano and
B. G. Glazer modulated the amplitudes of two carriers in quadrature
independently, so the points fill a square grid
(:func:`~engineeringkit.communications.qam_constellation`). Each rail is a
:math:`\sqrt M`-level amplitude decision, and the exact symbol-error rate
follows from one rail's error probability
(:func:`~engineeringkit.communications.ser_mqam_awgn`). At 16 points and
above, square QAM needs several decibels less power than M-PSK, which is
why every modem from the 1970s to Wi-Fi and 5G uses it.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.communications import from_db, minimum_distance_detect, qam_constellation, ser_mpsk_awgn, ser_mqam_awgn

# %%
# Constellation and simulated error rate
# --------------------------------------

rng = np.random.default_rng(1962)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
points = qam_constellation(16)
axes[0].plot(points.real, points.imag, "o")
for label, point in enumerate(points):
    axes[0].annotate(f"{label:04b}", (point.real, point.imag), xytext=(0, 9), textcoords="offset points", ha="center", fontsize=8)
axes[0].set(aspect="equal", title="16-QAM, Gray labelled", xlabel="in-phase", ylabel="quadrature")
assert np.isclose(np.mean(np.abs(points) ** 2), 1)

ebn0_db = np.arange(0, 19, 2.0)
for M in (4, 16, 64):
    esn0 = from_db(ebn0_db) * np.log2(M)
    theory = ser_mqam_awgn(esn0, M)
    axes[1].semilogy(ebn0_db, theory, label=f"{M}-QAM")
    pts = qam_constellation(M)
    n = 40_000
    labels = rng.integers(0, M, n)
    measured = []
    for g in esn0:
        noise = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) * np.sqrt(1 / (2 * g))
        measured.append(np.mean(minimum_distance_detect(pts[labels] + noise, pts) != labels))
    measured = np.array(measured)
    axes[1].semilogy(ebn0_db[measured > 0], measured[measured > 0], "k.", ms=4)
    assert np.allclose(measured[theory > 1e-2], theory[theory > 1e-2], rtol=0.1)
axes[1].semilogy(ebn0_db, ser_mpsk_awgn(from_db(ebn0_db) * 4, 16), "C1--", label="16-PSK")
axes[1].set(xlabel="$E_b/N_0$ [dB]", ylabel="symbol-error rate", ylim=(1e-5, 1), title="QAM in AWGN (dots: simulation)")
axes[1].grid(True, which="both", alpha=0.3)
axes[1].legend()
fig.tight_layout()

# %%
# QAM against PSK at 16 points
# ----------------------------
#
# The minimum distance of unit-energy 16-QAM is :math:`2\sqrt{1/10} = 0.632`;
# that of 16-PSK is :math:`2\sin(\pi/16) = 0.390`, a ratio of 4.2 dB.

gap = 20 * np.log10(2 * np.sqrt(0.1) / (2 * np.sin(np.pi / 16)))
target = 1e-4
grid = np.linspace(10, 30, 2001)
need_qam = grid[np.argmax(ser_mqam_awgn(from_db(grid), 16) < target)]
need_psk = grid[np.argmax(ser_mpsk_awgn(from_db(grid), 16) < target)]
print(f"distance gap {gap:.2f} dB; Es/N0 for SER 1e-4: 16-QAM {need_qam:.2f} dB, 16-PSK {need_psk:.2f} dB")
assert abs((need_psk - need_qam) - gap) < 0.3

plt.show()
