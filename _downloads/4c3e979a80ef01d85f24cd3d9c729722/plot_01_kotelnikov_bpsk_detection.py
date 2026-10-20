r"""
Kotel'nikov's optimum receiver: BPSK in white Gaussian noise
============================================================

Vladimir Kotel'nikov's 1947 thesis, *The Theory of Optimum Noise Immunity*,
asked which receiver makes the fewest errors when known signals are buried in
white Gaussian noise. For two antipodal signals the answer is a correlator (or
matched filter) followed by a sign decision, with bit-error rate
:math:`Q(\sqrt{2E_b/N_0})` (:func:`~engineeringkit.communications.ber_bpsk_awgn`,
:func:`~engineeringkit.communications.q_function`). A Monte Carlo simulation of
:func:`~engineeringkit.communications.bpsk` through noise reproduces the curve.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.communications import ber_bpsk_awgn, bpsk, q_function
from engineeringkit.communications.visualizers import plot_ber_curve

# %%
# Monte Carlo against theory
# --------------------------

rng = np.random.default_rng(1947)
ebn0_db = np.arange(0, 9)
simulated = []
for x_db in ebn0_db:
    ebn0 = 10 ** (x_db / 10)
    bits = rng.integers(0, 2, 400_000)
    received = bpsk(bits) + rng.normal(0, np.sqrt(1 / (2 * ebn0)), bits.size)
    simulated.append(np.mean((received > 0) != bits))
theory = ber_bpsk_awgn(10 ** (ebn0_db / 10))
for x, s, th in zip(ebn0_db, simulated, theory, strict=True):
    print(f"Eb/N0 = {x} dB: simulated {s:.2e}, theory {th:.2e}")
assert np.allclose(simulated[:6], theory[:6], rtol=0.1)
assert np.allclose(theory, q_function(np.sqrt(2 * 10 ** (ebn0_db / 10))))

fine = np.linspace(0, 10, 200)
ax = plot_ber_curve(fine, ber_bpsk_awgn(10 ** (fine / 10)), label="Q(sqrt(2 Eb/N0))")
ax.semilogy(ebn0_db, simulated, "o", label="simulation")
ax.legend()
ax.set_title("Coherent BPSK in AWGN")

plt.show()
