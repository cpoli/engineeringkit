r"""
The Kalman-Bucy filter (1961)
=============================

A state estimator must trade trust in the model against trust in the
measurement. Rudolf Kalman and Richard Bucy (1961) found the gain that
minimises the steady estimation error for a linear system driven by white
process noise of spectral density :math:`Q` and measured through white noise
of density :math:`R`. The error covariance solves the Riccati equation of
the LQR transposed (:func:`~engineeringkit.control_systems.kalman_gain`).

For a random walk :math:`\dot x = w` measured as :math:`y = x + v`, the
optimal gain is :math:`L = \sqrt{Q/R}` and the error variance is
:math:`P = \sqrt{QR}`. Any other gain :math:`\ell` gives the larger variance
:math:`(Q + \ell^2 R)/(2\ell)`: too small a gain lags the state and too large
a gain passes on measurement noise.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import lfilter

from engineeringkit.control_systems import kalman_gain

q, r = 1.0, 0.04
L, P, estimator_poles = kalman_gain([[0.0]], [[1.0]], [[q]], [[r]])
L, P = float(L[0, 0]), float(P[0, 0])
print(f"L = {L:.4f} 1/s, P = {P:.4f}, estimator pole {estimator_poles[0].real:.4f} 1/s")
assert np.isclose(L, np.sqrt(q / r)) and np.isclose(P, np.sqrt(q * r))

# %%
# Simulation
# ----------
#
# Continuous white noise of density :math:`Q` sampled every :math:`\Delta t`
# has variance :math:`Q/\Delta t`. The estimator
# :math:`\dot{\hat x} = \ell(y - \hat x)` is stepped with Euler's method.

rng = np.random.default_rng(1961)
dt, n = 1e-3, 1_000_000
x = np.cumsum(rng.normal(0, np.sqrt(q * dt), n))
y = x + rng.normal(0, np.sqrt(r / dt), n)


def error_variance(gain):
    estimate = lfilter([0, gain * dt], [1, gain * dt - 1], y)
    return np.var((x - estimate)[10_000:]), estimate


gains = L * np.array([0.25, 0.5, 1.0, 2.0, 4.0])
measured = np.array([error_variance(g)[0] for g in gains])
theory = (q + gains**2 * r) / (2 * gains)
for g, m, th in zip(gains, measured, theory):
    print(f"gain {g:5.2f}: error variance {m:.4f} (theory {th:.4f})")
assert np.allclose(measured, theory, rtol=0.1)
assert np.argmin(measured) == 2  # the Kalman gain is best

_, estimate = error_variance(L)
k = slice(200_000, 210_000)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(dt * np.arange(n)[k], y[k], color="0.8", lw=0.5, label="measurement (unfiltered)")
ax1.plot(dt * np.arange(n)[k], x[k], "k-", label="state")
ax1.plot(dt * np.arange(n)[k], estimate[k], "C0", label="Kalman-Bucy estimate")
ax1.set(xlabel="time [s]", ylabel="x", ylim=(x[k].min() - 1, x[k].max() + 1))
ax1.legend()
ell = np.logspace(np.log10(gains[0]) - 0.2, np.log10(gains[-1]) + 0.2, 200)
ax2.loglog(ell, (q + ell**2 * r) / (2 * ell), "k-", lw=1, label=r"$(Q + \ell^2R)/2\ell$")
ax2.loglog(gains, measured, "o", label="simulated")
ax2.axvline(L, color="C3", ls="--", label=f"L = {L:.0f} 1/s")
ax2.set(xlabel=r"estimator gain $\ell$ [1/s]", ylabel="error variance")
ax2.legend()
fig.tight_layout()

plt.show()
