r"""
Durbin's recursion for autoregressive fitting (1960)
====================================================

The Yule-Walker equations have a Toeplitz matrix, and Norman Levinson
(1947) had solved such systems, for Wiener's filter, by a recursion on the
order. James Durbin (1960) recast it for fitting autoregressions: the
order-:math:`m` model is the order-:math:`(m-1)` one plus a correction
weighted by a *reflection coefficient* :math:`k_m`. It costs :math:`p^2`
operations instead of :math:`p^3`, gives every order on the way, so the
prediction error can be watched as the order grows, and :math:`|k_m| < 1`
certifies a stable model
(:func:`~engineeringkit.signal_processing.levinson_durbin`). Linear
predictive coding of speech is built on it.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg, signal

from engineeringkit.signal_processing import levinson_durbin

# %%
# A fourth-order process
# ----------------------
#
# Two resonances. The biased sample autocorrelation of 20 000 samples
# feeds the recursion up to order 10.

poles = [0.95 * np.exp(1j * 0.6), 0.85 * np.exp(1j * 1.9)]
a_true = np.real(np.poly(poles + [p.conjugate() for p in poles]))
rng = np.random.default_rng(1960)
x = signal.lfilter([1.0], a_true, rng.normal(size=20000))
n = x.size
r = np.array([x[: n - k] @ x[k:] for k in range(11)]) / n
model = levinson_durbin(r, 10)

# %%
# The recursion solves the Toeplitz system
# ----------------------------------------

direct = linalg.solve(linalg.toeplitz(r[:10]), -r[1:11])
assert np.allclose(model.coefficients[1:], direct)
print("order 4 from the recursion:", np.round(levinson_durbin(r, 4).coefficients, 3))
print("true:                      ", np.round(a_true, 3))

# %%
# Prediction error against order
# ------------------------------
#
# :math:`E_m = (1 - k_m^2)E_{m-1}` falls until the true order and then
# stays flat, because the extra reflection coefficients are nearly zero.

k = model.reflection_coefficients
E = model.prediction_errors
fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
axes[0].plot(np.arange(11), E / E[0], "o-")
axes[0].axvline(4, color="0.6", ls="--")
axes[0].set(xlabel="order m", ylabel="$E_m/E_0$", yscale="log", title="Prediction error")
axes[1].stem(np.arange(1, 11), k, basefmt=" ")
axes[1].axhspan(-1, 1, color="0.95")
axes[1].set(xlabel="m", ylabel="$k_m$", ylim=(-1.1, 1.1), title="Reflection coefficients")
print("reflection coefficients:", np.round(k, 3))

# %%
# Cost
# ----
#
# The recursion takes about :math:`p^2` multiplications, Gaussian
# elimination about :math:`p^3/3`.

orders = np.array([4, 16, 64, 256])
axes[2].loglog(orders, orders**2, "o-", label="Levinson-Durbin, p²")
axes[2].loglog(orders, orders**3 / 3, "s-", label="elimination, p³/3")
axes[2].set(xlabel="order p", ylabel="multiplications", title="Cost")
axes[2].legend(fontsize=8)
fig.tight_layout()
assert np.all(np.abs(k) < 1)
assert np.abs(k[4:]).max() < 0.03 and abs(E[4] / E[-1] - 1) < 0.01
assert np.allclose(levinson_durbin(r, 4).coefficients, a_true, atol=0.02)

plt.show()
