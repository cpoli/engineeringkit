r"""
Wiener's optimal filter (1949)
==============================

Asked during the war to predict the path of an aircraft from noisy radar
data, Norbert Wiener posed filtering as an optimization: find the linear
filter that makes the mean squared error between its output and the
wanted signal as small as possible. His 1942 report, bound in yellow and
known to engineers as the "yellow peril" for its mathematics, was published
in 1949. Norman Levinson (1947) recast it for discrete FIR filters, where
the answer is the solution of the normal equations :math:`Rw = p`
(:func:`~engineeringkit.signal_processing.wiener_filter`). The example
recovers a slowly varying signal from additive white noise.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.signal_processing import filter_signal, wiener_filter

# %%
# Signal and noise
# ----------------
#
# The wanted signal :math:`s` is an AR(1) process,
# :math:`s_n = a s_{n-1} + e_n`, observed as :math:`x = s + v` with white
# noise :math:`v` of the same power.

rng = np.random.default_rng(1949)
a, n = 0.95, 200_000
s = filter_signal([np.sqrt(1 - a**2)], [1.0, -a], rng.normal(size=n))  # unit variance
x = s + rng.normal(0, 1.0, n)

# %%
# The optimal FIR filter
# ----------------------

taps = 30
w = wiener_filter(x, s, taps)
estimate = filter_signal(w, [1.0], x)
mse = np.mean((estimate - s) ** 2)
print(f"noisy input MSE {np.mean((x - s) ** 2):.3f}, Wiener output MSE {mse:.3f}")

# %%
# Against theory
# --------------
#
# Here the correlations are known exactly: :math:`r_{xx}(k) = a^{|k|} + \delta_k`
# and :math:`r_{sx}(k) = a^{|k|}`. Solving :math:`Rw = p` with them gives the
# theoretical weights and the minimum error :math:`1 - p^T w`. The filter
# estimated from data reaches it.

k = np.arange(taps)
R = a ** abs(k[:, None] - k[None, :]) + np.eye(taps)
p = a**k
w_theory = np.linalg.solve(R, p)
mse_theory = 1 - p @ w_theory
print(f"theoretical minimum MSE {mse_theory:.3f}")
assert np.allclose(w, w_theory, atol=0.01)
assert abs(mse - mse_theory) < 0.01

# %%
# No other filter does better
# ---------------------------
#
# A moving average of any length, the obvious smoother, has a larger
# error than the Wiener filter.

averages = {M: np.mean((filter_signal(np.ones(M) / M, [1.0], x) - s) ** 2) for M in (2, 4, 8, 16, 30)}
for M, err in averages.items():
    print(f"moving average of {M:2d}: MSE {err:.3f}")
assert min(averages.values()) > mse

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
window = slice(1000, 1300)
axes[0].plot(x[window], color="0.75", lw=0.8, label="observed")
axes[0].plot(s[window], "k", lw=1.5, label="signal")
axes[0].plot(estimate[window], label="Wiener estimate")
axes[0].set(xlabel="sample")
axes[0].legend(fontsize=8)
axes[1].stem(k, w, basefmt=" ", label="from data")
axes[1].plot(k, w_theory, "r.", label="theory")
axes[1].set(xlabel="tap", ylabel="weight", title=f"{taps}-tap Wiener filter")
axes[1].legend()
fig.tight_layout()

plt.show()
