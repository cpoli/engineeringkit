r"""
Widrow and Hoff's LMS adaptive filter (1960)
============================================

At Stanford in 1959, Bernard Widrow and his student Ted Hoff looked for a
way to train an adaptive linear neuron, ADALINE, without knowing the
statistics that the Wiener solution needs. Their least-mean-squares rule
nudges the weights after every sample against the gradient of the
instantaneous squared error: :math:`w \leftarrow w + \mu e_n u_n`
(:func:`~engineeringkit.signal_processing.lms_filter`). It needs only
:math:`2L` multiplications per sample, converges on average to the Wiener
filter, and follows it if the system changes. Echo cancellers, channel
equalizers and noise-cancelling headphones run on it.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.signal_processing import filter_signal, lms_filter, wiener_filter

# %%
# Identifying an unknown system
# -----------------------------
#
# White noise drives an unknown 8-tap FIR "plant" whose output is measured
# with a little noise. The LMS filter sees only the input and the measured
# output.

rng = np.random.default_rng(1960)
plant = np.array([0.6, -0.4, 0.25, 0.15, -0.1, 0.05, 0.03, -0.02])
n, taps = 6000, plant.size
x = rng.normal(size=n)
d = filter_signal(plant, [1.0], x) + rng.normal(0, 0.01, n)

result = lms_filter(x, d, taps, mu=0.01)
print("LMS weights:   ", np.round(result.weights, 3))
print("plant:         ", plant)
print("summary:", result.summary())
assert np.allclose(result.weights, plant, atol=0.01)

# %%
# It converges to the Wiener solution
# -----------------------------------

w_wiener = wiener_filter(x, d, taps)
assert np.allclose(result.weights, w_wiener, atol=0.01)

# %%
# Step size: speed against precision
# ----------------------------------
#
# For white input of unit power the weights converge when
# :math:`\mu < 2/L`, with a time constant of about :math:`1/\mu` samples.
# A larger step learns faster but leaves more excess error.

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
final = {}
for mu in (0.003, 0.01, 0.05):
    run = lms_filter(x, d, taps, mu)
    misalignment = np.sum((run.weight_history - plant) ** 2, axis=1) / np.sum(plant**2)
    final[mu] = run.summary()["final_mse"]
    axes[0].semilogy(misalignment, label=f"$\\mu$ = {mu}")
    print(f"mu = {mu:5.3f}: misalignment after 500 samples {misalignment[500]:.1e}, final MSE {final[mu]:.2e}")
assert final[0.003] < final[0.05]
axes[0].set(xlabel="sample", ylabel="weight misalignment", title="Learning curves")
axes[0].legend()

# %%
# Tracking a change
# -----------------
#
# Halfway through, the plant's first tap jumps from 0.6 to -0.3. The filter
# relearns it.

d_changed = d.copy()
d_changed[n // 2 :] += (-0.3 - 0.6) * x[n // 2 :]
tracking = lms_filter(x, d_changed, taps, mu=0.01)
assert abs(tracking.weights[0] + 0.3) < 0.01
tracking.plot(ax=axes[1])
axes[1].set(title="Squared error, plant changes at n = 3000")
fig.tight_layout()

plt.show()
