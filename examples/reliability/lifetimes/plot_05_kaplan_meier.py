r"""
Kaplan and Meier's product-limit estimator (1958)
=================================================

Field data are incomplete: some units are still working when the data are
collected, others were removed for reasons other than failure. Dropping
them biases the life estimate low; counting them as failures biases it too.
Kaplan and Meier (1958) multiplied the conditional survival
:math:`1 - d_j/n_j` at each failure, where :math:`n_j` counts only units
still under observation
(:func:`~engineeringkit.reliability.kaplan_meier`). The example censors a
Weibull sample at random and recovers the true survival curve.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.reliability import kaplan_meier, weibull_reliability

# %%
# Censored field data
# -------------------
#
# 200 pumps with Weibull lives (:math:`\eta` = 5000 h, :math:`\beta` =
# 1.8). Each enters service at a random date, so its observation ends at a
# uniformly distributed time up to 8000 h (time unit: hours).

rng = np.random.default_rng(1958)
eta, beta, n = 5000.0, 1.8, 200
life = eta * rng.weibull(beta, n)
window = rng.uniform(0, 8000, n)
times = np.minimum(life, window)
failed = (life <= window).astype(int)
km = kaplan_meier(times, failed)
print(km.summary())

# %%
# Compare with the truth and with naive estimates
# -----------------------------------------------
#
# "Ignore the censored units" treats only the failures as the sample.

t = np.linspace(0, 6000, 300)
naive = np.array([np.mean(times[failed == 1] > x) for x in t])
error = np.abs(km.survival_at(t) - weibull_reliability(t, eta, beta))
band = 2 * km.std_error
fig, ax = km.plot(label="Kaplan-Meier")
ax.fill_between(km.times, km.survival - band, km.survival + band, step="post", alpha=0.2, label=r"$\pm 2$ s.e. (Greenwood)")
ax.plot(t, weibull_reliability(t, eta, beta), "k--", label="true Weibull")
ax.plot(t, naive, color="C3", label="failures only (biased)")
ax.set(xlabel="time [h]", xlim=(0, 6000))
ax.legend()
fig.tight_layout()
naive_error = np.abs(naive - weibull_reliability(t, eta, beta))
print(f"largest error of Kaplan-Meier: {error.max():.3f}; of the naive estimate: {naive_error.max():.3f}")

assert error.max() < 0.08
assert naive_error.max() > 2 * error.max()

plt.show()
