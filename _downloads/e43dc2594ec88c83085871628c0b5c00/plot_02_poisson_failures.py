r"""
Poisson's law of rare failures
==============================

Poisson (1837) derived the distribution of the number of rare events among
many trials, :math:`P(k) = (\lambda T)^k e^{-\lambda T}/k!`
(:func:`~engineeringkit.reliability.poisson_failures_probability`). A fleet
of units failing independently at a constant rate :math:`\lambda` has
:math:`\lambda T` expected failures in total operating time :math:`T`
(:func:`~engineeringkit.reliability.expected_failures`), and the count is
Poisson. The example simulates failures as exponential gaps and checks the
counts, then sizes a spares kit.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.reliability import expected_failures, poisson_failures_probability

# %%
# Simulate failure counts
# -----------------------
#
# 50 units, each failing at :math:`2\times10^{-5}` per hour, run for a year
# of 8760 h; every failed unit is replaced at once. Time unit: hours.

rng = np.random.default_rng(1837)
lam, units, hours = 2e-5, 50, 8760.0
mean = expected_failures(lam, units * hours)
trials = 20_000
gaps = rng.exponential(1 / lam, size=(trials, units, 8))
counts = (np.cumsum(gaps, axis=2) <= hours).sum(axis=(1, 2))
k = np.arange(0, 25)
pmf = np.array([poisson_failures_probability(int(i), lam, units * hours) for i in k])
observed = np.bincount(counts, minlength=k.size)[: k.size] / trials
print(f"expected failures per year: {mean:.2f}; simulated mean {counts.mean():.2f}")
assert np.isclose(pmf.sum(), 1, atol=1e-6)
assert np.abs(observed - pmf).max() < 0.01

# %%
# How many spares?
# ----------------
#
# Stock enough spares that the year's failures exceed them with
# probability under 1%.

cumulative = np.cumsum(pmf)
spares = int(np.argmax(cumulative >= 0.99))
print(f"spares for 99% cover: {spares} (P = {cumulative[spares]:.4f})")
fig, ax = plt.subplots(figsize=(6.5, 3.5))
ax.bar(k, observed, alpha=0.5, label="simulated fleet")
ax.plot(k, pmf, "o-", color="C1", label="Poisson")
ax.axvline(spares + 0.5, color="gray", ls="--", label=f"{spares} spares")
ax.set(xlabel="failures in one year", ylabel="probability")
ax.legend()
fig.tight_layout()

assert cumulative[spares] >= 0.99 > cumulative[spares - 1]

plt.show()
