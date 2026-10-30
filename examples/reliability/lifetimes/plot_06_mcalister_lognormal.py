r"""
McAlister's lognormal law (1879)
================================

Francis Galton asked why the geometric mean of a set of observations is
sometimes a better centre than the arithmetic one, and Donald McAlister
(1879) answered with the law now called lognormal: if a quantity is the
product of many small independent factors, its logarithm is a sum, and so
is normally distributed. Fatigue cracks, corrosion pits and the
electromigration of metal lines grow in such multiplicative steps, which is
why the lognormal is the life model for these mechanisms
(:func:`~engineeringkit.reliability.lognormal_reliability`). Its hazard
(:func:`~engineeringkit.reliability.lognormal_hazard`) rises to a peak and
then falls, unlike a Weibull hazard with :math:`\beta > 1`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from engineeringkit.reliability import lognormal_cdf, lognormal_hazard, lognormal_mttf, weibull_hazard

# %%
# A product of many factors gives a lognormal life
# ------------------------------------------------
#
# The fatigue life of each of 20000 parts is a nominal 1000 h multiplied by
# 30 independent factors, each between 0.8 and 1.25 (surface finish,
# inclusions, local stress, residual stress, ...). No factor is lognormal,
# but their product is close to it, as McAlister argued.

rng = np.random.default_rng(1879)
life = 1000.0 * np.prod(rng.uniform(0.8, 1.25, size=(20000, 30)), axis=1)
mu, sigma = np.log(life).mean(), np.log(life).std()
print(f"ln(life): mean {mu:.3f}, standard deviation {sigma:.3f}")
print(f"MTTF: sample {life.mean():.0f} h, lognormal formula {lognormal_mttf(mu, sigma):.0f} h")

t = np.linspace(0, 4000, 200)
fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
axes[0].plot(np.sort(life), np.arange(1, life.size + 1) / life.size, label="simulated")
axes[0].plot(t, lognormal_cdf(t, mu, sigma), "--", label="lognormal")
axes[0].set(xlabel="life [h]", ylabel="failure probability", xlim=(0, 4000))
axes[0].legend()

# %%
# A hazard that rises and falls
# -----------------------------

h_time = np.linspace(1, 2000, 400)
axes[1].plot(h_time, lognormal_hazard(h_time, np.log(500), 0.8), label=r"lognormal, $\sigma$ = 0.8")
axes[1].plot(h_time, weibull_hazard(h_time, 650, 1.6), label=r"Weibull, $\beta$ = 1.6")
axes[1].set(xlabel="time [h]", ylabel="hazard [1/h]")
axes[1].legend()
fig.tight_layout()

peak = h_time[np.argmax(lognormal_hazard(h_time, np.log(500), 0.8))]
print(f"lognormal hazard peaks at {peak:.0f} h and then falls")

assert stats.kstest(np.log(life), stats.norm(mu, sigma).cdf).statistic < 0.02
assert np.isclose(lognormal_mttf(mu, sigma), life.mean(), rtol=0.01)
assert 100 < peak < 1900 and lognormal_hazard(2000.0, np.log(500), 0.8) < lognormal_hazard(peak, np.log(500), 0.8)

plt.show()
