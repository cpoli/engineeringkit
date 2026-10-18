r"""
Lundberg and Palmgren's rolling-bearing life
============================================

Identical ball bearings run under identical loads fail at very different
times, because fatigue cracks start at randomly placed weak points below
the raceway. Gustaf Lundberg and Arvid Palmgren (1947) of SKF treated the
lives as Weibull-distributed, with a slope of about 10/9 for ball bearings.
They defined the rating life :math:`L_{10}`, which 90% of bearings reach,
and showed that it varies as :math:`(C/P)^p` with the basic load rating
:math:`C` and the applied load :math:`P`. The exponent is :math:`p = 3`
for balls and :math:`10/3` for rollers
(:func:`~engineeringkit.mechanical_engineering.bearing_l10_life`). The
example draws a simulated batch of bearings and checks both statements.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.mechanical_engineering import bearing_l10_life

# %%
# A batch of bearings
# -------------------
#
# A Weibull life of slope :math:`e` survives with probability
# :math:`R = 0.9^{(L/L_{10})^e}`. Draw 20 000 lives and count how many pass
# the rated life.

rng = np.random.default_rng(1947)
C, P, slope = 30e3, 5e3, 10 / 9
L10 = bearing_l10_life(C, P)
lives = L10 * (np.log(1 / rng.random(20_000)) / np.log(1 / 0.9)) ** (1 / slope)
survived = np.mean(lives >= L10)
median = np.median(lives)
print(f"L10 = {L10:.3e} rev; fraction surviving L10 = {survived:.4f}")
print(f"median life {median:.3e} rev = {median / L10:.2f} x L10 (Weibull: {(np.log(2) / np.log(1 / 0.9)) ** (1 / slope):.2f})")
assert abs(survived - 0.9) < 0.01

# %%
# The cube law
# ------------
#
# Doubling the load of a ball bearing cuts its life eightfold; a roller
# bearing loses a factor of :math:`2^{10/3} \approx 10`.

loads = np.geomspace(1e3, 30e3, 50)
ball, roller = bearing_l10_life(C, loads), bearing_l10_life(C, loads, p=10 / 3)
print(f"ball: L10(2P)/L10(P) = {bearing_l10_life(C, 2 * P) / L10:.4f}")
assert np.isclose(bearing_l10_life(C, 2 * P) / L10, 1 / 8)
assert np.allclose(np.diff(np.log(ball)) / np.diff(np.log(loads)), -3)

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
ordered = np.sort(lives)
axes[0].semilogx(ordered, 1 - np.arange(1, ordered.size + 1) / ordered.size, label="simulated batch")
axes[0].axvline(L10, color="k", ls=":", label="L10")
axes[0].axhline(0.9, color="0.6", ls="--")
axes[0].set(xlabel="revolutions", ylabel="fraction surviving", title="Scatter of fatigue lives")
axes[0].legend()
axes[1].loglog(loads / 1e3, ball, label="ball, p = 3")
axes[1].loglog(loads / 1e3, roller, label="roller, p = 10/3")
axes[1].set(xlabel="equivalent load P [kN]", ylabel="L10 [revolutions]", title=f"C = {C / 1e3:.0f} kN")
axes[1].legend()
fig.tight_layout()

plt.show()
