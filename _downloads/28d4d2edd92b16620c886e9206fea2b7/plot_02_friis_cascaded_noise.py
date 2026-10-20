r"""
Friis's noise formula for cascaded amplifiers
=============================================

Harald Friis (1944) defined a receiver's *noise factor* :math:`F` as the
ratio of input to output signal-to-noise ratio and showed that in a chain of
stages :math:`F = F_1 + (F_2-1)/G_1 + (F_3-1)/(G_1G_2) + \cdots`
(:func:`~engineeringkit.communications.cascaded_noise_factor`). Later stages
are suppressed by the gain ahead of them, so a low-noise amplifier belongs
first. The example checks the formula against a sample-by-sample simulation
of noisy amplifiers and compares two orderings of the same parts.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.communications import cascaded_noise_factor, from_db, noise_figure

# %%
# Simulate the chain
# ------------------
#
# Each stage adds input-referred noise of power :math:`(F-1)N_0`, then
# multiplies by its gain.

rng = np.random.default_rng(1944)
stages = [("LNA", from_db(20), from_db(1.0)), ("mixer", from_db(-7), from_db(8.0)), ("IF amp", from_db(30), from_db(4.0))]
gains = np.array([g for _, g, _ in stages])
factors = np.array([f for _, _, f in stages])
n0, n = 1.0, 400_000
noise = rng.normal(0, np.sqrt(n0), n)
for _, g, f in stages:
    noise = np.sqrt(g) * (noise + rng.normal(0, np.sqrt((f - 1) * n0), n))
simulated = np.var(noise) / (n0 * np.prod(gains))
formula = cascaded_noise_factor(gains, factors)
print(f"simulated F = {simulated:.3f}, Friis F = {formula:.3f} ({float(noise_figure(formula)):.2f} dB)")
assert np.isclose(simulated, formula, rtol=0.02)

# %%
# Put the low-noise amplifier first
# ---------------------------------

orders = {"LNA, mixer, IF": [0, 1, 2], "mixer, LNA, IF": [1, 0, 2], "mixer, IF, LNA": [1, 2, 0]}
figures = {name: float(noise_figure(cascaded_noise_factor(gains[i], factors[i]))) for name, i in orders.items()}
for name, nf in figures.items():
    print(f"{name:16s}: NF = {nf:5.2f} dB")
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.barh(list(figures), list(figures.values()))
ax.set(xlabel="system noise figure [dB]", title="Same parts, different order")
fig.tight_layout()

plt.show()
