r"""
Frahm's vibration absorber and Den Hartog's tuning (1909)
=========================================================

Hermann Frahm (patent filed 1909) stopped ships' hulls and machines from
vibrating at their running speed by attaching a small spring-mass system
tuned to that speed. At the tuned frequency the absorber swings so that
its spring force cancels the applied force, and the machine stands still.
The cure works at one frequency only: the system now has two resonances, one
on each side. Jacob Pieter Den Hartog and John Ormondroyd (1928) added
damping to the absorber. Every response curve then passes through two
fixed points, whatever the damping, and the best design tunes the absorber
to make them equal and puts the peaks there
(:func:`~engineeringkit.mechanical_engineering.tuned_mass_damper`). This
is the tuned mass damper of tall buildings.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.mechanical_engineering import tuned_mass_damper

# %%
# Frahm's undamped absorber
# -------------------------
#
# Tuned to the primary frequency (:math:`f = 1`) with mass ratio
# :math:`\mu`, the primary mass stops at :math:`r = 1`. The two new
# resonances are at :math:`r^2 = 1 + \mu/2 \pm \sqrt{\mu + \mu^2/4}`.

m, k, mu = 500.0, 2e6, 0.1
frahm = tuned_mass_damper(m, k, mu, frequency_ratio=1.0, damping_ratio=0.0)
split = np.sqrt(1 + mu / 2 + np.array([-1, 1]) * np.sqrt(mu + mu**2 / 4))
print(f"amplification at r = 1: {frahm.amplification(1.0):.2e}; new resonances at r = {np.round(split, 4).tolist()}")
print(f"absorber stroke at r = 1: {abs(frahm.response(1.0)[1]):.1f} x static deflection")
assert frahm.amplification(1.0) < 1e-12
assert np.all(frahm.amplification(split) > 1e9)

# %%
# Den Hartog's fixed points
# -------------------------
#
# With :math:`f = 1/(1 + \mu)` the fixed points are at
# :math:`r^2 = (1 \pm \sqrt{\mu/(2 + \mu)})/(1 + \mu)` and both have height
# :math:`\sqrt{1 + 2/\mu}`. The optimum damping puts the peaks on them.

fixed = np.sqrt((1 + np.array([-1, 1]) * np.sqrt(mu / (2 + mu))) / (1 + mu))
r = np.linspace(0.6, 1.4, 4001)
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
for zeta in (0.02, 0.1, 0.4):
    design = tuned_mass_damper(m, k, mu, damping_ratio=zeta)
    assert np.allclose(design.amplification(fixed), np.sqrt(1 + 2 / mu), rtol=1e-9)
    axes[0].semilogy(r, design.amplification(r), label=f"zeta = {zeta}")
optimum = tuned_mass_damper(m, k, mu)
peak = optimum.amplification(r).max()
print(f"optimum: f = {optimum.frequency_ratio:.4f}, zeta = {optimum.damping_ratio:.4f}")
print(f"peak {peak:.3f}, fixed-point height sqrt(1 + 2/mu) = {np.sqrt(1 + 2 / mu):.3f}")
assert np.isclose(peak, np.sqrt(1 + 2 / mu), rtol=0.01)
axes[0].semilogy(r, optimum.amplification(r), "k", lw=2, label="Den Hartog optimum")
axes[0].plot(fixed, optimum.amplification(fixed), "ko")
axes[0].set(xlabel="forcing / primary frequency", ylabel="amplification", title=f"mu = {mu}, f = 1/(1 + mu)")
axes[0].legend(fontsize=8)
optimum.plot(ax=axes[1])
axes[1].set(title="With and without the optimum absorber")
fig.tight_layout()

plt.show()
