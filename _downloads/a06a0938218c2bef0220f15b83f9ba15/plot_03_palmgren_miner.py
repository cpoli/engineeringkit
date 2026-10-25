r"""
Palmgren and Miner's cumulative damage (1945)
=============================================

Real loads vary in amplitude. Arvid Palmgren (1924), estimating the lives of
ball bearings, assumed that each cycle at amplitude :math:`\sigma_i` uses up
a fraction :math:`1/N_i` of the life, and Milton Miner (1945) confirmed the
rule on riveted aluminium joints for aircraft: failure comes when
:math:`D = \sum n_i/N_i` reaches 1
(:func:`~engineeringkit.materials.miner_damage`). The example applies it to a
block load spectrum and reduces the spectrum to one equivalent amplitude.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import basquin_cycles, miner_damage

# %%
# A load spectrum
# ---------------
#
# One block of service is 10,000 cycles: a few large ones and many small.
# The S-N curve gives the life at each amplitude; the block's damage is
# the sum of the fractions.

sigma_f, b = 900e6, -0.09
amplitude = np.array([450e6, 380e6, 320e6, 260e6])
per_block = np.array([20.0, 300.0, 2680.0, 7000.0])
lives = basquin_cycles(amplitude, sigma_f, b)
damage = miner_damage(per_block, lives)
print(f"damage per block {damage:.4f}; life {1 / damage:.1f} blocks = {per_block.sum() / damage:.3g} cycles")
for s, n, N in zip(amplitude, per_block, lives):
    print(f"  {s / 1e6:.0f} MPa: {n:6.0f} cycles of {N:9.3g} -> {100 * n / N / damage:4.1f}% of the damage")

# %%
# The equivalent constant amplitude
# ---------------------------------
#
# With Basquin's law the damage of a cycle grows as
# :math:`\sigma^{k}`, :math:`k = -1/b`. The whole spectrum then does the
# same damage as :math:`\sum n_i` cycles at
# :math:`\sigma_{eq} = (\sum n_i\sigma_i^k / \sum n_i)^{1/k}`. Miner's sum
# ignores the order: high-low and low-high sequences give the same answer,
# although tests show high-low is usually worse.

k = -1 / b
sigma_eq = (np.sum(per_block * amplitude**k) / per_block.sum()) ** (1 / k)
equivalent = miner_damage([per_block.sum()], [basquin_cycles(sigma_eq, sigma_f, b)])
print(f"equivalent amplitude {sigma_eq / 1e6:.1f} MPa")
assert np.isclose(equivalent, damage)
assert np.isclose(miner_damage(per_block[::-1], lives[::-1]), damage)

fig, ax = plt.subplots(figsize=(6, 3.8))
ax.bar([f"{s / 1e6:.0f}" for s in amplitude], per_block / lives / damage)
ax.set(xlabel="stress amplitude [MPa]", ylabel="share of damage", title=f"Damage per block, equivalent amplitude {sigma_eq / 1e6:.0f} MPa")
fig.tight_layout()

plt.show()
