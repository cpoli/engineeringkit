r"""
Shannon's entropy and the binary symmetric channel
==================================================

Claude Shannon's "A mathematical theory of communication" (1948) measured
information by entropy, :math:`H = -\sum p\log_2 p`. For a binary source
:func:`~engineeringkit.communications.entropy_binary` peaks at one bit when
both symbols are equally likely. The same paper defined channel capacity. A
binary symmetric channel that flips each bit with probability :math:`p` has
capacity :math:`1 - H(p)`
(:func:`~engineeringkit.communications.binary_symmetric_capacity`): even a
very noisy channel can carry information reliably at a reduced rate, and a
channel with :math:`p = 1/2` carries nothing.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.communications import binary_symmetric_capacity, entropy_binary

# %%
# Entropy as the limit of lossless compression
# --------------------------------------------
#
# The empirical entropy of a long biased bit stream matches :math:`H(p)`.

rng = np.random.default_rng(1948)
for p in [0.5, 0.1, 0.01]:
    bits = rng.random(200_000) < p
    blocks = bits[: bits.size // 8 * 8].reshape(-1, 8) @ (1 << np.arange(8))
    counts = np.bincount(blocks, minlength=256) / blocks.size
    empirical = -np.sum(counts[counts > 0] * np.log2(counts[counts > 0])) / 8
    print(f"p = {p:4}: H(p) = {entropy_binary(p):.4f} bit, empirical (8-bit blocks) = {empirical:.4f} bit")
    assert np.isclose(empirical, entropy_binary(p), atol=0.01)

# %%
# Capacity of a bit-flipping channel
# ----------------------------------
#
# Mutual information :math:`I = H(Y) - H(Y|X)` over input distributions
# peaks at the uniform input, where it equals :math:`1 - H(p)`.

p = np.linspace(0, 1, 201)
q = np.linspace(0, 1, 201)
crossover = 0.11
mutual = entropy_binary(q * (1 - crossover) + (1 - q) * crossover) - entropy_binary(crossover)
assert np.isclose(mutual.max(), binary_symmetric_capacity(crossover))
fig, axes = plt.subplots(1, 3, figsize=(11, 3.3))
axes[0].plot(p, entropy_binary(p))
axes[0].set(xlabel="P(1)", ylabel="H [bit]", title="Binary entropy")
axes[1].plot(q, mutual)
axes[1].set(xlabel="input P(1)", ylabel="I(X;Y) [bit]", title=f"BSC, p = {crossover}")
axes[2].plot(p, binary_symmetric_capacity(p))
axes[2].set(xlabel="crossover probability p", ylabel="capacity [bit/use]", title="C = 1 - H(p)")
fig.tight_layout()

plt.show()
