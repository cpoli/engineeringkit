r"""
Elias's convolutional codes (1955)
==================================

Peter Elias's "Coding for noisy channels" showed that a code need not cut
the data into blocks. A shift register of length :math:`K` slides along the
bit stream, and each input bit produces :math:`n` output bits, sums modulo
two of the register taps (:func:`~engineeringkit.communications.convolutional_encode`).
Elias proved that such codes are at least as good as block codes. Their
strength is the *free distance*, the least weight of any path that leaves
the all-zero state and returns to it
(:func:`~engineeringkit.communications.convolutional_free_distance`).
For the standard rate-1/2 codes it grows with the constraint length.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.communications import convolutional_encode, convolutional_free_distance

# %%
# The (7, 5) encoder
# ------------------
#
# Generators :math:`1 + D + D^2` and :math:`1 + D^2`, octal 7 and 5.

data = [1, 0, 1, 1, 0, 0, 1]
code = convolutional_encode(data)
print("data:", data)
print("code:", " ".join("".join(map(str, pair)) for pair in code.reshape(-1, 2)))
impulse = convolutional_encode([1])
print("impulse response:", impulse.tolist(), "weight", impulse.sum())
assert impulse.sum() == convolutional_free_distance((0o7, 0o5)) == 5

# The code is linear: the encoding of a sum is the sum of the encodings.
rng = np.random.default_rng(1955)
a, b = rng.integers(0, 2, 50), rng.integers(0, 2, 50)
assert np.array_equal(convolutional_encode(a ^ b), convolutional_encode(a) ^ convolutional_encode(b))

# %%
# Free distance grows with the constraint length
# ----------------------------------------------
#
# The best known rate-1/2 generators (Proakis, Table 8.3-1).

best = {3: (0o5, 0o7), 4: (0o15, 0o17), 5: (0o23, 0o35), 6: (0o53, 0o75), 7: (0o133, 0o171)}
known = {3: 5, 4: 6, 5: 7, 6: 8, 7: 10}
distances = {K: convolutional_free_distance(g) for K, g in best.items()}
for K, d in distances.items():
    print(f"K = {K}: generators {tuple(oct(g) for g in best[K])}, d_free = {d}")
assert distances == known

# Distances between codewords of 4-bit messages are never below d_free.
words = [convolutional_encode([(k >> i) & 1 for i in range(4)]) for k in range(16)]
pairwise = np.array([[np.count_nonzero(u != v) for v in words] for u in words])
assert pairwise[~np.eye(16, dtype=bool)].min() == 5

fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
axes[0].stem(np.arange(impulse.size), impulse)
axes[0].set(xlabel="code bit", title="Impulse response of the (7, 5) code: weight 5")
axes[1].plot(list(distances), list(distances.values()), "o-")
axes[1].set(xlabel="constraint length K", ylabel="free distance", title="Best rate-1/2 codes")
fig.tight_layout()

plt.show()
