r"""
Hamming's single-error-correcting code (1950)
=============================================

Frustrated by weekend computer runs aborted by a single bad relay, Richard
Hamming (1950) designed codes that *correct* errors rather than merely
detect them. In the (7,4) code three parity bits sit at positions 1, 2 and 4.
When one bit flips, the failing parity checks, read as a binary number, give
its position (:func:`~engineeringkit.communications.hamming74_encode`,
:func:`~engineeringkit.communications.hamming74_decode`). Any two codewords
differ in at least three places (:func:`~engineeringkit.communications.hamming_distance`),
which is what makes single-error correction possible.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.communications import hamming74_decode, hamming74_encode, hamming_distance

# %%
# The syndrome points at the error
# --------------------------------

word = hamming74_encode([1, 0, 1, 1])
print("codeword:", word.tolist())
for position in range(7):
    corrupted = word.copy()
    corrupted[position] ^= 1
    result = hamming74_decode(corrupted)
    print(f"flip bit {position + 1}: syndrome {result.syndromes[0]} -> data {result.data.tolist()}")
    assert result.syndromes[0] == position + 1 and result.data.tolist() == [1, 0, 1, 1]

codebook = [hamming74_encode([(k >> i) & 1 for i in range(4)]) for k in range(16)]
distances = np.array([[hamming_distance(a, b) for b in codebook] for a in codebook])
print("minimum distance:", distances[~np.eye(16, dtype=bool)].min())

# %%
# Coded versus uncoded error rate on a bit-flipping channel
# ---------------------------------------------------------
#
# A block fails only if two or more of its seven bits flip:
# :math:`1 - (1-p)^7 - 7p(1-p)^6`, compared with :math:`1-(1-p)^4` uncoded.

rng = np.random.default_rng(1950)
p = np.logspace(-3, -1, 7)
simulated = []
for flip in p:
    data = rng.integers(0, 2, 4 * 50_000)
    received = hamming74_encode(data) ^ (rng.random(7 * 50_000) < flip)
    blocks_wrong = np.any((hamming74_decode(received).data != data).reshape(-1, 4), axis=1)
    simulated.append(blocks_wrong.mean())
coded = 1 - (1 - p) ** 7 - 7 * p * (1 - p) ** 6
assert np.allclose(simulated[3:], coded[3:], rtol=0.15)
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
axes[0].imshow(distances, cmap="viridis")
axes[0].set(title="Distances between the 16 codewords", xlabel="codeword", ylabel="codeword")
axes[1].loglog(p, 1 - (1 - p) ** 4, label="uncoded 4-bit block")
axes[1].loglog(p, coded, label="Hamming(7,4), theory")
axes[1].loglog(p, np.maximum(simulated, 1e-6), "o", label="Hamming(7,4), simulation")
axes[1].set(xlabel="bit-flip probability", ylabel="block error rate")
axes[1].legend()
fig.tight_layout()

plt.show()
