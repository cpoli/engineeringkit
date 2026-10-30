r"""
Viterbi's decoding algorithm (1967)
===================================

Decoding a convolutional code by brute force means comparing the received
sequence with all :math:`2^L` possible codewords. Andrew Viterbi saw that
two paths entering the same encoder state at the same time have the same
future, so only the better of them, the *survivor*, need be kept. The work
then grows linearly with :math:`L`
(:func:`~engineeringkit.communications.viterbi_decode`, returning a
:class:`~engineeringkit.communications.ViterbiDecodeResult`). Viterbi
introduced the algorithm as a device for proving error bounds. Jim Omura
and David Forney showed in the early 1970s that it is maximum-likelihood
decoding, and it became the decoder of deep-space links, GSM telephones and
dial-up modems. Soft decisions, the received values before slicing, gain
about 2 dB over hard bits.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.communications import ber_bpsk_awgn, bpsk, convolutional_encode, from_db, viterbi_decode

# %%
# Viterbi equals brute-force maximum likelihood
# ---------------------------------------------

rng = np.random.default_rng(1967)
codebook = np.array([convolutional_encode([(k >> i) & 1 for i in range(10)]) for k in range(1024)])
for _ in range(20):
    received = rng.integers(0, 2, codebook.shape[1])
    best = np.min(np.count_nonzero(codebook != received, axis=1))
    assert viterbi_decode(received).path_metric == best
print("20 random words: the survivor metric always equals the exhaustive minimum over 1024 codewords")

# %%
# Coding gain of the K = 7 code in white Gaussian noise
# -----------------------------------------------------

generators = (0o171, 0o133)
ebn0_db = np.arange(0, 7.5, 1.0)
hard, soft = [], []
data = rng.integers(0, 2, 20_000)
code = convolutional_encode(data, generators)
for g in from_db(ebn0_db):
    sigma = np.sqrt(1 / (2 * 0.5 * g))  # unit symbol energy, rate 1/2
    y = bpsk(code) + sigma * rng.standard_normal(code.size)
    hard.append(np.mean(viterbi_decode((y > 0).astype(int), generators).data != data))
    soft.append(np.mean(viterbi_decode(y, generators, soft=True).data != data))
hard, soft = np.array(hard), np.array(soft)
for e, h, s in zip(ebn0_db, hard, soft):
    print(f"Eb/N0 {e:.0f} dB: uncoded {ber_bpsk_awgn(from_db(e)):.1e}, hard {h:.1e}, soft {s:.1e}")
at4 = list(ebn0_db).index(4.0)
assert soft[at4] < hard[at4] < ber_bpsk_awgn(from_db(4.0))

fig, ax = plt.subplots(figsize=(6.5, 4))
fine = np.linspace(0, 10, 101)
ax.semilogy(fine, ber_bpsk_awgn(from_db(fine)), label="uncoded BPSK")
ax.semilogy(ebn0_db[hard > 0], hard[hard > 0], "o-", label="K = 7, hard decisions")
ax.semilogy(ebn0_db[soft > 0], soft[soft > 0], "s-", label="K = 7, soft decisions")
ax.set(xlabel="$E_b/N_0$ [dB]", ylabel="bit-error rate", ylim=(1e-5, 0.5))
ax.grid(True, which="both", alpha=0.3)
ax.legend()
fig.tight_layout()

plt.show()
