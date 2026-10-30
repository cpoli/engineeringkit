r"""
Huffman's minimum-redundancy code (1952)
========================================

In Robert Fano's information theory class at MIT, David Huffman was offered
a term paper instead of a final exam: find the most efficient binary code.
Shannon and Fano's top-down splitting was good but not always optimal.
Huffman worked bottom-up. He merged the two least likely symbols into one,
and repeated until a single node remained. The resulting prefix code
(:func:`~engineeringkit.communications.huffman_code`) has the least average
length of any symbol-by-symbol code, within one bit of the entropy
(:func:`~engineeringkit.communications.entropy`), and reaches it exactly
when every probability is a power of one half.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.communications import entropy, huffman_code

# %%
# Coding the letters of English
# -----------------------------
#
# Approximate single-letter frequencies of English text (space excluded).

letters = "ETAOINSHRDLCUMWFGYPBVKJXQZ"
freq = np.array(
    [12.7, 9.1, 8.2, 7.5, 7.0, 6.7, 6.3, 6.1, 6.0, 4.3, 4.0, 2.8, 2.8, 2.4, 2.4, 2.2, 2.0, 2.0, 1.9, 1.5, 1.0, 0.8, 0.15, 0.15, 0.1, 0.07]
)
p = freq / freq.sum()
code = huffman_code(p)
for letter, word in list(zip(letters, code.codewords))[:6]:
    print(f"{letter}: {word}")
summary = code.summary()
print(f"entropy {summary['entropy']:.3f} bit, Huffman {summary['average_length']:.3f} bit, fixed-length 5 bit")
assert code.entropy <= code.average_length < code.entropy + 1
assert np.isclose(sum(2.0**-code.lengths), 1)

message = "HUFFMAN"
bits = code.encode([letters.index(ch) for ch in message])
assert "".join(letters[k] for k in code.decode(bits)) == message
print(f"'{message}' takes {bits.size} bits instead of {5 * len(message)}")

# %%
# A dyadic source is coded exactly at its entropy
# -----------------------------------------------

dyadic = huffman_code([0.5, 0.25, 0.125, 0.0625, 0.0625])
print("dyadic codewords:", dyadic.codewords)
assert np.isclose(dyadic.average_length, entropy([0.5, 0.25, 0.125, 0.0625, 0.0625]))

fig, axes = plt.subplots(1, 2, figsize=(12, 3.8), gridspec_kw={"width_ratios": [3, 1]})
code.plot(ax=axes[0])
axes[0].set_xticks(np.arange(len(letters)), list(letters), rotation=0)
axes[0].set(xlabel="letter", title="Huffman lengths follow $-\\log_2 p$")
axes[1].bar(["entropy", "Huffman", "fixed"], [code.entropy, code.average_length, 5.0], color=["C2", "C0", "C7"])
axes[1].set(ylabel="bit per letter")
fig.tight_layout()

plt.show()
