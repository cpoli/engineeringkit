"""Result containers for engineeringkit.communications.

Most communication models return plain arrays or floats. Channel decoding
returns more than the data alone (which words were corrected, and where, or
the metric of the surviving trellis path), and a Huffman code is a table that
encodes and decodes, so these get frozen dataclasses with the usual
``summary()``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import special

__all__ = ["HammingDecodeResult", "HuffmanCode", "ViterbiDecodeResult"]


@dataclass(frozen=True)
class HammingDecodeResult:
    """Output of :func:`engineeringkit.communications.hamming74_decode`."""

    data: np.ndarray
    """ndarray of int: Decoded data bits, four per codeword."""

    syndromes: np.ndarray
    """ndarray of int: Per-word syndrome, the 1-based position of the corrected bit (0 if none)."""

    corrected: np.ndarray
    """ndarray of int: Codewords after correction."""

    def summary(self) -> dict:
        """Number of words decoded and number in which a bit was corrected."""
        return {"words": int(self.syndromes.size), "corrected_words": int(np.count_nonzero(self.syndromes))}


@dataclass(frozen=True)
class ViterbiDecodeResult:
    """Output of :func:`engineeringkit.communications.viterbi_decode`."""

    data: np.ndarray
    """ndarray of int: Decoded data bits, without the zero tail."""

    path_metric: float
    """float: Metric of the surviving path: the Hamming distance between the
    received bits and the re-encoded decision (hard decisions), or the squared
    Euclidean distance to its ``±1`` symbols (soft decisions)."""

    def summary(self) -> dict:
        """Number of decoded bits and the surviving path's metric."""
        return {"bits": int(self.data.size), "path_metric": float(self.path_metric)}


@dataclass(frozen=True)
class HuffmanCode:
    """Output of :func:`engineeringkit.communications.huffman_code`: a prefix-free code table."""

    probabilities: np.ndarray
    """ndarray of float: Symbol probabilities, in the order given."""

    codewords: tuple
    """tuple of str: Binary codeword of each symbol, such as ``"010"``."""

    @property
    def lengths(self) -> np.ndarray:
        """Codeword length of each symbol [bit]."""
        return np.array([len(word) for word in self.codewords])

    @property
    def average_length(self) -> float:
        r"""Expected codeword length :math:`\bar L = \sum p_i \ell_i` [bit per symbol]."""
        return float(self.probabilities @ self.lengths)

    @property
    def entropy(self) -> float:
        r"""Source entropy :math:`H = -\sum p_i \log_2 p_i` [bit per symbol], the lower bound on :attr:`average_length`."""
        return float(-np.sum(special.xlogy(self.probabilities, self.probabilities)) / np.log(2))

    def encode(self, symbols) -> np.ndarray:
        """Concatenate the codewords of a sequence of symbol indices.

        Parameters
        ----------
        symbols : array_like of int
            Indices into :attr:`probabilities`.

        Returns
        -------
        ndarray of int
            The encoded bit stream.
        """
        s = np.asarray(symbols)
        if s.ndim != 1 or not s.size or not np.issubdtype(s.dtype, np.integer) or np.any((s < 0) | (s >= len(self.codewords))):
            raise ValueError("symbols must be a nonempty vector of valid symbol indices")
        return np.array([int(bit) for index in s for bit in self.codewords[index]], dtype=np.int64)

    def decode(self, bits) -> np.ndarray:
        """Split a bit stream back into symbol indices; the prefix property makes this unambiguous.

        Parameters
        ----------
        bits : array_like
            Binary stream made of whole codewords.

        Returns
        -------
        ndarray of int
            The decoded symbol indices.
        """
        b = np.asarray(bits)
        if b.ndim != 1 or not b.size or not np.all((b == 0) | (b == 1)):
            raise ValueError("bits must be a nonempty binary vector")
        table = {word: index for index, word in enumerate(self.codewords)}
        symbols, word = [], ""
        for bit in b.astype(int):
            word += str(bit)
            if word in table:
                symbols.append(table[word])
                word = ""
        if word:
            raise ValueError("bits end in the middle of a codeword")
        return np.array(symbols, dtype=np.int64)

    def summary(self) -> dict:
        """Entropy, average length and efficiency ``H / L``."""
        return {
            "symbols": len(self.codewords),
            "entropy": self.entropy,
            "average_length": self.average_length,
            "efficiency": self.entropy / self.average_length,
        }

    def plot(self, ax=None):
        r"""Bar chart of codeword lengths next to the ideal lengths :math:`-\log_2 p_i`; return ``(fig, ax)``.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure is created if omitted.

        Returns
        -------
        tuple
            ``(fig, ax)``.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        index = np.arange(len(self.codewords))
        positive = self.probabilities > 0
        ax.bar(index - 0.2, self.lengths, 0.4, label="Huffman length")
        ax.bar(index[positive] + 0.2, -np.log2(self.probabilities[positive]), 0.4, label=r"$-\log_2 p$")
        ax.set_xticks(index, list(self.codewords), rotation=45)
        ax.set(xlabel="codeword", ylabel="length [bit]")
        ax.legend()
        return ax.figure, ax
