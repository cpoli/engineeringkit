r"""Error-correcting codes: Hamming distance and the Hamming(7,4) code.

Hamming placed three parity bits at positions 1, 2 and 4 of a seven-bit
word so that the parity-check failures, read as a binary number, spell out
the position of any single flipped bit.
"""

from __future__ import annotations

import numpy as np

from engineeringkit.communications.core.base import HammingDecodeResult

__all__ = ["hamming_distance", "hamming74_encode", "hamming74_decode"]

# Columns of the parity-check matrix are the binary expansions of 1..7.
_H = np.array([[(position >> bit) & 1 for position in range(1, 8)] for bit in range(3)])
_DATA_POSITIONS = np.array([2, 4, 5, 6])
_PARITY_POSITIONS = np.array([0, 1, 3])


def _bits(values, name, block):
    b = np.asarray(values)
    if b.ndim != 1 or not b.size or b.size % block or not np.all((b == 0) | (b == 1)):
        raise ValueError(f"{name} must be a nonempty binary vector whose length is a multiple of {block}")
    return b.astype(np.int64)


def hamming_distance(a, b):
    """Number of positions at which two equal-length binary words differ.

    Parameters
    ----------
    a, b : array_like
        Binary words of the same shape.

    Returns
    -------
    int
        Hamming distance.

    Notes
    -----
    See :doc:`/guides/communications`.

    Examples
    --------
    >>> hamming_distance([1, 0, 1, 1], [1, 1, 1, 0])
    2
    """
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or not a.size or not np.all((a == 0) | (a == 1)) or not np.all((b == 0) | (b == 1)):
        raise ValueError("matching nonempty binary arrays are required")
    return int(np.count_nonzero(a != b))


def hamming74_encode(data):
    """Encode data bits with the Hamming(7,4) single-error-correcting code.

    Parameters
    ----------
    data : array_like
        Binary data; its length must be a multiple of four.

    Returns
    -------
    ndarray of int
        Codewords ``[p1, p2, d1, p3, d2, d3, d4]`` concatenated, 7/4 as long
        as the input. Any two codewords differ in at least three positions.

    Notes
    -----
    See :doc:`/guides/communications`.

    References
    ----------
    R. W. Hamming, "Error detecting and error correcting codes," *Bell Syst.
    Tech. J.* 29 (1950) 147-160.

    Examples
    --------
    >>> hamming74_encode([1, 0, 1, 1]).tolist()
    [0, 1, 1, 0, 0, 1, 1]
    """
    d = _bits(data, "data", 4).reshape(-1, 4)
    words = np.zeros((len(d), 7), dtype=np.int64)
    words[:, _DATA_POSITIONS] = d
    words[:, _PARITY_POSITIONS] = (words @ _H.T) % 2
    return words.reshape(-1)


def hamming74_decode(codewords) -> HammingDecodeResult:
    """Decode Hamming(7,4) codewords, correcting up to one flipped bit per word.

    Parameters
    ----------
    codewords : array_like
        Received binary words; length must be a multiple of seven.

    Returns
    -------
    HammingDecodeResult
        Decoded data bits, the syndrome of each word (the 1-based position
        of the corrected bit, or 0), and the corrected codewords.

    Notes
    -----
    Two or more errors in one word are silently mis-corrected; the code's
    minimum distance of three guarantees only single-error correction. See
    :doc:`/guides/communications`.

    Examples
    --------
    >>> received = hamming74_encode([1, 0, 1, 1])
    >>> received[4] ^= 1
    >>> result = hamming74_decode(received)
    >>> result.data.tolist(), result.syndromes.tolist()
    ([1, 0, 1, 1], [5])
    """
    words = _bits(codewords, "codewords", 7).reshape(-1, 7).copy()
    syndromes = ((words @ _H.T) % 2) @ np.array([1, 2, 4])
    rows = np.flatnonzero(syndromes)
    words[rows, syndromes[rows] - 1] ^= 1
    return HammingDecodeResult(words[:, _DATA_POSITIONS].reshape(-1), syndromes, words.reshape(-1))
