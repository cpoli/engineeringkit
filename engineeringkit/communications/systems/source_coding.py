r"""Source coding: the entropy of a discrete source and Huffman's minimum-redundancy code.

Shannon's source-coding theorem says no uniquely decodable code can average
fewer than :math:`H` bits per symbol. Huffman (1952) built the prefix-free
code that comes closest, by repeatedly merging the two least likely symbols.
"""

from __future__ import annotations

import heapq

import numpy as np
from scipy import special

from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import vector as _vector
from engineeringkit.communications.core.base import HuffmanCode

__all__ = ["entropy", "huffman_code"]


def _probabilities(probabilities):
    p = _vector("probabilities", probabilities, min_length=2)
    _nonnegative(probabilities=p)
    if not np.isclose(p.sum(), 1.0, rtol=0, atol=1e-9):
        raise ValueError("probabilities must sum to one")
    return p


def entropy(probabilities):
    r"""Entropy of a discrete memoryless source, :math:`H = -\sum_i p_i\log_2 p_i`.

    Parameters
    ----------
    probabilities : array_like
        Symbol probabilities [-], at least two, nonnegative and summing to
        one.

    Returns
    -------
    float
        Entropy [bit per symbol], between zero and :math:`\log_2 M`.

    Notes
    -----
    :func:`entropy_binary` is the two-symbol case. Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/communications`.

    References
    ----------
    C. E. Shannon, "A mathematical theory of communication," *Bell Syst.
    Tech. J.* 27 (1948) 379-423, 623-656.

    Examples
    --------
    >>> entropy([0.5, 0.25, 0.125, 0.125])
    1.75
    """
    p = _probabilities(probabilities)
    return float(-np.sum(special.xlogy(p, p)) / np.log(2))


def huffman_code(probabilities) -> HuffmanCode:
    r"""Huffman's optimal binary prefix code for a discrete memoryless source.

    The two least probable entries are merged into one, with a ``0`` and a
    ``1`` prepended to their codewords, until a single entry remains. Ties
    are broken by order of creation, so the result is deterministic.

    Parameters
    ----------
    probabilities : array_like
        Symbol probabilities [-], at least two, nonnegative and summing to
        one.

    Returns
    -------
    HuffmanCode
        Codeword table with :attr:`~HuffmanCode.average_length`,
        :attr:`~HuffmanCode.entropy`, ``encode`` and ``decode``. The average
        length satisfies :math:`H \le \bar L < H + 1`.

    Notes
    -----
    See :doc:`/guides/communications`.

    References
    ----------
    D. A. Huffman, "A method for the construction of minimum-redundancy
    codes," *Proc. IRE* 40 (1952) 1098-1101.

    Examples
    --------
    >>> code = huffman_code([0.4, 0.2, 0.2, 0.1, 0.1])
    >>> round(code.average_length, 6)
    2.2
    """
    p = _probabilities(probabilities)
    words = [""] * p.size
    heap = [(float(q), k, [k]) for k, q in enumerate(p)]
    heapq.heapify(heap)
    order = p.size
    while len(heap) > 1:
        q0, _, group0 = heapq.heappop(heap)
        q1, _, group1 = heapq.heappop(heap)
        for k in group0:
            words[k] = "0" + words[k]
        for k in group1:
            words[k] = "1" + words[k]
        heapq.heappush(heap, (q0 + q1, order, group0 + group1))
        order += 1
    return HuffmanCode(p, tuple(words))
