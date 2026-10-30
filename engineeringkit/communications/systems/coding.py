r"""Channel codes: Hamming(7,4), cyclic redundancy checks and convolutional codes.

Hamming placed three parity bits at positions 1, 2 and 4 of a seven-bit
word so that the parity-check failures, read as a binary number, spell out
the position of any single flipped bit. Peterson and Brown's cyclic codes
append the remainder of a polynomial division so that a receiver detects
errors by dividing again. Elias's convolutional codes pass the data through
a shift register, and Viterbi's algorithm decodes them by keeping one
survivor path per register state.
"""

from __future__ import annotations

import heapq

import numpy as np

from engineeringkit._validation import integer as _integer
from engineeringkit._validation import vector as _vector
from engineeringkit.communications.core.base import HammingDecodeResult, ViterbiDecodeResult

__all__ = [
    "hamming_distance",
    "hamming74_encode",
    "hamming74_decode",
    "CRC16_CCITT",
    "crc_remainder",
    "crc_encode",
    "crc_check",
    "convolutional_encode",
    "convolutional_free_distance",
    "viterbi_decode",
]

CRC16_CCITT = 0x11021
"""Generator polynomial :math:`x^{16} + x^{12} + x^5 + 1` of CRC-16-CCITT (ITU-T V.41), leading term included."""

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


def _generator_degree(generator):
    _integer("generator", generator, 2)
    return int(generator).bit_length() - 1


def _mod2_remainder(bits, generator, augment):
    degree = _generator_degree(generator)
    register = 0
    for bit in np.concatenate((bits, np.zeros(degree if augment else 0, dtype=np.int64))):
        register = (register << 1) | int(bit)
        if register >> degree:
            register ^= generator
    return np.array([(register >> k) & 1 for k in range(degree - 1, -1, -1)], dtype=np.int64)


def crc_remainder(bits, generator):
    r"""Cyclic-redundancy-check bits: the remainder of :math:`x^r m(x)` divided by :math:`g(x)` over GF(2).

    Parameters
    ----------
    bits : array_like
        Message bits, most significant (highest power of :math:`x`) first.
    generator : int
        Generator polynomial :math:`g(x)` of degree :math:`r \ge 1` written
        as an integer with its leading term, such as ``0b1011`` for
        :math:`x^3 + x + 1` or :data:`CRC16_CCITT`.

    Returns
    -------
    ndarray of int
        The :math:`r` check bits, most significant first.

    Notes
    -----
    Plain polynomial division: register preset to zero, no bit reflection
    and no final XOR (the "XMODEM" convention for CRC-16-CCITT). See
    :doc:`/guides/communications`.

    References
    ----------
    W. W. Peterson and D. T. Brown, "Cyclic codes for error detection,"
    *Proc. IRE* 49 (1961) 228-235.

    Examples
    --------
    >>> crc_remainder([1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0], 0b1011).tolist()
    [1, 0, 0]
    """
    return _mod2_remainder(_bits(bits, "bits", 1), generator, augment=True)


def crc_encode(bits, generator):
    """Append the CRC check bits to a message, giving a codeword divisible by the generator.

    Parameters
    ----------
    bits : array_like
        Message bits, most significant first.
    generator : int
        Generator polynomial with its leading term; see :func:`crc_remainder`.

    Returns
    -------
    ndarray of int
        The message followed by its check bits.

    Notes
    -----
    See :doc:`/guides/communications`.

    Examples
    --------
    >>> crc_encode([1, 0, 1, 1], 0b111).tolist()
    [1, 0, 1, 1, 0, 1]
    """
    message = _bits(bits, "bits", 1)
    return np.concatenate((message, _mod2_remainder(message, generator, augment=True)))


def crc_check(codeword, generator):
    """True when a received codeword leaves no remainder, i.e. no error was detected.

    Parameters
    ----------
    codeword : array_like
        Received bits, message then check bits.
    generator : int
        Generator polynomial with its leading term; see :func:`crc_remainder`.

    Returns
    -------
    bool
        ``True`` if :math:`g(x)` divides the received polynomial. A degree-
        :math:`r` generator with a nonzero constant term detects every burst
        of :math:`r` or fewer errors.

    Notes
    -----
    See :doc:`/guides/communications`.

    Examples
    --------
    >>> crc_check([1, 0, 1, 1, 0, 1], 0b111), crc_check([1, 1, 1, 1, 0, 1], 0b111)
    (True, False)
    """
    return not np.any(_mod2_remainder(_bits(codeword, "codeword", 1), generator, augment=False))


def _trellis(generators):
    g = tuple(generators) if np.ndim(generators) == 1 else ()
    if len(g) < 2:
        raise ValueError("generators must list at least two generator polynomials")
    for k, value in enumerate(g):
        _integer(f"generators[{k}]", value, 1)
    K = max(int(value).bit_length() for value in g)
    if K < 2:
        raise ValueError("the constraint length must be at least 2")
    states = 1 << (K - 1)
    register = (np.arange(2)[None, :] << (K - 1)) | np.arange(states)[:, None]
    outputs = np.stack([np.vectorize(lambda r, gj=int(gj): bin(r & gj).count("1") & 1)(register) for gj in g], axis=-1)
    return K, len(g), outputs, register >> 1


def convolutional_encode(bits, generators=(0o7, 0o5)):
    r"""Rate-:math:`1/n` feedforward convolutional encoder, terminated with :math:`K-1` zero tail bits.

    Each generator is an integer whose binary digits are the shift-register
    taps, most significant for the current input bit, so the octal pair
    ``(0o7, 0o5)`` is the classic constraint-length-3 code with
    :math:`g_0 = 1 + D + D^2` and :math:`g_1 = 1 + D^2`.

    Parameters
    ----------
    bits : array_like
        Data bits.
    generators : sequence of int, optional
        At least two generator polynomials; the longest sets the constraint
        length :math:`K`.

    Returns
    -------
    ndarray of int
        Code bits, :math:`n` per input bit (outputs of the generators in
        turn), for the data and the tail: :math:`n(L + K - 1)` in all.

    Notes
    -----
    See :doc:`/guides/communications`.

    References
    ----------
    P. Elias, "Coding for noisy channels," *IRE Convention Record* 3, part 4
    (1955) 37-46.

    Examples
    --------
    >>> convolutional_encode([1, 0, 1, 1]).tolist()
    [1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1]
    """
    K, _, _, _ = _trellis(generators)
    u = np.concatenate((_bits(bits, "bits", 1), np.zeros(K - 1, dtype=np.int64)))
    taps = [[(int(g) >> (K - 1 - i)) & 1 for i in range(K)] for g in generators]
    return np.stack([np.convolve(u, t)[: u.size] % 2 for t in taps], axis=1).reshape(-1)


def convolutional_free_distance(generators=(0o7, 0o5)):
    """Free distance of a convolutional code: the least weight of a path that leaves and rejoins the zero state.

    Parameters
    ----------
    generators : sequence of int, optional
        Generator polynomials as in :func:`convolutional_encode`.

    Returns
    -------
    int
        Free distance :math:`d_{free}`. A maximum-likelihood decoder corrects
        any pattern of fewer than :math:`d_{free}/2` errors spread across
        one such span.

    Notes
    -----
    Found by Dijkstra's shortest-path search over the trellis with Hamming
    weights as edge lengths. See :doc:`/guides/communications`.

    Examples
    --------
    >>> convolutional_free_distance((0o7, 0o5)), convolutional_free_distance((0o171, 0o133))
    (5, 10)
    """
    _, _, outputs, next_state = _trellis(generators)
    weights = outputs.sum(axis=-1)
    heap = [(int(weights[0, 1]), int(next_state[0, 1]))]
    settled: set[int] = set()
    while True:
        weight, state = heapq.heappop(heap)
        if state == 0:
            return weight
        if state in settled:
            continue
        settled.add(state)
        for u in (0, 1):
            heapq.heappush(heap, (weight + int(weights[state, u]), int(next_state[state, u])))


def viterbi_decode(received, generators=(0o7, 0o5), soft=False) -> ViterbiDecodeResult:
    r"""Maximum-likelihood decoding of a terminated convolutional code by the Viterbi algorithm.

    At each step every trellis state keeps only the best of the paths
    entering it (its *survivor*), so the work grows linearly with the
    message length instead of exponentially.

    Parameters
    ----------
    received : array_like
        Received code sequence from :func:`convolutional_encode`: hard bits
        (``soft=False``) or real soft values with bit 1 sent as ``+1`` and
        bit 0 as ``-1``, as in :func:`bpsk` (``soft=True``). The length must
        be a multiple of the number of generators and cover the tail.
    generators : sequence of int, optional
        Generator polynomials used by the encoder.
    soft : bool, optional
        Use the squared Euclidean distance (soft decisions, optimal in white
        Gaussian noise) instead of the Hamming distance.

    Returns
    -------
    ViterbiDecodeResult
        Decoded data bits and the surviving path's metric.

    Notes
    -----
    The encoder must have been terminated in the zero state, as
    :func:`convolutional_encode` does. See :doc:`/guides/communications`.

    References
    ----------
    A. J. Viterbi, "Error bounds for convolutional codes and an
    asymptotically optimum decoding algorithm," *IEEE Trans. Inf. Theory*
    13 (1967) 260-269.

    Examples
    --------
    >>> code = convolutional_encode([1, 0, 1, 1])
    >>> code[[2, 7]] ^= 1
    >>> viterbi_decode(code).data.tolist()
    [1, 0, 1, 1]
    """
    K, n, outputs, _ = _trellis(generators)
    if soft:
        r = _vector("received", received)
        if r.size % n:
            raise ValueError("received must hold a whole number of code symbols")
        reference = 2.0 * outputs - 1
    else:
        r = _bits(received, "received", n).astype(float)
        reference = outputs.astype(float)
    steps = r.size // n
    if steps < K:
        raise ValueError("received is shorter than one data bit plus the tail")
    low = (1 << (K - 2)) - 1
    destination = np.arange(1 << (K - 1))
    u = destination >> (K - 2)
    parents = np.stack((((destination & low) << 1), ((destination & low) << 1) | 1))
    metrics = np.full(destination.size, np.inf)
    metrics[0] = 0.0
    survivors = np.zeros((steps, destination.size), dtype=np.int64)
    for t, symbol in enumerate(r.reshape(steps, n)):
        distance = (np.sum((reference - symbol) ** 2, axis=-1) if soft else np.sum(reference != symbol, axis=-1)).astype(float)
        candidates = metrics[parents] + distance[parents, u]
        choice = np.argmin(candidates, axis=0)
        survivors[t] = parents[choice, destination]
        metrics = candidates[choice, destination]
        if t >= steps - (K - 1):
            metrics[u == 1] = np.inf
    state, decided = 0, np.zeros(steps, dtype=np.int64)
    for t in range(steps - 1, -1, -1):
        decided[t] = state >> (K - 2)
        state = survivors[t, state]
    return ViterbiDecodeResult(decided[: steps - (K - 1)], float(metrics[0]))
