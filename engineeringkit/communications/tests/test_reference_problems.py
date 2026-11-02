"""Worked examples from S. Haykin, *Digital Communications*, Wiley, 1988 (ISBN 978-0-471-62947-4).

Example numbers, pages and printed answers are as reproduced in the FOSSEE
Python Textbook Companion for this book (Wiley India reprint pages),
https://tbc-python.fossee.in/book-details/977/. Each answer is checked to
the precision printed there.
"""

import numpy as np

from engineeringkit import communications as c


def test_haykin_example_2_2_second_order_extension():
    """Haykin (1988), Example 2.2, p. 19: H = 1.5 bits, and 3 bits for the second-order extension."""
    p = np.array([0.25, 0.25, 0.5])
    assert abs(c.entropy(p) - 1.5) < 1e-10
    assert abs(c.entropy(np.outer(p, p).ravel()) - 3.0) < 1e-10


def test_haykin_example_2_3_huffman_code():
    """Haykin (1988), Example 2.3, p. 20: L = 2.2 bits, H = 2.12 bits, length variance 0.16."""
    code = c.huffman_code([0.4, 0.2, 0.2, 0.1, 0.1])
    assert code.lengths.tolist() == [2, 2, 2, 3, 3]
    assert abs(code.average_length - 2.2) < 1e-10
    assert round(code.entropy, 2) == 2.12
    variance = code.probabilities @ (code.lengths - code.average_length) ** 2
    assert abs(variance - 0.16) < 1e-10


def test_haykin_example_5_3_pcm_channel_bandwidth_and_snr():
    """Haykin (1988), Example 5.3, p. 192: 8-bit PCM of a 4 kHz message needs B = 32 kHz and gives 40.8 dB.

    The book's (SNR)o = 6n - 7.2 dB loads the quantizer to four times the
    rms signal; pcm_sqnr_db is for a full-scale sine (peak = sqrt(2) rms),
    so the loading costs 10 log10(16/2) dB. The 0.09 dB left over is the
    book's rounding of 6.02 dB per bit to 6.
    """
    n, W = 8, 4000.0
    bit_rate = n * c.nyquist_symbol_rate(W)
    assert abs(bit_rate / c.nyquist_symbol_rate(1.0) - 32000.0) < 1e-9
    snr = c.pcm_sqnr_db(n) - c.db(16 / 2)
    assert abs(snr - 40.8) < 0.1


def test_haykin_example_8_3_cyclic_hamming_code():
    """Haykin (1988), Example 8.3, p. 389: g(X) = 1 + X + X^3 and message 1001 give parity 011, codeword 0111001.

    Haykin lists bits lowest power first; the package writes them highest
    power first, so each printed word appears reversed.
    """
    g = 0b1011
    assert c.crc_remainder([1, 0, 0, 1], g).tolist()[::-1] == [0, 1, 1]
    assert c.crc_encode([1, 0, 0, 1], g).tolist()[::-1] == [0, 1, 1, 1, 0, 0, 1]


def test_haykin_example_8_5_cyclic_code_syndrome():
    """Haykin (1988), Example 8.5, p. 396: codeword 0111001 has syndrome 000; with its middle bit in error, 110.

    The syndrome is the remainder of the received word on g(X), which for a
    single error at X^3 is the remainder of X^3, i.e. crc_remainder([1]).
    """
    g = 0b1011
    assert c.crc_check([1, 0, 0, 1, 1, 1, 0], g)
    assert not c.crc_check([1, 0, 0, 0, 1, 1, 0], g)
    assert c.crc_remainder([1], g).tolist()[::-1] == [1, 1, 0]


def test_haykin_example_8_7_convolutional_encoder():
    """Haykin (1988), Example 8.7, p. 401: g1 = 1 + D + D^2, g2 = 1 + D^2 encode 10011 as 11 10 11 11 01 01 11."""
    encoded = c.convolutional_encode([1, 0, 0, 1, 1], (0o7, 0o5))
    assert encoded.tolist() == [1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1]
