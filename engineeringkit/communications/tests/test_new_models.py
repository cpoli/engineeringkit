"""Analytical checks for Carson's rule, Hamming(7,4), OFDM and array-valued inputs."""

import numpy as np
import pytest

from engineeringkit import communications as c


def test_carson_bandwidth_broadcast_fm():
    assert c.carson_bandwidth(75e3, 15e3) == 180e3
    with pytest.raises(ValueError):
        c.carson_bandwidth(-1, 15e3)


def test_carson_rule_captures_98_percent_of_fm_power():
    from scipy.special import jv

    beta = 5.0
    n_inside = int(beta + 1)  # sidebands within 2(delta f + W)
    power = jv(0, beta) ** 2 + 2 * sum(jv(n, beta) ** 2 for n in range(1, n_inside + 1))
    assert power > 0.98


def test_hamming74_all_codewords_have_minimum_distance_three():
    data = np.array([[(k >> i) & 1 for i in range(4)] for k in range(16)])
    words = [c.hamming74_encode(d) for d in data]
    distances = [c.hamming_distance(a, b) for i, a in enumerate(words) for b in words[i + 1 :]]
    assert min(distances) == 3


def test_hamming74_corrects_every_single_error():
    rng = np.random.default_rng(0)
    data = rng.integers(0, 2, 40)
    code = c.hamming74_encode(data)
    for position in range(7):
        received = code.copy()
        received[position::7] ^= 1
        result = c.hamming74_decode(received)
        assert np.array_equal(result.data, data)
        assert np.all(result.syndromes == position + 1)
        assert np.array_equal(result.corrected, code)
        assert result.summary() == {"words": 10, "corrected_words": 10}
    assert c.hamming74_decode(code).summary()["corrected_words"] == 0


@pytest.mark.parametrize("bad", [[1, 0, 1], [1, 0, 2, 1], [], [[1, 0, 1, 1]]])
def test_hamming74_rejects_bad_input(bad):
    with pytest.raises(ValueError):
        c.hamming74_encode(bad)
    with pytest.raises(ValueError):
        c.hamming74_decode(bad)


def test_ofdm_cyclic_prefix_diagonalises_a_multipath_channel():
    rng = np.random.default_rng(1)
    n, cp = 64, 8
    symbols = c.qpsk(rng.integers(0, 2, 2 * n * 5))
    channel = np.array([1.0, 0.5 - 0.2j, 0.25j])
    tx = c.ofdm_modulate(symbols, n, cp)
    rx = np.convolve(tx, channel)[: tx.size]
    response = np.fft.fft(channel, n)
    equalised = c.ofdm_demodulate(rx, n, cp).reshape(-1, n) / response
    assert np.allclose(equalised.reshape(-1), symbols)
    assert np.mean(np.abs(tx[cp : cp + n]) ** 2) == pytest.approx(np.mean(np.abs(symbols[:n]) ** 2))


@pytest.mark.parametrize(
    "args",
    [([1, 1, 1], 4, 0), ([1, 1, 1, 1], 1, 0), ([1, 1, 1, 1], 4, 4), ([np.nan, 1, 1, 1], 4, 0)],
)
def test_ofdm_rejects_bad_input(args):
    with pytest.raises(ValueError):
        c.ofdm_modulate(*args)
    with pytest.raises(ValueError):
        c.ofdm_demodulate(*args)


def test_array_inputs_and_new_validation():
    assert c.entropy_binary(np.array([0.0, 0.5, 1.0])) == pytest.approx([0, 1, 0])
    assert c.noise_figure(np.array([1.0, 10.0])) == pytest.approx([0, 10])
    for fn, args in [
        (c.from_db, (np.nan,)),
        (c.from_dbm, (np.inf,)),
        (c.q_function, (np.nan,)),
        (c.link_budget_db, (30, 0, 0, np.nan)),
        (c.am_signal, ([1], [0], -1)),
        (c.dsb_sc, ([np.nan], [0], 1)),
    ]:
        with pytest.raises(ValueError):
            fn(*args)


def test_qpsk_halves_bandwidth_at_equal_bit_error_rate():
    bits = np.array([0, 1, 1, 0, 1, 1, 0, 0])
    assert c.qpsk(bits).size == c.bpsk(bits).size // 2
    assert np.allclose(np.abs(c.qpsk(bits)) ** 2, 1)
