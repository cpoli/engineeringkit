import numpy as np
import pytest

from engineeringkit import communications as c


def test_fm_constant_message_nonuniform_grid_and_initial_phase():
    times = np.array([0, 0.001, 0.003, 0.01, 0.02])
    modulation = c.fm_signal(np.ones(5) * 2, times, 10, 3)
    np.testing.assert_allclose(modulation, np.cos(2 * np.pi * 16 * times), atol=1e-14)
    assert modulation[0] == 1


def test_linear_message_phase_integral():
    times = np.linspace(0, 0.1, 101)
    np.testing.assert_allclose(c.fm_signal(2 * times, times, 10, 2), np.cos(2 * np.pi * 10 * times + 4 * np.pi * times**2), atol=1e-14)


def test_unsigned_modulation_does_not_wrap():
    bits = np.array([0, 1, 1, 0], dtype=np.uint8)
    np.testing.assert_allclose(c.bpsk(bits), [-1, 1, 1, -1])
    symbols = c.qpsk(bits)
    np.testing.assert_allclose(abs(symbols), 1)
    assert symbols[0].real < 0 and symbols[0].imag > 0


def test_link_budget_agrees_with_linear_friis():
    power = c.friis_received_power(1, 2, 3, 100, 2e9)
    budget = c.link_budget_db(30, c.db(2), c.db(3), c.free_space_path_loss_db(100, 2e9))
    assert np.isclose(c.from_dbm(budget), power)
    assert np.isclose(c.free_space_path_loss(200, 2e9) / c.free_space_path_loss(100, 2e9), 4)


def test_noise_and_channel_limits():
    assert np.isclose(c.cascaded_noise_factor([10, 10], [2, 3]), 2.2)
    assert c.shannon_capacity(1e6, 0) == 0
    assert c.shannon_capacity(1e6, 1) == 1e6
    assert c.binary_symmetric_capacity(0.5) == 0
    assert c.binary_symmetric_capacity(0) == 1
    assert np.isclose(c.ber_bpsk_awgn(1), 0.07864960352514257)


@pytest.mark.parametrize("bits", [[2], [np.nan], []])
def test_nonbinary_modulation(bits):
    with pytest.raises(ValueError):
        c.bpsk(bits)


def test_fm_and_noise_shape_validation():
    with pytest.raises(ValueError):
        c.fm_signal([1, 1], [0, 0], 10, 2)
    with pytest.raises(ValueError):
        c.cascaded_noise_factor([10], [2, 3])
    with pytest.raises(ValueError):
        c.hamming_distance([0, 1], [0])


def test_pcm_levels_agree_with_the_signal_processing_quantizer():
    from engineeringkit.signal_processing import quantize

    x = np.linspace(-1.3, 1.3, 1001)
    for bits in (1, 3, 8, 12):
        np.testing.assert_allclose(c.pcm_decode(c.pcm_encode(x, bits, 1.3), bits, 1.3), quantize(x, -1.3, 1.3, bits), atol=1e-12)


def test_viterbi_equals_brute_force_maximum_likelihood():
    rng = np.random.default_rng(7)
    codebook = {k: c.convolutional_encode([(k >> i) & 1 for i in range(8)]) for k in range(256)}
    for _ in range(25):
        received = rng.integers(0, 2, codebook[0].size)
        best = min(int(np.count_nonzero(word != received)) for word in codebook.values())
        result = c.viterbi_decode(received)
        assert result.path_metric == best
        assert np.count_nonzero(c.convolutional_encode(result.data) != received) == best


def test_huffman_average_length_is_the_sum_of_merged_probabilities():
    # Every merge adds one bit to all the symbols beneath it, so L equals the sum of the merged weights.
    import heapq

    rng = np.random.default_rng(3)
    p = rng.dirichlet(np.ones(9))
    heap, merged = list(p), 0.0
    heapq.heapify(heap)
    while len(heap) > 1:
        weight = heapq.heappop(heap) + heapq.heappop(heap)
        merged += weight
        heapq.heappush(heap, weight)
    assert c.huffman_code(p).average_length == pytest.approx(merged)


def test_mpsk_and_qam_symbol_error_rates_by_simulation():
    rng = np.random.default_rng(1962)
    n = 100_000
    for points, ser in [(c.psk_constellation(8), c.ser_mpsk_awgn), (c.qam_constellation(16), c.ser_mqam_awgn)]:
        esn0 = 10 ** (13 / 10)
        labels = rng.integers(0, points.size, n)
        noise = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) * np.sqrt(1 / (2 * esn0))
        measured = np.mean(c.minimum_distance_detect(points[labels] + noise, points) != labels)
        assert measured == pytest.approx(ser(esn0, points.size), rel=0.1)


def test_rayleigh_ber_and_selection_outage_by_simulation():
    rng = np.random.default_rng(1880)
    n, mean = 400_000, 10.0
    h = (rng.standard_normal((2, n)) + 1j * rng.standard_normal((2, n))) / np.sqrt(2)
    bits = rng.integers(0, 2, n)
    noise = (rng.standard_normal((2, n)) + 1j * rng.standard_normal((2, n))) * np.sqrt(1 / (2 * mean))
    received = h * c.bpsk(bits) + noise
    single = np.mean((np.real(np.conj(h[0]) * received[0]) > 0) != bits)
    combined = np.mean((np.real(c.maximal_ratio_combine(received, h)) > 0) != bits)
    assert single == pytest.approx(c.ber_bpsk_rayleigh(mean), rel=0.05)
    assert combined == pytest.approx(c.ber_bpsk_mrc_rayleigh(mean, 2), rel=0.1)
    best = np.max(np.abs(h) ** 2, axis=0) * mean
    assert np.mean(best < 2.0) == pytest.approx(c.selection_outage_probability(2.0, mean, 2), rel=0.05)
