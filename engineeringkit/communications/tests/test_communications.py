import numpy as np
import pytest
from scipy import integrate
from scipy.integrate import trapezoid
from scipy.special import erfc

import engineeringkit.communications as m


def test_capacity_and_db():
    assert np.isclose(m.from_db(m.db(10)), 10)
    assert np.isclose(m.shannon_capacity(1, 1), 1)


def test_ber():
    assert m.ber_bpsk_awgn(10) < m.ber_bpsk_awgn(1)


# ---------------------------------------------------------------------------
# Curriculum models (ROADMAP D2): analytical checks.
# ---------------------------------------------------------------------------


def test_matched_filter_reaches_the_cauchy_schwarz_bound():
    rng = np.random.default_rng(43)
    pulse = rng.standard_normal(32) + 1j * rng.standard_normal(32)
    h = m.matched_filter(pulse)
    energy, N0 = np.sum(np.abs(pulse) ** 2), 0.3
    peak = np.convolve(pulse, h)[pulse.size - 1]
    assert np.isclose(peak, energy)

    def output_snr(filt):
        return np.abs(np.convolve(pulse, filt)[pulse.size - 1]) ** 2 / (N0 / 2 * np.sum(np.abs(filt) ** 2))

    assert np.isclose(output_snr(h), m.matched_filter_snr(energy, N0))
    for _ in range(20):
        assert output_snr(rng.standard_normal(32) + 1j * rng.standard_normal(32)) < m.matched_filter_snr(energy, N0)


CONSTELLATIONS = [m.psk_constellation(8), m.psk_constellation(16, np.pi / 16), m.qam_constellation(16), m.qam_constellation(64)]


@pytest.mark.parametrize("points", CONSTELLATIONS)
def test_constellations_have_unit_energy_and_gray_neighbours(points):
    assert np.isclose(np.mean(np.abs(points) ** 2), 1)
    distance = np.abs(points[:, None] - points[None, :])
    nearest = np.isclose(distance, distance[distance > 1e-9].min())
    labels = np.arange(points.size)
    flips = np.vectorize(lambda a, b: bin(a ^ b).count("1"))(labels[:, None], labels[None, :])
    assert np.all(flips[nearest] == 1)


def test_four_point_constellations_coincide_with_qpsk():
    np.testing.assert_allclose(np.sort_complex(m.psk_constellation(4, np.pi / 4)), np.sort_complex(m.qam_constellation(4)), atol=1e-15)
    np.testing.assert_allclose(np.sort_complex(m.qam_constellation(4)), np.sort_complex(m.qpsk([0, 0, 0, 1, 1, 0, 1, 1])), atol=1e-15)


def test_symbol_error_rates_match_closed_forms():
    g = np.array([0.0, 0.5, 2.0, 10.0])
    np.testing.assert_allclose(m.ser_mpsk_awgn(g, 2), m.ber_bpsk_awgn(g), rtol=1e-9)
    q = m.q_function(np.sqrt(2 * g) * np.sin(np.pi / 4))
    np.testing.assert_allclose(m.ser_mpsk_awgn(g, 4), 2 * q - q**2, rtol=1e-9, atol=1e-15)
    np.testing.assert_allclose(m.ser_mqam_awgn(g, 4), 2 * q - q**2, rtol=1e-12)
    assert m.ser_mpsk_awgn(0.0, 8) == pytest.approx(7 / 8)
    assert m.ser_mqam_awgn(0.0, 16) == pytest.approx(1 - 1 / 16)
    high = 10 ** (20 / 10)
    assert m.ser_mpsk_awgn(high, 8) == pytest.approx(2 * m.q_function(np.sqrt(2 * high) * np.sin(np.pi / 8)), rel=1e-3)


def test_minimum_distance_detection_without_noise_returns_the_labels():
    points = m.qam_constellation(16)
    labels = np.arange(16).reshape(4, 4)
    np.testing.assert_array_equal(m.minimum_distance_detect(points[labels], points), labels)


def test_raised_cosine_satisfies_nyquist_criterion():
    T = 1e-3
    for a in (0.0, 0.25, 0.5, 1.0):
        k = np.arange(-20, 21)
        np.testing.assert_allclose(m.raised_cosine_pulse(k * T, T, a), (k == 0).astype(float), atol=1e-12)
        f = np.linspace(-0.5 / T, 0.5 / T, 101)
        folded = sum(m.raised_cosine_spectrum(f + n / T, T, a) for n in range(-3, 4))
        np.testing.assert_allclose(folded, T, rtol=1e-12)
    assert m.raised_cosine_pulse(T / 0.5, T, 0.25) == pytest.approx(np.pi / 4 * np.sinc(2), abs=1e-15)
    assert m.raised_cosine_pulse(T / 0.5 + 1e-12, T, 0.25) == pytest.approx(np.pi / 4 * np.sinc(2), abs=1e-9)
    assert m.raised_cosine_bandwidth(1 / T, 0.5) == pytest.approx(750)


def test_raised_cosine_pulse_is_the_inverse_transform_of_its_spectrum():
    T, a = 1.0, 0.35
    f = np.linspace(-1, 1, 20001)
    for t in (0.3, 1.0 / (2 * a), 2.7):
        transform = trapezoid(m.raised_cosine_spectrum(f, T, a) * np.cos(2 * np.pi * f * t), f)
        assert transform == pytest.approx(m.raised_cosine_pulse(t, T, a), abs=1e-6)


def test_root_raised_cosine_has_unit_energy_and_squares_to_raised_cosine():
    T = 1.0
    dt = T / 100
    t = np.arange(-60 * T, 60 * T + dt / 2, dt)
    for a in (0.25, 0.5, 1.0):
        h = m.root_raised_cosine_pulse(t, T, a)
        assert np.sum(h**2) * dt == pytest.approx(1, abs=1e-4)
        full = np.convolve(h, h)[t.size // 2 : t.size // 2 + t.size] * dt
        centre = np.abs(t) < 10 * T
        np.testing.assert_allclose(full[centre], m.raised_cosine_pulse(t[centre], T, a), atol=2e-4)
    quarter = 1 / (4 * 0.5)
    assert m.root_raised_cosine_pulse(quarter, T, 0.5) == pytest.approx(m.root_raised_cosine_pulse(quarter + 1e-7, T, 0.5), abs=1e-5)
    assert m.root_raised_cosine_pulse(0.0, T, 0.5) == pytest.approx(1 - 0.5 + 2 / np.pi)
    np.testing.assert_allclose(m.root_raised_cosine_pulse(t[:5], T, 0.0), np.sinc(t[:5]))


def test_pcm_round_trip_error_is_within_half_a_step():
    x = np.linspace(-2, 2, 401)
    for N in (1, 4, 8):
        y = m.pcm_decode(m.pcm_encode(x, N, 2.0), N, 2.0)
        q = 4.0 / 2**N
        assert np.max(np.abs(y - x)) <= q / 2 + 1e-12
        assert np.unique(y).size == min(2**N, x.size)
    assert m.pcm_encode([5.0, -5.0], 2).tolist() == [1, 1, 0, 0]


def test_pcm_sqnr_follows_six_db_per_bit():
    t = np.arange(200_000)
    sine = np.sin(2 * np.pi * 0.1234567 * t)
    shortfall = []
    for N in (6, 8, 10, 12):
        noise = m.pcm_decode(m.pcm_encode(sine, N), N) - sine
        shortfall.append(m.pcm_sqnr_db(N) - 10 * np.log10(0.5 / np.mean(noise**2)))
    np.testing.assert_allclose(shortfall, 2.3 * 2.0 ** (-np.array([6, 8, 10, 12]) / 2), rtol=0.15)
    assert m.pcm_sqnr_db(9) - m.pcm_sqnr_db(8) == pytest.approx(20 * np.log10(2))


def test_mu_law_is_invertible_and_odd():
    x = np.linspace(-3, 3, 61)
    y = m.mu_law_compress(x, 100, 3)
    np.testing.assert_allclose(m.mu_law_expand(y, 100, 3), x, atol=1e-12)
    np.testing.assert_allclose(y, -m.mu_law_compress(-x, 100, 3))
    assert m.mu_law_compress(3.0, 100, 3) == pytest.approx(3)
    assert m.mu_law_compress(0.03, 100, 3) > 10 * 0.03


def test_delta_modulation_tracks_within_one_step_without_overload():
    fs, step = 8000.0, 0.05
    t = np.arange(800) / fs
    A, f = 1.0, 0.9 * m.delta_slope_limit(step, fs) / (2 * np.pi)
    x = A * np.sin(2 * np.pi * f * t)
    bits = m.delta_modulate(x, step)
    staircase = m.delta_demodulate(bits, step)
    assert np.max(np.abs(staircase[:-1] - x[1:])) <= 2 * step
    overload = m.delta_demodulate(m.delta_modulate(3 * x, step), step)
    assert np.max(np.abs(overload[:-1] - 3 * x[1:])) > 10 * step
    np.testing.assert_allclose(m.delta_demodulate(m.delta_modulate(x, step, 0.5), step, 0.5)[0], 0.5 - step)


def test_entropy_generalises_the_binary_entropy():
    assert m.entropy(np.full(8, 1 / 8)) == pytest.approx(3)
    assert m.entropy([0.11, 0.89]) == pytest.approx(m.entropy_binary(0.11))
    assert m.entropy([1.0, 0.0]) == 0


def test_huffman_code_is_optimal_prefix_free_and_reversible():
    dyadic = m.huffman_code([0.5, 0.25, 0.125, 0.125])
    assert dyadic.lengths.tolist() == [1, 2, 3, 3]
    assert dyadic.average_length == pytest.approx(dyadic.entropy)
    rng = np.random.default_rng(1952)
    for size in (2, 5, 12, 40):
        p = rng.dirichlet(np.ones(size))
        code = m.huffman_code(p)
        assert code.entropy <= code.average_length + 1e-12 < code.entropy + 1
        assert sum(2.0**-code.lengths) == pytest.approx(1)  # Kraft equality: a full binary tree
        words = code.codewords
        assert not any(a != b and b.startswith(a) for a in words for b in words)
        symbols = rng.integers(0, size, 300)
        np.testing.assert_array_equal(code.decode(code.encode(symbols)), symbols)
    summary = m.huffman_code([0.4, 0.2, 0.2, 0.1, 0.1]).summary()
    assert summary["average_length"] == pytest.approx(2.2)
    assert summary["efficiency"] == pytest.approx(m.entropy([0.4, 0.2, 0.2, 0.1, 0.1]) / 2.2)


def test_crc_detects_every_single_error_and_every_short_burst():
    rng = np.random.default_rng(1961)
    message = rng.integers(0, 2, 64)
    word = m.crc_encode(message, m.CRC16_CCITT)
    assert word.size == 80 and m.crc_check(word, m.CRC16_CCITT)
    for length in range(1, 17):
        for start in range(0, word.size - length + 1, 7):
            burst = word.copy()
            burst[start] ^= 1
            burst[start + length - 1] ^= 1 if length > 1 else 0
            burst[start + 1 : start + length - 1] ^= rng.integers(0, 2, max(length - 2, 0))
            assert not m.crc_check(burst, m.CRC16_CCITT)


def test_crc16_ccitt_check_value():
    message = np.unpackbits(np.frombuffer(b"123456789", dtype=np.uint8))
    check = m.crc_remainder(message, m.CRC16_CCITT)
    assert int("".join(map(str, check)), 2) == 0x31C3  # CRC-16/XMODEM catalogue check value


def test_convolutional_code_is_linear_and_has_its_free_distance():
    rng = np.random.default_rng(1955)
    a, b = rng.integers(0, 2, 30), rng.integers(0, 2, 30)
    np.testing.assert_array_equal(m.convolutional_encode(a ^ b), m.convolutional_encode(a) ^ m.convolutional_encode(b))
    assert m.convolutional_encode([1]).sum() == m.convolutional_free_distance() == 5
    assert m.convolutional_free_distance((0o171, 0o133)) == 10
    assert m.convolutional_free_distance((0o15, 0o17)) == 6
    assert m.convolutional_encode(a, (0o15, 0o17, 0o13)).size == 3 * (30 + 3)


def test_viterbi_corrects_errors_spread_beyond_the_free_distance():
    rng = np.random.default_rng(1967)
    data = rng.integers(0, 2, 120)
    code = m.convolutional_encode(data)
    assert np.array_equal(m.viterbi_decode(code).data, data)
    received = code.copy()
    received[::20] ^= 1  # one error in every 20 code bits: each alone within any free-distance span
    result = m.viterbi_decode(received)
    np.testing.assert_array_equal(result.data, data)
    assert result.path_metric == np.count_nonzero(received != code)
    assert result.summary() == {"bits": 120, "path_metric": result.path_metric}
    soft = 2.0 * code - 1 + 0.45 * rng.standard_normal(code.size)
    np.testing.assert_array_equal(m.viterbi_decode(soft, soft=True).data, data)
    np.testing.assert_array_equal(m.viterbi_decode(m.convolutional_encode(data, (0b11, 0b10)), (0b11, 0b10)).data, data)


def test_rayleigh_density_and_fading_error_rates():
    r = np.linspace(0, 12, 20001)
    assert trapezoid(m.rayleigh_pdf(r, 4.0), r) == pytest.approx(1, abs=1e-7)
    assert trapezoid(r**2 * m.rayleigh_pdf(r, 4.0), r) == pytest.approx(4, abs=1e-6)
    g = np.array([1.0, 10.0, 1000.0])
    np.testing.assert_allclose(m.ber_bpsk_mrc_rayleigh(g, 1), m.ber_bpsk_rayleigh(g))
    assert m.ber_bpsk_rayleigh(1e4) == pytest.approx(1 / (4e4), rel=1e-3)
    assert m.ber_bpsk_rayleigh(0.0) == 0.5
    ratio = m.ber_bpsk_mrc_rayleigh(1e3, 3) / m.ber_bpsk_mrc_rayleigh(1e4, 3)
    assert ratio == pytest.approx(1e3, rel=0.01)  # diversity order three


def test_mrc_ber_equals_the_average_over_the_combined_snr():
    # The combined Eb/N0 of L i.i.d. Rayleigh branches is gamma-distributed with shape L.
    from scipy.stats import gamma

    for L, mean in [(2, 3.0), (4, 10.0)]:
        average, _ = integrate.quad(lambda s: 0.5 * erfc(np.sqrt(s)) * gamma.pdf(s, L, scale=mean), 0, np.inf)
        assert m.ber_bpsk_mrc_rayleigh(mean, L) == pytest.approx(average, rel=1e-7)


def test_selection_outage_and_maximal_ratio_combining():
    assert m.selection_outage_probability(1.0, 10.0) == pytest.approx(1 - np.exp(-0.1))
    assert m.selection_outage_probability(1.0, 10.0, 3) == pytest.approx((1 - np.exp(-0.1)) ** 3)
    h = np.array([[0.3 + 0.4j, 1.0], [-1j, 0.0]])
    s = np.array([1 + 1j, -1 + 0j])
    np.testing.assert_allclose(m.maximal_ratio_combine(h * s, h), s)
    np.testing.assert_allclose(m.maximal_ratio_combine([[2.0, -2.0]], [2.0]), [1, -1])


def test_antenna_aperture_relations():
    f = 3e9
    lam = m.wavelength(f)
    assert m.effective_aperture(1.0, f) == pytest.approx(lam**2 / (4 * np.pi))
    assert m.aperture_gain(m.effective_aperture(40.0, f), f) == pytest.approx(40)
    # Friis in effective areas: Pr/Pt = Ae_t Ae_r / (d lam)^2.
    Gt, Gr, d = 30.0, 20.0, 5e3
    areas = m.effective_aperture(Gt, f) * m.effective_aperture(Gr, f)
    assert m.friis_received_power(1.0, Gt, Gr, d, f) == pytest.approx(areas / (d * lam) ** 2)


def test_half_wave_dipole_against_numerical_integration():
    theta = np.linspace(0, np.pi, 20001)
    F = m.half_wave_dipole_pattern(theta)
    assert m.directivity(theta, F**2) == pytest.approx(m.half_wave_dipole_directivity(), rel=1e-6)
    assert 10 * np.log10(m.half_wave_dipole_directivity()) == pytest.approx(2.15, abs=0.005)
    integral, _ = integrate.quad(lambda th: np.cos(np.pi / 2 * np.cos(th)) ** 2 / np.sin(th), 1e-12, np.pi - 1e-12)
    assert m.half_wave_dipole_resistance() == pytest.approx(m.FREE_SPACE_IMPEDANCE / (2 * np.pi) * integral, rel=1e-6)
    assert m.half_wave_dipole_resistance() == pytest.approx(73.08, abs=0.01)
    assert m.half_wave_dipole_resistance(377.0) / m.half_wave_dipole_resistance(376.730313668) == pytest.approx(377.0 / 376.730313668)
    assert m.half_wave_dipole_pattern(np.pi / 2) == 1.0 and m.half_wave_dipole_pattern(np.pi) == 0.0
    # Its gain times the effective aperture formula matches the classic 0.13 lambda^2.
    assert m.effective_aperture(m.half_wave_dipole_directivity(), 1e8) / m.wavelength(1e8) ** 2 == pytest.approx(0.1306, abs=1e-4)


def test_am_and_fm_noise_performance():
    gamma = 1000.0
    assert m.am_output_snr(gamma) == pytest.approx(gamma / 3)
    assert m.am_output_snr(gamma, mu=0.0) == 0
    assert m.fm_output_snr(gamma, 1.0) / m.am_output_snr(gamma) == pytest.approx(4.5)
    assert m.fm_output_snr(gamma, 5.0) / m.fm_output_snr(gamma, 2.5) == pytest.approx(4)
    assert m.fm_threshold_snr(5.0, 10.0) == pytest.approx(120)
    assert m.baseband_snr(2e-12, 4e-21, 5e3) == pytest.approx(1e5)
    bt = m.carson_bandwidth(5.0 * 15e3, 15e3)
    assert m.fm_threshold_snr(5.0) == pytest.approx(10 * bt / 15e3)
