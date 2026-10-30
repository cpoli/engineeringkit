import numpy as np
import pytest

from engineeringkit import communications as c


def test_units_noise_rates_and_symbol_count():
    assert c.wavelength(299792458) == 1
    assert c.dbm(0.001) == 0
    assert np.isclose(c.thermal_noise_power(1, 1), 1.380649e-23)
    assert c.noise_figure(10) == 10
    assert c.hartley_capacity(100, 4) == 400
    assert c.spectral_efficiency(400, 100) == 4
    assert c.symbol_rate(400, 4) == 200
    assert c.nyquist_symbol_rate(100) == 200
    assert c.ber_qpsk_awgn(1) == c.ber_bpsk_awgn(1)
    assert c.q_function(0) == 0.5
    assert c.hamming_distance([0, 1, 1], [0, 0, 1]) == 1


def test_am_and_suppressed_carrier_waveforms():
    t = np.array([0, 0.25, 0.5])
    message = np.ones(3) * 0.5
    np.testing.assert_allclose(c.am_signal(message, t, 1, mu=0.5, carrier_amplitude=2), [2.5, 0, -2.5], atol=1e-15)
    np.testing.assert_allclose(c.dsb_sc(message, t, 1), [0.5, 0, -0.5], atol=1e-15)
    with pytest.raises(ValueError):
        c.qpsk([0, 1, 1])
    with pytest.raises(ValueError):
        c.cascaded_noise_factor([10], [0.5])


@pytest.mark.parametrize(
    "fn, args",
    [
        (c.matched_filter, ([],)),
        (c.matched_filter, ([[1.0]],)),
        (c.matched_filter_snr, (1.0, 0.0)),
        (c.minimum_distance_detect, ([], [1, -1])),
        (c.minimum_distance_detect, ([np.nan], [1, -1])),
        (c.ser_mpsk_awgn, (-1.0, 4)),
        (c.ser_mqam_awgn, (1.0, 8)),
        (c.ser_mqam_awgn, (1.0, 12)),
        (c.psk_constellation, (6,)),
        (c.psk_constellation, (8, np.inf)),
        (c.qam_constellation, (8,)),
        (c.qam_constellation, (2,)),
        (c.raised_cosine_pulse, (0.0, 0.0, 0.5)),
        (c.raised_cosine_spectrum, (0.0, 1.0, 1.5)),
        (c.root_raised_cosine_pulse, (np.nan, 1.0, 0.5)),
        (c.raised_cosine_bandwidth, (-1.0, 0.5)),
        (c.pcm_encode, ([0.1], 0)),
        (c.pcm_encode, ([0.1], True)),
        (c.pcm_encode, ([], 3)),
        (c.pcm_encode, ([0.1], 3, 0.0)),
        (c.pcm_decode, ([1, 0], 3)),
        (c.pcm_decode, ([1, 0, 2], 3)),
        (c.pcm_sqnr_db, (53,)),
        (c.mu_law_compress, (1.5,)),
        (c.mu_law_expand, (0.5, 0.0)),
        (c.delta_modulate, ([[0.0]], 0.1)),
        (c.delta_modulate, ([0.0], 0.0)),
        (c.delta_demodulate, ([], 0.1)),
        (c.delta_slope_limit, (0.1, 0.0)),
        (c.entropy, ([1.0],)),
        (c.entropy, ([0.5, 0.6],)),
        (c.huffman_code, ([-0.5, 1.5],)),
        (c.crc_remainder, ([1, 0], 1)),
        (c.crc_check, ([1, 2], 0b111)),
        (c.convolutional_encode, ([1, 0], 7)),
        (c.convolutional_encode, ([1, 0], (7,))),
        (c.convolutional_encode, ([1, 0], (7, 0))),
        (c.convolutional_encode, ([1, 0], (1, 1))),
        (c.convolutional_encode, ([1, 2],)),
        (c.viterbi_decode, ([1, 0, 1],)),
        (c.viterbi_decode, ([1.0, 0.5, 0.2], (7, 5), True)),
        (c.viterbi_decode, ([1, 1, 1, 0],)),
        (c.rayleigh_pdf, (-1.0,)),
        (c.ber_bpsk_rayleigh, (np.nan,)),
        (c.ber_bpsk_mrc_rayleigh, (1.0, 0)),
        (c.selection_outage_probability, (1.0, 0.0)),
        (c.maximal_ratio_combine, ([1.0, 2.0], [1.0])),
        (c.maximal_ratio_combine, ([[1.0, 2.0]], [1.0, 1.0])),
        (c.maximal_ratio_combine, ([[1.0, 2.0]], [[1.0, 1.0, 1.0]])),
        (c.maximal_ratio_combine, ([[1.0], [2.0]], [0.0, 0.0])),
        (c.effective_aperture, (0.0, 1e9)),
        (c.aperture_gain, (1.0, 1e9, 1.2)),
        (c.directivity, ([0, 1, 2], [1, 1, 1])),
        (c.directivity, ([0, np.pi], [1, 1])),
        (c.directivity, (np.linspace(0, np.pi, 5), [1, 1, 1, 1])),
        (c.directivity, (np.linspace(0, np.pi, 5), np.zeros(5))),
        (c.directivity, (np.linspace(0, np.pi, 5), -np.ones(5))),
        (c.half_wave_dipole_pattern, (np.nan,)),
        (c.half_wave_dipole_resistance, (0.0,)),
        (c.baseband_snr, (1.0, 0.0, 1.0)),
        (c.am_output_snr, (10.0, 1.5)),
        (c.fm_output_snr, (10.0, -1.0)),
        (c.fm_threshold_snr, (5.0, 0.0)),
    ],
)
def test_curriculum_models_reject_invalid_input(fn, args):
    with pytest.raises(ValueError):
        fn(*args)


def test_huffman_code_rejects_bad_streams():
    code = c.huffman_code([0.5, 0.25, 0.25])
    for bad in ([], [0, 3], np.array([0.5])):
        with pytest.raises(ValueError):
            code.encode(bad)
    for bad in ([], [0, 2], [1]):
        with pytest.raises(ValueError):
            code.decode(bad)


def test_curriculum_models_return_scalars_and_arrays():
    assert isinstance(c.raised_cosine_spectrum(0.1, 1.0, 0.5), float)
    assert c.raised_cosine_spectrum(np.zeros(3), 1.0, 0.0).shape == (3,)
    assert c.raised_cosine_spectrum(0.5, 1.0, 0.0) == 0.5
    assert isinstance(c.root_raised_cosine_pulse(0.3, 1.0, 0.5), float)
    assert c.raised_cosine_pulse(np.linspace(-3, 3, 7), 1.0, 0.0).shape == (7,)
    assert isinstance(c.ser_mpsk_awgn(1.0, 8), float)
    assert c.ser_mpsk_awgn(np.ones((2, 2)), 8).shape == (2, 2)
    assert c.half_wave_dipole_pattern(np.array([0.0, np.pi / 2])).tolist() == [0.0, 1.0]
