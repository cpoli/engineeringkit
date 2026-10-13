import numpy as np
import pytest

from engineeringkit import signal_processing as s


def test_parseval_energy_with_complex_signal():
    x = np.exp(2j * np.pi * np.arange(16) / 16) + 0.2 * np.arange(16)
    assert np.isclose(s.energy(x), s.energy(s.dft(x)) / len(x))
    assert np.isclose(s.rms(x) ** 2, s.power(x))


def test_sampling_reconstruction_at_sample_locations():
    samples = np.array([0.0, 1.0, -0.5, 2.0])
    fs = 10.0
    np.testing.assert_allclose(s.sinc_reconstruct(samples, fs, np.arange(4) / fs), samples, atol=1e-14)


def test_butterworth_cutoff_and_fir_dc_gain():
    b, a = s.butterworth_lowpass(3, 100, 1000)
    _, response = s.frequency_response(b, a, worN=[100], fs=1000)
    assert np.isclose(abs(response[0]), 1 / np.sqrt(2), rtol=1e-12)
    assert np.isclose(s.fir_lowpass(31, 100, 1000).sum(), 1)


def test_psd_integrates_to_variance():
    fs = 1000
    n = 1000
    x = np.sin(2 * np.pi * 50 * np.arange(n) / fs)
    frequencies, psd = s.power_spectral_density(x, fs)
    assert np.isclose(psd.sum() * (frequencies[1] - frequencies[0]), 0.5)


def test_snr_degenerate_and_shape_cases():
    assert s.snr_db([1], [0]) == np.inf
    assert s.snr_db([0], [1]) == -np.inf
    with pytest.raises(ValueError):
        s.snr_db([0], [0])
    with pytest.raises(ValueError):
        s.snr_db([1, 1], [1])
    with pytest.raises(ValueError):
        s.moving_average([1, 2], 3)
