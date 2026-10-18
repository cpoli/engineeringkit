"""Analytical checks for the models added in the signal_processing rebuild."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from scipy import signal
from scipy.special import sici

from engineeringkit import signal_processing as s


def test_radix2_fft_matches_numpy_and_counts_butterflies():
    rng = np.random.default_rng(0)
    for n in (1, 2, 8, 256):
        x = rng.normal(size=n) + 1j * rng.normal(size=n)
        X, count = s.radix2_fft(x, return_multiplications=True)
        np.testing.assert_allclose(X, np.fft.fft(x), atol=1e-10)
        assert count == (n // 2) * int(np.log2(n))
    np.testing.assert_allclose(s.radix2_fft([1.0, 2.0]), [3.0, -1.0])
    for bad in ([1.0, 2.0, 3.0], np.ones((2, 2)), [np.nan, 1.0]):
        with pytest.raises(ValueError):
            s.radix2_fft(bad)


def test_fourier_partial_sum_gibbs_overshoot():
    k = np.arange(1, 400)
    b = np.where(k % 2 == 1, 4 / (np.pi * k), 0.0)
    t = np.linspace(0, 0.01, 4001)
    overshoot = s.fourier_partial_sum(t, 1.0, np.zeros_like(b), b).max() - 1
    # The square wave jumps by 2; the limit is 2 * (Si(pi)/pi - 1/2).
    assert np.isclose(overshoot, 2 * (sici(np.pi)[0] / np.pi - 0.5), rtol=1e-3)
    assert np.isclose(s.fourier_partial_sum(0.25, 1.0, [0.0], [2.0], a0=1.0), 3.0)
    with pytest.raises(ValueError):
        s.fourier_partial_sum(0.0, 1.0, [1.0], [1.0, 2.0])
    with pytest.raises(ValueError):
        s.fourier_partial_sum(0.0, 0.0, [1.0], [1.0])


def test_haar_transform_preserves_energy_and_inverts():
    rng = np.random.default_rng(1)
    x = rng.normal(size=64)
    coefficients = s.haar_transform(x)
    assert [c.size for c in coefficients] == [1, 1, 2, 4, 8, 16, 32]
    assert np.isclose(sum(np.sum(c**2) for c in coefficients), np.sum(x**2))
    assert np.isclose(coefficients[0][0], x.sum() / 8)
    np.testing.assert_allclose(s.inverse_haar_transform(coefficients), x)
    two = s.haar_transform(x[:12], levels=2)
    assert [c.size for c in two] == [3, 3, 6]
    np.testing.assert_allclose(s.inverse_haar_transform(two), x[:12])
    for bad, kwargs in (([1.0], {}), ([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], {"levels": 2}), ([1 + 1j, 2.0], {}), ([1.0, 2.0], {"levels": 0})):
        with pytest.raises(ValueError):
            s.haar_transform(bad, **kwargs)
    with pytest.raises(ValueError):
        s.inverse_haar_transform([[1.0]])
    with pytest.raises(ValueError):
        s.inverse_haar_transform([[1.0], [1.0, 2.0]])


def test_wiener_filter_identifies_an_fir_system():
    rng = np.random.default_rng(2)
    x = rng.normal(size=20000)
    h = np.array([0.8, -0.4, 0.2])
    d = signal.lfilter(h, 1, x) + 0.01 * rng.normal(size=x.size)
    w = s.wiener_filter(x, d, 5)
    np.testing.assert_allclose(w, [*h, 0, 0], atol=0.01)
    with pytest.raises(ValueError):
        s.wiener_filter(np.zeros(10), np.zeros(10), 2)
    for args in ((x[:5], d[:4], 2), (x[:5], d[:5], 6), (x[:5] + 0j, d[:5], 2), (x[:5], d[:5], 0)):
        with pytest.raises(ValueError):
            s.wiener_filter(*args)


def test_lms_converges_to_the_wiener_solution():
    rng = np.random.default_rng(3)
    x = rng.normal(size=20000)
    h = np.array([0.5, 0.25, -0.1])
    d = signal.lfilter(h, 1, x)
    result = s.lms_filter(x, d, 3, mu=0.01)
    np.testing.assert_allclose(result.weights, h, atol=1e-6)
    np.testing.assert_allclose(result.error, d - result.output)
    assert result.weight_history.shape == (x.size, 3)
    summary = result.summary()
    assert summary["taps"] == 3 and summary["final_mse"] < 1e-12 < summary["initial_mse"]
    started = s.lms_filter(x[:10], d[:10], 3, 0.01, initial_weights=h)
    np.testing.assert_allclose(started.weights, h)
    fig, ax = result.plot()
    assert ax.get_yscale() == "log"
    plt.close(fig)
    with pytest.raises(ValueError):
        s.lms_filter(x[:10], d[:10], 3, 0.01, initial_weights=[1.0, 2.0])
    with pytest.raises(ValueError):
        s.lms_filter(x[:10], d[:10], 3, 0.0)


def test_goertzel_matches_dft_bins_and_off_bin_sums():
    rng = np.random.default_rng(4)
    x = rng.normal(size=205)
    fs = 8000.0
    X = np.fft.fft(x)
    for k in (0, 18, 31):
        assert np.isclose(s.goertzel(x, k * fs / x.size, fs), X[k])
    f = np.array([697.0, 1209.0])
    direct = np.exp(-2j * np.pi * np.outer(f, np.arange(x.size)) / fs) @ x
    np.testing.assert_allclose(s.goertzel(x, f, fs), direct)
    with pytest.raises(ValueError):
        s.goertzel(np.ones((2, 2)), 1.0, 4.0)


def test_welch_psd_matches_scipy_and_integrates_to_the_variance():
    rng = np.random.default_rng(5)
    x = rng.normal(scale=2.0, size=8192)
    f, psd = s.welch_psd(x, 100.0, segment_length=256)
    f_ref, psd_ref = signal.welch(x, fs=100.0, nperseg=256, detrend=False)
    np.testing.assert_allclose(f, f_ref)
    np.testing.assert_allclose(psd, psd_ref)
    assert np.isclose(np.sum(psd) * (f[1] - f[0]), 4.0, rtol=0.05)
    _, odd = s.welch_psd(x[:999], 10.0, segment_length=99, overlap=0, window="boxcar")
    _, odd_ref = signal.welch(x[:999], fs=10.0, nperseg=99, noverlap=0, window="boxcar", detrend=False)
    np.testing.assert_allclose(odd, odd_ref)
    for kwargs in ({"segment_length": 1}, {"segment_length": 9000}, {"segment_length": 8, "overlap": 8}):
        with pytest.raises(ValueError):
            s.welch_psd(x, 100.0, **kwargs)


def test_stft_frames_and_summary_and_plot():
    fs = 1000.0
    t = np.arange(2000) / fs
    x = np.where(t < 1, np.sin(2 * np.pi * 100 * t), np.sin(2 * np.pi * 300 * t))
    result = s.stft(x, fs, 200, hop=50, window=("gaussian", 30))
    summary = result.summary()
    assert (summary["frames"], summary["bins"]) == (37, 101)
    assert np.isclose(summary["time_step"], 0.05) and np.isclose(summary["frequency_resolution"], 5.0)
    peaks = result.frequency[np.argmax(np.abs(result.values), axis=0)]
    assert peaks[0] == 100 and peaks[-1] == 300
    np.testing.assert_allclose(result.values[:, 0], np.fft.rfft(x[:200] * signal.get_window(("gaussian", 30), 200)))
    assert s.stft(x[:400], fs, 400).summary()["time_step"] == 0.0
    fig, ax = result.plot()
    assert ax.get_ylabel() == "frequency [Hz]"
    plt.close(fig)
    for kwargs in ({"window_length": 1}, {"window_length": 3000}, {"window_length": 8, "hop": 0}):
        with pytest.raises(ValueError):
            s.stft(x, fs, **kwargs)


def test_spectrum_summary_and_single_bin():
    spectrum = s.Spectrum.from_signal(np.cos(2 * np.pi * 2 * np.arange(16) / 16), 16.0)
    assert spectrum.summary() == {"bins": 9, "resolution": 1.0, "peak_frequency": 2.0}
    assert s.Spectrum.from_signal([1.0], 1.0).summary()["resolution"] == 0.0


def test_kaiser_order_formulas_and_design_meets_the_attenuation():
    for attenuation in (60.0, 40.0):
        numtaps, beta = s.kaiser_order(attenuation, 50.0, 1000.0)
        reference_taps, reference_beta = signal.kaiserord(attenuation, 2 * 50.0 / 1000.0)
        assert numtaps == reference_taps and np.isclose(beta, reference_beta)
        h = s.fir_lowpass(numtaps, 200.0, 1000.0, window=("kaiser", beta))
        f, H = s.frequency_response(h, 1, worN=np.linspace(225.0, 500.0, 500), fs=1000.0)
        assert 20 * np.log10(np.abs(H).max()) < -attenuation + 1.0
    assert s.kaiser_order(15.0, 50.0, 1000.0)[1] == 0.0
    with pytest.raises(ValueError):
        s.kaiser_order(60.0, 600.0, 1000.0)


def test_second_order_sections_survive_where_polynomials_fail():
    fs, cutoff = 1000.0, 5.0
    sos = s.butterworth_lowpass(12, cutoff, fs, output="sos")
    assert sos.shape == (6, 6)
    f, H = s.frequency_response(sos, None, worN=[cutoff], fs=fs)
    assert np.isclose(abs(H[0]), 1 / np.sqrt(2), rtol=1e-6)
    b, a = s.butterworth_lowpass(12, cutoff, fs)
    # Rounding the order-12 polynomial pushes a pole outside the unit circle.
    assert np.abs(np.roots(a)).max() > 1
    step = s.filter_signal(sos, None, np.ones(5000))
    assert np.isclose(step[-1], 1.0, atol=1e-6)
    hp = s.butterworth_highpass(4, 50.0, fs, output="sos")
    assert abs(s.frequency_response(hp, None, worN=[0.0], fs=fs)[1][0]) < 1e-12
    for bad in (np.ones((2, 5)), np.ones((0, 6)), np.full((1, 6), np.nan)):
        with pytest.raises(ValueError):
            s.filter_signal(bad, None, [1.0])
    with pytest.raises(ValueError):
        s.butterworth_lowpass(2, 600.0, 1000.0)
    with pytest.raises(ValueError):
        s.butterworth_lowpass(2, 100.0, 1000.0, output="zpk")
