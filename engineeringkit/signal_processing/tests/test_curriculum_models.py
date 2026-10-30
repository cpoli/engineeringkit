"""Analytical and benchmark checks for the models added by ROADMAP item D10 (curriculum gaps)."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from scipy import fft as sfft
from scipy import signal

from engineeringkit import signal_processing as s


def _gain_db(b, a, f, fs):
    return 20 * np.log10(np.abs(s.frequency_response(b, a, worN=np.atleast_1d(f), fs=fs)[1]))


def test_chebyshev_and_elliptic_designs_meet_their_edges():
    fs = 1000.0
    b, a = s.chebyshev1_lowpass(5, 0.5, 100.0, fs)
    f = np.linspace(0, 100.0, 2001)
    gain = _gain_db(b, a, f, fs)
    # Equiripple passband: the gain stays in [-R, 0] and touches both ends.
    assert gain.max() < 1e-9 and np.isclose(gain.min(), -0.5, atol=1e-6)
    assert np.isclose(gain[0], 0.0, atol=1e-9) and np.isclose(gain[-1], -0.5)
    b2, a2 = s.chebyshev2_lowpass(5, 50.0, 150.0, fs)
    stop = _gain_db(b2, a2, np.linspace(150.0, 499.0, 2000), fs)
    assert stop.max() < -50 + 1e-6 and np.isclose(stop[0], -50.0)
    assert np.isclose(_gain_db(b2, a2, 0.0, fs)[0], 0.0, atol=1e-9)
    sos = s.elliptic_lowpass(6, 0.5, 70.0, 100.0, fs, output="sos")
    assert np.isclose(_gain_db(sos, None, 100.0, fs)[0], -0.5)
    beyond = _gain_db(sos, None, np.linspace(100.0, 499.0, 4000), fs)
    # Once the gain first reaches -70 dB it never comes back above it: equiripple stopband.
    assert beyond[np.argmax(beyond <= -70) :].max() < -70 + 1e-6
    for call in (
        lambda: s.chebyshev1_lowpass(4, 0.0, 100.0, fs),
        lambda: s.chebyshev2_lowpass(4, -1.0, 100.0, fs),
        lambda: s.elliptic_lowpass(4, 1.0, 1.0, 100.0, fs),
        lambda: s.elliptic_lowpass(4, 1.0, 40.0, 600.0, fs),
    ):
        with pytest.raises(ValueError):
            call()


def test_iir_order_matches_scipy_and_the_order_meets_the_specification():
    fs = 8000.0
    references = {"butterworth": signal.buttord, "chebyshev": signal.cheb1ord, "elliptic": signal.ellipord}
    for f_pass, f_stop, ripple, attenuation in ((1000.0, 1500.0, 1.0, 60.0), (300.0, 350.0, 0.1, 40.0), (2000.0, 3500.0, 3.0, 20.0)):
        orders = {}
        for family, reference in references.items():
            orders[family] = s.iir_order(f_pass, f_stop, ripple, attenuation, fs, family)
            assert orders[family] == reference(f_pass, f_stop, ripple, attenuation, fs=fs)[0]
        assert orders["elliptic"] <= orders["chebyshev"] <= orders["butterworth"]
    n = s.iir_order(1000.0, 1500.0, 1.0, 60.0, fs, "elliptic")
    sos = s.elliptic_lowpass(n, 1.0, 60.0, 1000.0, fs, output="sos")
    assert _gain_db(sos, None, np.linspace(1500.0, 3999.0, 2000), fs).max() < -60 + 1e-6
    for args in ((150.0, 100.0, 1.0, 60.0, 1000.0), (100.0, 600.0, 1.0, 60.0, 1000.0), (100.0, 150.0, 60.0, 60.0, 1000.0)):
        with pytest.raises(ValueError):
            s.iir_order(*args)
    with pytest.raises(ValueError):
        s.iir_order(100.0, 150.0, 1.0, 60.0, 1000.0, "bessel")


def test_bilinear_transform_with_prewarping():
    fs = 1000.0
    b, a = signal.butter(4, 2 * np.pi * 100.0, analog=True)
    bz, az = s.bilinear_transform(b, a, fs)
    ref = signal.bilinear(b, a, fs=fs)
    np.testing.assert_allclose(bz, ref[0], rtol=1e-12)
    np.testing.assert_allclose(az, ref[1], rtol=1e-12)
    # Without prewarping the -3 dB point moves down to (fs/pi) arctan(pi f/fs).
    warped = fs / np.pi * np.arctan(np.pi * 100.0 / fs)
    assert np.isclose(_gain_db(bz, az, warped, fs)[0], -10 * np.log10(2), atol=1e-9)
    # Prewarping at the cutoff puts it back, and reproduces SciPy's digital Butterworth.
    bp, ap = s.bilinear_transform(b, a, fs, prewarp=100.0)
    assert np.isclose(_gain_db(bp, ap, 100.0, fs)[0], -10 * np.log10(2), atol=1e-9)
    bd, ad = signal.butter(4, 100.0, fs=fs)
    np.testing.assert_allclose(bp, bd, atol=1e-14)
    np.testing.assert_allclose(ap, ad, rtol=1e-10)
    # Prewarping designs the analog prototype at the warped frequency.
    fa = s.prewarp_frequency(100.0, fs)
    b2, a2 = signal.butter(4, 2 * np.pi * fa, analog=True)
    np.testing.assert_allclose(s.bilinear_transform(b2, a2, fs)[1], ad, rtol=1e-10)
    np.testing.assert_allclose(s.prewarp_frequency([0.0, 250.0], 1000.0), [0.0, 1000.0 / np.pi])
    # A stable analog pole maps inside the unit circle.
    assert np.all(np.abs(np.roots(s.bilinear_transform([1.0], [1.0, 0.1, 1e4], fs)[1])) < 1)
    for args, kwargs in (
        (([1.0, 1.0, 1.0], [1.0, 1.0], fs), {}),
        (([1.0], [0.0, 1.0], fs), {}),
        (([1.0], [1.0, 1.0], fs), {"prewarp": 600.0}),
        (([1.0], [1.0, 1.0], fs), {"prewarp": 0.0}),
    ):
        with pytest.raises(ValueError):
            s.bilinear_transform(*args, **kwargs)
    with pytest.raises(ValueError):
        s.prewarp_frequency(500.0, 1000.0)


def test_parks_mcclellan_equiripple_beats_the_window_method():
    fs = 1000.0
    h = s.equiripple_lowpass(51, 100.0, 150.0, fs, weight=10.0)
    f, H = s.frequency_response(h, 1, worN=8192, fs=fs)
    passband = np.abs(np.abs(H[f <= 100.0]) - 1)
    stopband = np.abs(H[f >= 150.0])
    # Weight 10 makes the stopband deviation a tenth of the passband one.
    assert np.isclose(passband.max() / stopband.max(), 10.0, rtol=0.02)
    # The window design of the same length has a worse stopband for a no-better passband.
    hw = s.fir_lowpass(51, 125.0, fs, window="hamming")
    _, Hw = s.frequency_response(hw, 1, worN=8192, fs=fs)
    assert np.abs(Hw[f >= 150.0]).max() > stopband.max()
    n = s.equiripple_order(0.01, 0.001, 50.0, fs)
    h = s.equiripple_lowpass(n + 4, 100.0, 150.0, fs, weight=10.0)
    _, H = s.frequency_response(h, 1, worN=8192, fs=fs)
    assert np.abs(H[f >= 150.0]).max() < 0.001
    for call in (
        lambda: s.equiripple_lowpass(41, 150.0, 100.0, fs),
        lambda: s.equiripple_lowpass(2, 100.0, 150.0, fs),
        lambda: s.equiripple_order(0.0, 0.01, 50.0, fs),
        lambda: s.equiripple_order(0.01, 0.01, 600.0, fs),
    ):
        with pytest.raises(ValueError):
            call()


def test_circular_convolution_aliases_and_zero_padding_removes_it():
    rng = np.random.default_rng(6)
    x, h = rng.normal(size=12), rng.normal(size=5)
    linear = np.convolve(x, h)
    circular = s.circular_convolution(x, h)
    np.testing.assert_allclose(np.real(np.fft.ifft(np.fft.fft(x) * np.fft.fft(h, 12))), circular, atol=1e-12)
    aliased = linear[:12].copy()
    aliased[:4] += linear[12:]
    np.testing.assert_allclose(circular, aliased, atol=1e-12)
    np.testing.assert_allclose(s.circular_convolution(x, h, n=16), linear, atol=1e-12)
    np.testing.assert_allclose(s.circular_convolution([1j, 0.0], [1.0, 1.0]), [1j, 1j])
    for args in (([1.0, 2.0], [1.0], 1), (np.ones((2, 2)), [1.0], None), ([np.nan], [1.0], None)):
        with pytest.raises(ValueError):
            s.circular_convolution(*args)


def test_overlap_add_equals_linear_convolution():
    rng = np.random.default_rng(7)
    x, h = rng.normal(size=1000), rng.normal(size=31)
    for block in (None, 1, 64, 1000, 5000):
        y = s.overlap_add(x, h, block_length=block)
        assert np.isrealobj(y)
        np.testing.assert_allclose(y, np.convolve(x, h), atol=1e-10)
    z = x[:50] + 1j * x[50:100]
    np.testing.assert_allclose(s.overlap_add(z, h), np.convolve(z, h), atol=1e-10)
    with pytest.raises(ValueError):
        s.overlap_add(x, h, block_length=0)


def test_decimation_and_interpolation():
    fs = 1000.0
    n = np.arange(2000)
    slow, fast = np.sin(2 * np.pi * 20 * n / fs), np.sin(2 * np.pi * 400 * n / fs)
    np.testing.assert_array_equal(s.downsample(n, 3, phase=1), n[1::3])
    # Plain downsampling by 4 folds 400 Hz onto 100 Hz; decimating removes it.
    plain = s.downsample(slow + fast, 4)
    decimated = s.decimate(slow + fast, 4)
    reference = s.downsample(slow, 4)
    inner = slice(30, -30)
    assert np.abs(plain - reference)[inner].max() > 0.9
    assert np.abs(decimated - reference)[inner].max() < 0.01
    np.testing.assert_allclose(s.decimate(slow, 1), slow)
    up = s.upsample([1.0, 2.0], 2)
    np.testing.assert_array_equal(up, [1.0, 0.0, 2.0, 0.0])
    # Interpolation keeps the original samples and fills in the band-limited curve.
    y = s.interpolate(slow, 5, numtaps=81)
    np.testing.assert_allclose(y[::5], slow, atol=1e-12)
    fine = np.sin(2 * np.pi * 20 * np.arange(y.size) / (5 * fs))
    assert np.abs(y - fine)[200:-200].max() < 0.01
    np.testing.assert_allclose(s.interpolate(slow, 1), slow)
    for call in (
        lambda: s.downsample(n, 2, phase=2),
        lambda: s.downsample(np.ones((2, 2)), 2),
        lambda: s.decimate(slow, 0),
        lambda: s.decimate(slow, 2, numtaps=40),
        lambda: s.interpolate(slow, 2, numtaps=0),
    ):
        with pytest.raises(ValueError):
            call()


def test_hilbert_transform_analytic_signal_envelope_and_frequency():
    fs = 1000.0
    t = np.arange(1000) / fs
    x = np.cos(2 * np.pi * 50 * t)
    np.testing.assert_allclose(s.hilbert_transform(x), np.sin(2 * np.pi * 50 * t), atol=1e-12)
    np.testing.assert_allclose(s.analytic_signal(x), signal.hilbert(x), atol=1e-12)
    np.testing.assert_allclose(s.analytic_signal(x[:999]), signal.hilbert(x[:999]), atol=1e-12)
    # H(H x) = -x and x is orthogonal to H x.
    np.testing.assert_allclose(s.hilbert_transform(s.hilbert_transform(x)), -x, atol=1e-12)
    assert abs(x @ s.hilbert_transform(x)) < 1e-9
    am = (1 + 0.5 * np.cos(2 * np.pi * 5 * t)) * np.cos(2 * np.pi * 100 * t)
    np.testing.assert_allclose(s.envelope(am), 1 + 0.5 * np.cos(2 * np.pi * 5 * t), atol=1e-12)
    np.testing.assert_allclose(s.instantaneous_frequency(x, fs), 50.0, atol=1e-9)
    h = s.hilbert_fir(31)
    np.testing.assert_allclose(h, -h[::-1])
    f, H = s.frequency_response(h, 1, worN=np.array([0.0, 100.0, 250.0, 400.0]), fs=fs)
    assert abs(H[0]) < 1e-12 and np.allclose(np.abs(H[1:]), 1, atol=0.02)
    for call in (lambda: s.analytic_signal([1.0]), lambda: s.analytic_signal([1j, 1.0]), lambda: s.hilbert_fir(30)):
        with pytest.raises(ValueError):
            call()
    with pytest.raises(ValueError):
        s.instantaneous_frequency(x, 0.0)


def test_dct_is_orthonormal_matches_scipy_and_compacts_energy():
    rng = np.random.default_rng(8)
    x = rng.normal(size=32)
    np.testing.assert_allclose(s.dct(x), sfft.dct(x, norm="ortho"), atol=1e-12)
    np.testing.assert_allclose(s.idct(s.dct(x)), x, atol=1e-12)
    assert np.isclose(np.sum(s.dct(x) ** 2), np.sum(x**2))
    ramp = np.linspace(0, 1, 32)
    dct_share = np.sort(s.dct(ramp) ** 2)[::-1][:4].sum() / np.sum(ramp**2)
    dft_share = np.sort(np.abs(np.fft.fft(ramp)) ** 2 / 32)[::-1][:4].sum() / np.sum(ramp**2)
    assert dct_share > 0.999 > dft_share
    for bad in ([1j], np.ones((2, 2))):
        with pytest.raises(ValueError):
            s.dct(bad)


def _ar2(n, a, seed):
    e = np.random.default_rng(seed).normal(size=n)
    return signal.lfilter([1.0], a, e)


def test_levinson_durbin_solves_yule_walker_and_reports_reflections():
    r = np.array([2.0, 1.2, 0.3, -0.2])
    model = s.levinson_durbin(r, 3)
    a = model.coefficients
    R = np.array([[r[abs(i - k)] for k in range(3)] for i in range(3)])
    np.testing.assert_allclose(R @ a[1:], -r[1:])
    assert np.isclose(model.noise_variance, r @ a)
    np.testing.assert_allclose(model.prediction_errors[1:], model.prediction_errors[:-1] * (1 - model.reflection_coefficients**2))
    assert model.order == 3 and model.summary()["stable"]
    assert np.isclose(model.reflection_coefficients[0], -0.6)
    with pytest.raises(ValueError):
        s.levinson_durbin([1.0, 1.5], 1)
    for args in (([1.0, 0.5], 2), ([0.0, 0.5], 1), ([1j, 0.5], 1)):
        with pytest.raises(ValueError):
            s.levinson_durbin(*args)


def test_yule_walker_recovers_an_ar2_process_and_its_spectrum():
    a = np.array([1.0, -1.5, 0.9])
    x = _ar2(50000, a, 9)
    model = s.yule_walker(x, 2)
    np.testing.assert_allclose(model.coefficients, a, atol=0.01)
    assert np.isclose(model.noise_variance, 1.0, rtol=0.03)
    f, P = model.psd(fs=2.0, n=4001)
    # The model spectrum integrates to the record's variance (it matches r_0).
    assert np.isclose(np.trapezoid(P, f), np.var(x), rtol=1e-3)
    peak = f[np.argmax(P)]
    assert np.isclose(peak, np.arccos(1.5 * 1.9 / (4 * 0.9)) / np.pi, atol=0.005)
    summary = model.summary()
    assert summary["order"] == 2 and summary["prediction_gain"] > 5
    fig, ax = model.plot(fs=2.0)
    assert ax.get_ylabel() == "PSD [dB]"
    plt.close(fig)
    assert s.yule_walker(x + 5.0, 2, demean=False).coefficients[1] < -0.99
    with pytest.raises(ValueError):
        model.psd(n=1)
    with pytest.raises(ValueError):
        s.yule_walker(np.ones(10), 2)
    with pytest.raises(ValueError):
        s.yule_walker(x[:2], 2)


def test_real_cepstrum_finds_an_echo():
    rng = np.random.default_rng(10)
    pulse = signal.lfilter([1.0], [1.0, -0.6], rng.normal(size=64))
    x = np.concatenate([pulse, np.zeros(448)])
    echo = x + 0.5 * np.roll(x, 120)
    c = s.real_cepstrum(echo)
    assert np.argmax(c[20:256]) + 20 == 120
    assert np.isclose(c[120] - s.real_cepstrum(x)[120], 0.25, atol=0.02)
    c_padded = s.real_cepstrum(echo, n_fft=1024)
    assert c_padded.size == 1024 and np.allclose(c_padded[1:512], c_padded[1023:512:-1])
    with pytest.raises(ValueError):
        s.real_cepstrum([1.0, -1.0])
    with pytest.raises(ValueError):
        s.real_cepstrum(echo, n_fft=10)
