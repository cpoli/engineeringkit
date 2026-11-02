"""Worked examples from the signal-processing textbooks, checked to the precision the books print.

OS: A. V. Oppenheim and R. W. Schafer, *Discrete-Time Signal Processing*, 3rd ed., Pearson (2010).
PM: J. G. Proakis and D. G. Manolakis, *Digital Signal Processing: Principles, Algorithms, and
Applications*, 4th ed., Pearson Prentice Hall (2007).
"""

import numpy as np

from engineeringkit import signal_processing as s


def test_proakis_manolakis_example_1_4_3_aliased_frequencies():
    # 3cos(2000 pi t) + 5sin(6000 pi t) + 10cos(12000 pi t) sampled at 5 kHz: the
    # 3 kHz term aliases to 2 kHz and the 6 kHz term to 1 kHz.
    assert [s.alias_frequency(f, 5000.0) for f in (1000.0, 3000.0, 6000.0)] == [1000.0, 2000.0, 1000.0]


def test_oppenheim_schafer_example_4_2_aliasing_of_a_sinusoid():
    # cos(16000 pi t) sampled at 6 kHz gives the same samples as cos(4000 pi t).
    assert s.alias_frequency(8000.0, 6000.0) == 2000.0
    n = np.arange(12)
    assert np.allclose(np.cos(16000 * np.pi * n / 6000), np.cos(2 * np.pi * n / 3))


def test_proakis_manolakis_example_2_3_2_convolution():
    assert s.convolution([1, 2, 1, -1], [1, 2, 3, 1]).tolist() == [1, 4, 8, 8, 3, -2, -1]


def test_proakis_manolakis_example_2_6_1_crosscorrelation():
    r = s.correlation([2, -1, 3, 7, 1, 2, -3], [1, -1, 2, -2, 4, 1, -2, 5])
    assert r.tolist() == [10, -9, 19, 36, -14, 33, 0, 7, 13, -18, 16, -7, 5, -3]


def test_proakis_manolakis_example_7_2_1_circular_convolution():
    assert s.circular_convolution([2, 1, 2, 1], [1, 2, 3, 4]).tolist() == [14, 16, 14, 16]


def test_oppenheim_schafer_example_7_3_bilinear_butterworth():
    # 0.89125 <= |H| for |w| <= 0.2 pi, |H| <= 0.17783 for 0.3 pi <= |w|; Td = 1, so fs = 1.
    ripple, attenuation = -20 * np.log10(0.89125), -20 * np.log10(0.17783)
    assert s.iir_order(0.1, 0.15, ripple, attenuation, 1.0) == 6  # N = 5.305
    # Omega_c = 0.766 meets the stopband edge exactly and exceeds the passband.
    cutoff = np.arctan(0.766 / 2) / np.pi
    b, a = s.butterworth_lowpass(6, cutoff, 1.0)
    _, response = s.frequency_response(b, a, worN=[0.1, 0.15], fs=1.0)
    assert abs(response[0]) > 0.89125
    assert np.isclose(abs(response[1]), 0.17783, rtol=3e-3)


def test_oppenheim_schafer_example_7_4_design_comparisons():
    # 0.99 <= |H| <= 1.01 for |w| <= 0.4 pi, |H| <= 0.001 for 0.6 pi <= |w|.
    ripple, attenuation = 20 * np.log10(1.01 / 0.99), 20 * np.log10(1.01 / 0.001)
    orders = [s.iir_order(0.2, 0.3, ripple, attenuation, 1.0, family) for family in ("butterworth", "chebyshev", "elliptic")]
    assert orders == [14, 8, 6]


def test_oppenheim_schafer_example_7_8_kaiser_lowpass():
    # wp = 0.4 pi, ws = 0.6 pi, delta = 0.001: A = 60, beta = 5.653, M = 37.
    numtaps, beta = s.kaiser_order(60.0, 0.1, 1.0)
    assert round(beta, 3) == 5.653
    assert numtaps - 1 == 37


def test_oppenheim_schafer_section_7_8_1_parks_mcclellan_lowpass():
    # delta1 = 0.01, delta2 = 0.001, K = 10: Kaiser's estimate M = 25.34, so 26; M = 26
    # misses delta2 and M = 27 meets it.
    assert s.equiripple_order(0.01, 0.001, 0.1, 1.0) - 1 == 26
    stopband = np.linspace(0.3, 0.5, 2001)
    peaks = []
    for order in (26, 27):
        h = s.equiripple_lowpass(order + 1, 0.2, 0.3, 1.0, weight=10.0)
        _, response = s.frequency_response(h, [1.0], worN=stopband, fs=1.0)
        peaks.append(np.abs(response).max())
    assert peaks[0] > 0.001 > peaks[1]
