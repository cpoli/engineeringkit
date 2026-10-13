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
