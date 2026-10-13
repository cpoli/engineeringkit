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
