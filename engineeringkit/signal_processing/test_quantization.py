import numpy as np
import pytest

from engineeringkit.signal_processing import quantization_step, quantize, snr_db


def test_adc_has_exact_number_of_levels_and_saturates():
    levels = np.unique(quantize(np.linspace(-2, 2, 10001), -1, 1, 3))
    np.testing.assert_allclose(levels, np.linspace(-0.875, 0.875, 8))
    assert quantize(-100, -1, 1, 3) == -0.875
    assert quantize(100, -1, 1, 3) == 0.875


def test_in_range_error_bounded_by_half_step():
    x = np.linspace(-1, 1, 1001)
    assert np.max(abs(quantize(x, -1, 1, 8) - x)) <= quantization_step(-1, 1, 8) / 2 + 1e-15


@pytest.mark.parametrize("bits", [0, -1, 1.5, True, 53])
def test_invalid_bits(bits):
    with pytest.raises(ValueError):
        quantization_step(0, 1, bits)


def test_complex_snr():
    assert np.isclose(snr_db([1j, 1], [0.1j, 0.1]), 20)
