import numpy as np
import pytest

from engineeringkit.units import convert


def test_exact_length_conversion_and_roundtrip():
    assert convert(1, "in", "m") == 0.0254
    assert np.isclose(convert(1, "ft", "in"), 12)
    assert np.isclose(convert(convert(123, "mm", "ft"), "ft", "mm"), 123)


def test_pressure_scales_and_dimension_guard():
    assert convert(1, "bar", "Pa") == 1e5
    assert convert(1, "MPa", "kPa") == 1000
    assert np.isclose(convert(1, "psi", "Pa"), 6894.757293168)
    with pytest.raises(ValueError):
        convert(1, "m", "Pa")
    with pytest.raises(ValueError):
        convert(1, "C", "K")


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_nonfinite_conversion(value):
    with pytest.raises(ValueError):
        convert(value, "m", "mm")
