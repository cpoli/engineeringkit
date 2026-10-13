import numpy as np
import pytest

from engineeringkit import constants
from engineeringkit._validation import interval


def test_si_constants():
    assert constants.STANDARD_GRAVITY == 9.80665
    assert constants.STANDARD_ATMOSPHERE == 101325
    assert np.isclose(constants.UNIVERSAL_GAS_CONSTANT, 8.31446261815324)


def test_ordered_interval_requires_real_values():
    with pytest.raises(ValueError, match="real"):
        interval("probability", 0.5 + 1j, 0, 1)
