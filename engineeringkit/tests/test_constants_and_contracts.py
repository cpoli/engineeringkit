import warnings

import numpy as np
import pytest

import engineeringkit as ek
from engineeringkit import constants
from engineeringkit._validation import ValidityWarning, interval, validity


def test_si_constants():
    assert constants.STANDARD_GRAVITY == 9.80665
    assert constants.STANDARD_ATMOSPHERE == 101325
    assert np.isclose(constants.UNIVERSAL_GAS_CONSTANT, 8.31446261815324)


def test_ordered_interval_requires_real_values():
    with pytest.raises(ValueError, match="real"):
        interval("probability", 0.5 + 1j, 0, 1)


def test_validity_warns_with_the_stated_range_and_first_offender():
    with pytest.warns(ValidityWarning, match=r"^model is meant for Re >= 10000; got Re = 500$"):
        validity("model", "Re", np.array([2e4, 500.0, 100.0]), 1e4)
    with pytest.warns(ValidityWarning, match=r"meant for Ra <= 1e\+09; got Ra = 1e\+10"):
        validity("model", "Ra", 1e10, high=1e9)
    with pytest.warns(ValidityWarning, match=r"meant for 0.6 <= Pr <= 160; got Pr = 0.5"):
        validity("model", "Pr", 0.5, 0.6, 160)


def test_validity_is_silent_inside_the_range_and_at_its_ends():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        validity("model", "Pr", np.array([0.6, 7.0, 160.0]), 0.6, 160)
        validity("model", "Re", 1e4, 1e4)


def test_validity_warning_points_at_the_caller_and_can_be_silenced():
    def public(Re):
        validity("public", "Re", Re, 1e4)

    with pytest.warns(ValidityWarning) as record:
        public(500.0)
    assert record[0].filename == __file__
    assert issubclass(ValidityWarning, UserWarning)
    assert ek.ValidityWarning is ValidityWarning
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        warnings.filterwarnings("ignore", category=ValidityWarning)
        public(500.0)
