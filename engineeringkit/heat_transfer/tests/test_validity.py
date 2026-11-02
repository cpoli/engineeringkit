"""Each heat_transfer correlation warns outside the range its docstring states, and only there."""

import warnings

import numpy as np
import pytest

from engineeringkit import ValidityWarning
from engineeringkit import heat_transfer as h

WATER_FILM = (960.0, 0.6, 2.8e-4, 0.68, 2.257e6, 373.15)


@pytest.mark.parametrize(
    "call, match",
    [
        (lambda: h.dittus_boelter_nusselt(500.0, 0.7), r"dittus_boelter_nusselt is meant for Re >= 10000; got Re = 500"),
        (lambda: h.dittus_boelter_nusselt(1e5, 0.01), r"0.6 <= Pr <= 160; got Pr = 0.01"),
        (lambda: h.dittus_boelter_nusselt(np.array([1e5, 2e3]), 0.7, heating=False), r"got Re = 2000"),
        (lambda: h.churchill_chu_nusselt(1e11, 0.71, laminar=True), r"Ra <= 1e\+09; got Ra = 1e\+11"),
        (lambda: h.flat_plate_nusselt(1e5, 0.01), r"flat_plate_nusselt is meant for 0.6 <= Pr <= 60"),
        (lambda: h.flat_plate_nusselt(1e5, 100.0, local=True), r"got Pr = 100"),
        (lambda: h.churchill_bernstein_nusselt(0.1, 0.7), r"Re\*Pr >= 0.2; got Re\*Pr = 0.07"),
        (lambda: h.zukauskas_nusselt(5e4, 0.01, 0.03, 0.03), r"zukauskas_nusselt is meant for 0.7 <= Pr <= 500"),
        (lambda: h.zukauskas_nusselt(5e4, 0.7, 0.02, 0.03), r"ST/SL >= 0.7; got ST/SL = 0.666667"),
        (lambda: h.hausen_nusselt(5000.0, 5.0, 0.01, 1.0), r"hausen_nusselt is meant for Re <= 2300"),
        (lambda: h.nusselt_film_condensation(*WATER_FILM, 273.15, 20.0), r"nusselt_film_condensation is meant for Re_delta <= 1800"),
    ],
)
def test_correlation_warns_outside_its_range(call, match):
    with pytest.warns(ValidityWarning, match=match):
        result = call()
    assert np.all(np.isfinite(result))


def test_correlations_are_silent_inside_their_ranges():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        h.dittus_boelter_nusselt(np.array([1e4, 1e6]), np.array([0.6, 160.0]))
        h.churchill_chu_nusselt(1e12, 0.71)  # the full form holds at every Ra
        h.churchill_chu_nusselt(1e9, 0.71, laminar=True)
        h.flat_plate_nusselt(np.array([1e3, 1e7]), np.array([0.6, 60.0]))
        h.churchill_bernstein_nusselt(1.0, 0.2)
        h.zukauskas_nusselt(5e4, np.array([0.7, 500.0]), 0.02, 0.03, "staggered")  # ST/SL limit is for aligned banks
        h.zukauskas_nusselt(5e4, 0.7, 0.021, 0.03)
        h.hausen_nusselt(2300.0, 5.0, 0.01, 1.0)
        h.nusselt_film_condensation(*WATER_FILM, 363.15, 0.5)


def test_warning_points_at_the_caller():
    with pytest.warns(ValidityWarning) as record:
        h.dittus_boelter_nusselt(500.0, 0.7)
    assert record[0].filename == __file__


def test_out_of_range_result_is_unchanged():
    with pytest.warns(ValidityWarning):
        assert h.dittus_boelter_nusselt(500.0, 0.7) == pytest.approx(0.023 * 500.0**0.8 * 0.7**0.4, rel=1e-14)
