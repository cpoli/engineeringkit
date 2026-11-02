"""Each fluid_engineering model warns outside the range its docstring states, and only there."""

import warnings

import numpy as np
import pytest

from engineeringkit import ValidityWarning
from engineeringkit import fluid_engineering as f


@pytest.mark.parametrize(
    "call, match",
    [
        (lambda: f.stokes_settling_velocity(1e-3, 2650.0, 998.0, 1e-3), r"stokes_settling_velocity is meant for Re <= 1"),
        (lambda: f.stokes_settling_velocity(np.array([5e-5, 1e-3]), 50.0, 998.0, 1e-3), r"Re <= 1"),  # a rising sphere
        (lambda: f.rehbock_discharge_coefficient(6.0, 1.0), r"rehbock_discharge_coefficient is meant for H/P <= 5; got H/P = 6"),
        (lambda: f.pipe_wave_speed(2.2e9, 1000.0, 0.2, 0.02, 200e9), r"pipe_wave_speed is meant for D/e >= 25; got D/e = 10"),
    ],
)
def test_model_warns_outside_its_range(call, match):
    with pytest.warns(ValidityWarning, match=match):
        result = call()
    assert np.all(np.isfinite(result))


def test_models_are_silent_inside_their_ranges():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        f.stokes_settling_velocity(np.array([5e-5, 1e-5]), 2650.0, 998.0, 1e-3)
        f.stokes_settling_velocity(1e-3, 2650.0, 0.0, 1e-3)  # in vacuum Re is zero
        f.rehbock_discharge_coefficient(np.array([0.0, 5.0]), 1.0)
        f.pipe_wave_speed(2.2e9, 1000.0, 0.5, 0.02, 200e9)
        f.colebrook_friction_factor(1000.0, 1e-4)  # the laminar branch is 64/Re by design


def test_out_of_range_result_is_unchanged():
    with pytest.warns(ValidityWarning):
        assert f.rehbock_discharge_coefficient(6.0, 1.0) == pytest.approx(0.611 + 0.075 * 6.0)
